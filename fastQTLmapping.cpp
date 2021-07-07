#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include <numeric>
#include <cmath>
#include <ctime>
#include <initializer_list>
#include <algorithm>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include "mkl.h"
#include <omp.h>
#include <iomanip>
#include "fastQTLmapping.h"

using namespace std;
using namespace meqtllib;

namespace {
const uint64_t kTransUG = 0x0303030303030303ull;  // 0303 = 0000001100000011
const uint64_t kTransSG = 0x0003000300030003ull;  // 0003 = 0000000000000011
uint64_t ConvertG64(const std::array<uint64_t, 4> &g64) {
  uint64_t geno64;
  geno64 = g64[3] & kTransUG;
  geno64 <<= 2;
  geno64 |= g64[2] & kTransUG;
  geno64 <<= 2;
  geno64 |= g64[1] & kTransUG;
  geno64 <<= 2;
  geno64 |= g64[0] & kTransUG;
  return geno64;
}
uint64_t ConvertG64(const std::array<uint64_t, 5> &g64) {
  uint64_t geno64;
  geno64 = g64[4] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[3] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[2] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[1] & kTransSG;
  geno64 <<= 3;
  geno64 |= g64[0] & kTransSG;
  return geno64;
}
}  // namespace

namespace snplib {
void SNP::ConvertGeno(size_t num_snps, size_t idx,
                      std::array<uint64_t, 4> &g64) {
  uint64_t g8[32] = {
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
      0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull};
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_g = geno_ + i * num_bytes_;
    g8[i] = static_cast<uint64_t>(tmp_g[idx]);
  }
  for (size_t i = 0; i < 4; ++i) {
    g64[i] = g8[i] + (g8[4 + i] << 8) + (g8[8 + i] << 16) + (g8[12 + i] << 24) +
             (g8[16 + i] << 32) + (g8[20 + i] << 40) + (g8[24 + i] << 48) +
             (g8[28 + i] << 56);
  }
}
void SNP::ConvertGeno(size_t num_snps, size_t idx,
                      std::array<uint64_t, 5> &g64) {
  uint64_t g8[20] = {0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
                     0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
                     0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull,
                     0x55ull, 0x55ull, 0x55ull, 0x55ull, 0x55ull};
  for (size_t i = 0; i < num_snps; ++i) {
    auto *tmp_g = geno_ + i * num_bytes_;
    g8[i] = static_cast<uint64_t>(tmp_g[idx]);
  }
  for (size_t i = 0; i < 5; ++i) {
    g64[i] =
        g8[i] + (g8[5 + i] << 16) + (g8[10 + i] << 32) + (g8[15 + i] << 48);
  }
}
SNP::SNP(const uint8_t *geno, size_t num_samples)
    : geno_(geno),
      num_samples_(num_samples),
      num_full_bytes_(num_samples / 4),
      num_samples_left_(num_samples_ % 4),
      num_bytes_(num_full_bytes_ + (num_samples_left_ != 0 ? 1 : 0)) {}
void SNP::TransposeGeno(size_t num_snps, size_t idx, uint64_t *geno64) {
  std::array<uint64_t, 4> g64;
  for (size_t i = 0; i < num_full_bytes_; ++i) {
    ConvertGeno(num_snps, i, g64);
    geno64[32 * (4 * i) + idx] = ConvertG64(g64);
    for (size_t j = 1; j < 4; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[32 * (4 * i + j) + idx] = ConvertG64(g64);
    }
  }
  if (num_samples_left_ > 0) {
    ConvertGeno(num_snps, num_full_bytes_, g64);
    geno64[32 * (4 * num_full_bytes_) + idx] = ConvertG64(g64);
    for (size_t j = 1; j < num_samples_left_; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[32 * (4 * num_full_bytes_ + j) + idx] = ConvertG64(g64);
    }
  }
}
void SNP::TransposeGeno(size_t num_snps, uint64_t *geno64) {
  std::array<uint64_t, 5> g64;
  for (size_t i = 0; i < num_full_bytes_; ++i) {
    ConvertGeno(num_snps, i, g64);
    geno64[4 * i] = ConvertG64(g64);
    for (size_t j = 1; j < 4; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[4 * i + j] = ConvertG64(g64);
    }
  }
  if (num_samples_left_ > 0) {
    ConvertGeno(num_snps, num_full_bytes_, g64);
    geno64[4 * num_full_bytes_] = ConvertG64(g64);
    for (size_t j = 1; j < num_samples_left_; ++j) {
      for (auto &iter : g64) {
        iter >>= 2;
      }
      geno64[4 * num_full_bytes_ + j] = ConvertG64(g64);
    }
  }
}

void UnpackGeno(char *bfileName, vector<vector<double> >& geno_data, vector<vector<int> >& NA_data, size_t num_samples, size_t num_snps) {
  std::ifstream inputFile;
  double *snp_mask_d = new double[num_samples];

  inputFile.open(bfileName, std::ios::in | std::ios::binary);
  std::istreambuf_iterator<char> inputFile_it(inputFile);
  inputFile_it++; inputFile_it++; inputFile_it++;
  std::vector<uint8_t> geno( (inputFile_it),
                          std::istreambuf_iterator<char>() );

  const std::array<double, 4> geno_table{2.0, 0.0, 1.0, 0.0}; // Homozygote A1, missing, Heterozygote, Homozygote A2
  const std::array<double, 4> mask_table{0.0, 1.0, 0.0, 0.0};
  SNP snp(geno.data(), num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    double *snp_geno_d = geno_data[i].data();    
    snp.UnpackGeno(geno_table, mask_table, snp_geno_d, snp_mask_d);
    for (size_t j = 0; j < num_samples; ++j) {
        if (snp_mask_d[j] > 0) {
            NA_data[i].push_back(j);
        }
    }
    snp += 1;
  }
  inputFile.close();
}
}  // namespace snplib


namespace meqtllib {
void calcBfileSize(char *bfileNameRoot, int &num_samples, int &num_snps)
{
    char bfileName[strlen(bfileNameRoot)+10];
    std::string s;
    std::ifstream inputFile;

    strcpy(bfileName, bfileNameRoot);
    strcat(bfileName, ".fam");
    inputFile.open(bfileName, std::ios::in);
    for (num_samples = 0; std::getline(inputFile, s); ++num_samples)
    ;
    inputFile.close();

    strcpy(bfileName, bfileNameRoot);
    strcat(bfileName, ".bim");
    inputFile.open(bfileName, std::ios::in);
    for (num_snps = 0; std::getline(inputFile, s); ++num_snps)
    ;
    inputFile.close();
}

void calcInputSize(char* omicsFileName, int& omicsNum, int& sampleSize)
{
    int i, rowsCount, columnsCount;
    ifstream inputFile;
    string s;
    istringstream is;
    string one_item;

    inputFile.open(omicsFileName);
    assert(inputFile.is_open());
    getline(inputFile, s);
    is.str(s);
    columnsCount = 0;
    while (is >> one_item){
        columnsCount++;
    }
    rowsCount = 1;
    while (getline(inputFile, s)){
        rowsCount++;
    }
    inputFile.close();
    omicsNum = columnsCount;
    sampleSize = rowsCount;
    s.clear(); is.str(""); is.clear();
}

//input 2 dimension omics data
void input2Ddouble(vector<vector<double> >& omicsData, char* fileName, vector<vector<int> >& NASignMark, string NASign, int omicsNum, int sampleSize)
{
    int i, j;
    ifstream inputFile;
    string s;
    istringstream is;
    string one_item;

    inputFile.open(fileName);
    for (i = 0; i < sampleSize; i++){
        getline(inputFile, s);
        is.str(s);
        for (j = 0; j < omicsNum; j++){
            is >> one_item;
            if (one_item == NASign){
                omicsData[j][i] = 0.0;
                NASignMark[j].push_back(i); // transpose input matrix
            } else {
                omicsData[j][i] = stod(one_item); // transpose input matrix
            }
        }
        s.clear(); is.str(""); is.clear();
    }
    inputFile.close();
}

int matrix_inv(vector<double>& a)
{
    
    int status;
    MKL_INT ipiv[VARNUM];
    MKL_INT info;
    MKL_INT varNum = VARNUM, inc = 1;
    
    info = LAPACKE_dgetrf (LAPACK_ROW_MAJOR, varNum, varNum, a.data(), varNum, ipiv);
    if (info != 0){
        return info;
    }
    info = LAPACKE_dgetri (LAPACK_ROW_MAJOR, varNum, a.data(), varNum, ipiv);
    if (info != 0){
        return info;
    }
    return 0;
}

linearFitRlt linearFit(int currentOmics1, int currentOmics2, 
                     vector<vector<double> >& omics1Data, vector<vector<double> >& omics2Data,
                     vector<vector<int> >& NASignMark1, vector<vector<int> >& NASignMark2,
                     vector<double>& omics1RowSum, vector<double>& omics1SqrRowSum, vector<double>& omics2RowSum, 
                     vector<double>& omics1VectorWithNA, vector<double>& omics2VectorWithNA,
                     vector<double>& b_hat_tmp, vector<double>& b_hat, vector<double>& Y_hat, vector<double>& e,
                     vector<double>& c)
{
    int i,j,k;
    int n = sampleSize;
    int m = covarNum;
    vector<int> NASignMarkTmp;
    int status = 0;
    MKL_INT varNum = VARNUM, inc = 1;

    gsl_set_error_handler_off(); // GSL error handler off

    linearFitRlt rlt;
    rlt.currentOmics1=currentOmics1 + 1; rlt.currentOmics2=currentOmics2 + 1;

    // omit NA
    for (auto l : NASignMark1[currentOmics1]){
        NASignMarkTmp.push_back(l);
    }
    for (auto l : NASignMark2[currentOmics2]){
        NASignMarkTmp.push_back(l);
    }
    sort(NASignMarkTmp.begin(), NASignMarkTmp.end());
    NASignMarkTmp.erase(unique(NASignMarkTmp.begin(), NASignMarkTmp.end()), NASignMarkTmp.end());
    // QC by missing-rate threshold
    if (NASignMarkTmp.size() >= missingRateThd * n) {
      rlt.status = 1;
      return rlt;
    }
    n = n - NASignMarkTmp.size();
    int df = n - (m + 1) - 1;
    bool NAflag = NASignMarkTmp.size();

// extract matrix X
    double omics1RowSumWithNA, omics1SqrRowSumWithNA, omics2RowSumWithNA;
    if (!NAflag){
    } else {
        omics1RowSumWithNA = omics1RowSum[currentOmics1]; omics1SqrRowSumWithNA = omics1SqrRowSum[currentOmics1];
        omics2RowSumWithNA = omics2RowSum[currentOmics2];
        k = 0; // subscript of NASignMark
        for (i = 0; i < n+NASignMarkTmp.size(); i++){ // be careful, i contains NA cases
            if (k < NASignMarkTmp.size()){
                if (i == NASignMarkTmp[k]){
                    omics1RowSumWithNA -= omics1Data[currentOmics1][i]; omics1SqrRowSumWithNA -= omics1Data[currentOmics1][i] * omics1Data[currentOmics1][i];
                    omics2RowSumWithNA -= omics2Data[currentOmics2][i];
                    k++;
                    continue;
                }
            }
            omics1VectorWithNA[i-k] = omics1Data[currentOmics1][i];
            omics2VectorWithNA[i-k] = omics2Data[currentOmics2][i];
        }
        for (i = n; i < n+NASignMarkTmp.size(); i++){
            omics1VectorWithNA[i] = 0;
            omics2VectorWithNA[i] = 0;
        }

        // QC by MAF threshold
        if (omics1RowSumWithNA * 0.5 / n < MAFThd) {
          rlt.status = 4;
          return rlt;
        }
    }

// error handle constant valconstante
    bool constFlagOmics1 = true, constFlagOmics2 = true;
    if (!NAflag){
        for (i = 1; i < n; i++){
            if (omics1Data[currentOmics1][i] != omics1Data[currentOmics1][i - 1]){
                constFlagOmics1 = false; break;
            }
        }
        for (i = 1; i < n; i++){
            if (omics2Data[currentOmics2][i] != omics2Data[currentOmics2][i - 1]){
                constFlagOmics2 = false; break;
            }
        }
    } else {
        for (i = 1; i < n; i++){
            if (omics1VectorWithNA[i] != omics1VectorWithNA[i-1]){
                constFlagOmics1 = false; break;
            }
        }
        for (i = 1; i < n; i++){
            if (omics2VectorWithNA[i] != omics2VectorWithNA[i-1]){
                constFlagOmics2 = false; break;
            }
        }
    }
    if (constFlagOmics1 || constFlagOmics2){
        rlt.status = 2;
        return rlt;
    }

// calc beta
    if (!NAflag){
        c[0] = n;
        c[1] = omics1RowSum[currentOmics1]; c[2] = omics1RowSum[currentOmics1];
        c[3] = omics1SqrRowSum[currentOmics1];
        status = matrix_inv(c);
// error handle matrix inversion error
        if (status){
            rlt.p = -9; rlt.status = 3;
            return rlt;
        }
        b_hat_tmp[0] = omics2RowSum[currentOmics2];
        b_hat_tmp[1] = cblas_ddot ((MKL_INT) sampleSize, omics1Data[currentOmics1].data(), 1, omics2Data[currentOmics2].data(), 1);
        cblas_dgemv (CblasRowMajor, CblasNoTrans, varNum, varNum, 1, c.data(), varNum, b_hat_tmp.data(), inc, 0, b_hat.data(), inc);
    } else {
        c[0] = n;
        c[1] = omics1RowSumWithNA; c[2] = omics1RowSumWithNA;
        c[3] = omics1SqrRowSumWithNA;
        status = matrix_inv(c); 
    // error handle matrix inversion error
        if (status){
            rlt.p = -9; rlt.status = 3;
            return rlt;
        }
        b_hat_tmp[0] = omics2RowSumWithNA;
        b_hat_tmp[1] = cblas_ddot ((MKL_INT) sampleSize, omics1VectorWithNA.data(), 1, omics2VectorWithNA.data(), 1);
        cblas_dgemv (CblasRowMajor, CblasNoTrans, varNum, varNum, 1, c.data(), varNum, b_hat_tmp.data(), inc, 0, b_hat.data(), inc);
    }

// test beta
    double t,p_b;
    double MSE;
    if (!NAflag){
        e.assign(sampleSize, b_hat[0]);
        cblas_daxpy ((MKL_INT) sampleSize, b_hat[1], omics1Data[currentOmics1].data(), inc, e.data(), inc);
        cblas_daxpy ((MKL_INT) sampleSize, -1, omics2Data[currentOmics2].data(),inc, e.data(), inc);
    } else {
        e.assign(sampleSize, b_hat[0]);
        cblas_daxpy ((MKL_INT) sampleSize, b_hat[1], omics1VectorWithNA.data(),inc, e.data(), inc);
        cblas_daxpy ((MKL_INT) sampleSize, -1, omics2VectorWithNA.data(), inc, e.data(), inc);
    }
    MSE = cblas_ddot ((MKL_INT) sampleSize, e.data(), 1, e.data(), 1);
    MSE /= df;

    t = b_hat[1]/sqrt(MSE*c[3]);
    p_b = gsl_cdf_tdist_Q(abs(t), df)*2;

    // handle error of out of range
    if (p_b < 1e-308){
        p_b = 1e-308;
    }
    
    rlt.b = (float)b_hat[1]; rlt.t = (float)t;  rlt.p = (float)-log10(p_b);
    rlt.se = (float)b_hat[1] / t; 
    float t2 = pow(t,2); rlt.r2 = (float)t2 / (df + t2);
    rlt.df = df;
    rlt.status = 0;

    return rlt;
}
}  // namespace meqtllib

int main(int argc, char **argv){
    char *omics1FileName, *omics2FileName, *outputFileName;
    int i, j, thread_count;
    string NASign;
    string model;
    ifstream inputFile;
    string one_item;
    long precision_config;
    float pFilter1Level, pFilter2Level;

    omics1FileName = argv[1];
    omics2FileName = argv[2];
    outputFileName = argv[3];
    NASign = argv[4];
    missingRateThd = atof(argv[5]);
    MAFThd = atof(argv[6]);
    pFilter1Level = -log10(atof(argv[7]));
    pFilter2Level = -log10(atof(argv[8]));
    thread_count = stol(argv[9]);

    // calculate bfile size
    calcBfileSize(omics1FileName, sampleSize, omics1Num);

    vector<vector<double> > omics1Data(omics1Num, vector<double> (sampleSize)); // SNPnum * samples
    vector<vector<int> > NASignMark1(omics1Num, vector<int>(0)); // NA mark for second omics layer, SNPnum * NAs

    // input bfile
    char bfileName[strlen(omics1FileName)+10];
    strcpy(bfileName, omics1FileName); strcat(bfileName, ".bed");
    snplib::UnpackGeno(bfileName, omics1Data, NASignMark1, sampleSize, omics1Num);

    // calculate input file size
    calcInputSize(omics2FileName, omics2Num, sampleSize);

    vector<vector<double> > omics2Data(omics2Num, vector<double> (sampleSize)); // Traitnum * samples
    vector<vector<int> > NASignMark2(omics2Num, vector<int>(0)); // NA mark for second omics layer, Traitnum * NAs

    //input omics2ylation data
    input2Ddouble(omics2Data, omics2FileName, NASignMark2, NASign, omics2Num, sampleSize);

    char outputLogFileName[strlen(outputFileName)+10];
    strcpy(outputLogFileName, outputFileName);
    strcat(outputLogFileName, ".log");

    ofstream outputLogFile;
    outputLogFile.open(outputLogFileName);
    outputLogFile << omics1FileName << endl;
    outputLogFile << omics2FileName << endl;
    outputLogFile << "omics 1 number : " << omics1Num << endl;
    outputLogFile << "omics 2 number : " << omics2Num << endl;
    outputLogFile << "sample number : " << sampleSize << endl;
    outputLogFile << endl;

    // preprocess the row sum for matrix c
    vector<double> omics1RowSum(omics1Num);
    vector<double> omics2RowSum(omics2Num);
    vector<double> omics1SqrRowSum(omics1Num);
    double omicsRowSumCurrent, omicsSqrRowSumCurrent;
    for (i = 0; i < omics1Num; i++){
        omicsRowSumCurrent = 0; omicsSqrRowSumCurrent = 0;
        for (j = 0; j < sampleSize; j++){
            omicsRowSumCurrent += omics1Data[i][j];
            omicsSqrRowSumCurrent += omics1Data[i][j] * omics1Data[i][j];
        }
        omics1RowSum[i] = omicsRowSumCurrent;
        omics1SqrRowSum[i] = omicsSqrRowSumCurrent;
    }
    for (i = 0; i < omics2Num; i++){
        omicsRowSumCurrent = 0; 
        for (j = 0; j < sampleSize; j++){
            omicsRowSumCurrent += omics2Data[i][j];
        }
        omics2RowSum[i] = omicsRowSumCurrent;
    }

    // main process
    double time_start_whole = omp_get_wtime(), time_end_whole;
    double time_start, time_end;
    int progressFlag = max((int)omics1Num / 10,1);
    
    // current result for each thread
    linearFitRlt rltArr[omics2Num];

    // header of result file for P < 1e-10
    char output1LevelFileName[strlen(outputFileName)+10];
    strcpy(output1LevelFileName, outputFileName);
    strcat(output1LevelFileName, ".strict");
    ofstream output1LevelFile;
    output1LevelFile.open(output1LevelFileName);
    output1LevelFile << "SNP.id\t" << "Trait.id\t" << "BETA\t" << "SE\t" << "R2\t" << "T\t" << "-lgP\t" << "DF\n";
    output1LevelFile.close();

    // header of result file for P < 1e-2
    char output2LevelFileName[strlen(outputFileName)+10];
    strcpy(output2LevelFileName, outputFileName);
    strcat(output2LevelFileName, ".loose");
    ofstream output2LevelFile;
    output2LevelFile.open(output2LevelFileName);
    output2LevelFile << "SNP.id\t" << "Trait.id\t" << "BETA\t" << "T\t" << "-lgP\t" << "DF\n";
    output2LevelFile.close();

#       pragma omp parallel for num_threads(thread_count) schedule(dynamic) \
        private(i, j, time_start, time_end, rltArr) \
        shared(omics1Data, omics2Data,  \
        omics1Num, omics2Num, covarNum, sampleSize, \
        omics1RowSum, omics2RowSum, omics1SqrRowSum, \
        NASignMark1, NASignMark2, model, \
        progressFlag, \
        precision_config, output1LevelFileName, output2LevelFileName)
        for (i = 0; i < omics1Num; i++){
            // time_start=omp_get_wtime();
            if (i % progressFlag == 0){
                time_end=omp_get_wtime();
                cout << (float)i / omics1Num * 100 << "% finish, time use:" << time_end - time_start_whole << "s" << endl;
            }

            // number of items in current result
            int sigNum = 0;

            // alloc tmporary matrixs and vectors at once
            vector<double> omics1VectorWithNA(sampleSize), omics2VectorWithNA(sampleSize), b_hat_tmp(2), b_hat(2), Y_hat(sampleSize), e(sampleSize), c(2 * 2);
            for (j = 0; j < omics2Num; j++){
                linearFitRlt rlt;
                rlt = linearFit(i, j, 
                      omics1Data, omics2Data, 
                      NASignMark1, NASignMark2, 
                      omics1RowSum, omics1SqrRowSum, omics2RowSum, 
                      omics1VectorWithNA, omics2VectorWithNA,
                      b_hat_tmp, b_hat, Y_hat, e, c);

                if (rlt.p >= pFilter2Level & rlt.status == 0){
                    rltArr[sigNum] = rlt;
                    sigNum++;
                }
            }

            // P< 1e-10, maker BETA SE R2 T P
            // 1e-10 < P < 1e-2, maker BETA T P
            // P > 1e-2, no output
#           pragma omp critical
            {
                ofstream output1LevelFile;
                output1LevelFile.open(output1LevelFileName, ios::app);
                output1LevelFile << setprecision(2);
                
                ofstream output2LevelFile;
                output2LevelFile.open(output2LevelFileName, ios::app);
                output2LevelFile << setprecision(2);

                for (j = 0; j < sigNum; j++){
                    if (rltArr[j].p >= pFilter1Level){
                        output1LevelFile << rltArr[j].currentOmics1 << "\t";
                        output1LevelFile << rltArr[j].currentOmics2 << "\t";
                        output1LevelFile << rltArr[j].b << "\t";
                        output1LevelFile << rltArr[j].se << "\t";
                        output1LevelFile << rltArr[j].r2 << "\t";
                        output1LevelFile << setiosflags(ios::fixed);
                        output1LevelFile << rltArr[j].t << "\t";
                        output1LevelFile << rltArr[j].p << "\t";
                        output1LevelFile << rltArr[j].df << "\n";
                        output1LevelFile << resetiosflags(ios::fixed);
                    } else {
                        output2LevelFile << rltArr[j].currentOmics1 << "\t";
                        output2LevelFile << rltArr[j].currentOmics2 << "\t";
                        output2LevelFile << rltArr[j].b << "\t";
                        output2LevelFile << setiosflags(ios::fixed);
                        output2LevelFile << rltArr[j].t << "\t";
                        output2LevelFile << rltArr[j].p << "\t";
                        output2LevelFile << rltArr[j].df << "\n";
                        output2LevelFile << resetiosflags(ios::fixed);
                    }
                }
                output1LevelFile.close();
                output2LevelFile.close();
            }

            // free temporary matrixs and vectors at once
            vector<double>().swap(omics1VectorWithNA); vector<double>().swap(omics2VectorWithNA);
            vector<double>().swap(b_hat); vector<double>().swap(b_hat_tmp);
            vector<double>().swap(Y_hat); vector<double>().swap(e);
            vector<double>().swap(c);

            // time_end = omp_get_wtime();
            // cout << i << "th SNP time use:" << time_end - time_start << "s" << endl;
            // time_start = omp_get_wtime();
        }

    // whole calc time stamp
    time_end_whole = omp_get_wtime();
    outputLogFile << "whole script time use:" << time_end_whole - time_start_whole << "s" << endl;
    outputLogFile << endl;
    outputLogFile.close();


    return 0;
}
