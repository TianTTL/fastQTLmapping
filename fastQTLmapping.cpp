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
void calcBfileSize(char *bfileNameRoot, int &num_samples, int &num_snps) {
    char bfileName[strlen(bfileNameRoot)+10];
    string s;
    ifstream inputFile;

    strcpy(bfileName, bfileNameRoot);
    strcat(bfileName, ".fam");
    inputFile.open(bfileName, ios::in);
    for (num_samples = 0; getline(inputFile, s); ++num_samples)
    ;
    inputFile.close();

    strcpy(bfileName, bfileNameRoot);
    strcat(bfileName, ".bim");
    inputFile.open(bfileName, ios::in);
    for (num_snps = 0; getline(inputFile, s); ++num_snps)
    ;
    inputFile.close();
}

void getBfileSNPid(char *bfileNameRoot, int num_snps, vector<string>& omicsName) {
    char bfileName[strlen(bfileNameRoot)+10];
    string s, one_item;
    istringstream is;
    ifstream inputFile;

    strcpy(bfileName, bfileNameRoot);
    strcat(bfileName, ".bim");
    inputFile.open(bfileName, ios::in);
    for (int i = 0; i < num_snps; i++) {
        is.str(s);
        is >> one_item; is >> one_item;
        omicsNum[i] = one_item;
        inputFile.close();
    }
}

void calcInputSize(char* omicsFileName, int& omicsNum, int& sampleSize) {
    int i, rowsCount, columnsCount;
    ifstream inputFile;
    string s, one_item;
    istringstream is;

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
void input2Ddouble(vector<vector<double> >& omicsData, char* fileName, vector<vector<int> >& NASignMark, string NASign, 
                   vector<string>& omicsName, int omicsNum, int sampleSize) {
    int i, j;
    ifstream inputFile;
    string s, one_item;
    istringstream is;

    inputFile.open(fileName);
    for (i = 0; i < sampleSize; i++){
        getline(inputFile, s);
        is.str(s);
        iss >> one_item;
        omicsName[i] = one_item;
        for (j = 0; j < omicsNum; j++){
            is >> one_item;
            if (one_item == NASign){
                omicsData[j][i] = 0.0;  // transpose input matrix, and set NA value to 0.0
                NASignMark[j].push_back(i);
            } else {
                omicsData[j][i] = stod(one_item); // transpose input matrix
            }
        }
        s.clear(); is.str(""); is.clear();
    }
    inputFile.close();
}

// normalize each row in omics
void preprocessing(vector<vector<double> >& omicsData, vector<double>& omicsRowSum, vector<double>& omicsRowSD, 
                   sampleSize, vector<vector<int> >& NASignMarkCurr) {
    double omicsRowSumCurr, omicsRowSDCurr, omicsRowMeanCurr;
    int sampleSizeCurr;
    // preprocess the row sum for omics
    for (i = 0; i < omicsNum; i++){
        omicsRowSumCurr = 0;
        for (j = 0; j < sampleSize; j++){
            omicsRowSumCurr += omicsData[i][j];
        }
        omicsRowSum[i] = omicsRowSumCurr;
    }
    // preprocess the row SD for omics
    for (i = 0; i < omicsNum; i++){
        sampleSizeCurr = sampleSize - NASignMarkCurr[i].size(); // only non-missing locus are involved in SD calculation
        omicsRowMeanCurr = omicsRowSum[i] / sampleSizeCurr;
        omicsRowSDCurr = 0;
        for (j = 0; j < sampleSize; j++){
            if (find(NASignMarkCurr[i].begin(), NASignMarkCurr[i].end(), j) == NASignMarkCurr[i].end()) {
                omicsRowSDCurr += pow(omicsData[i][j] - omicsRowMeanCurr, 2);
            }
        }
        omicsRowSD[i] = sqrt(omicsRowSumCurr / (sampleSizeCurr - 1));
    }
    // normalize the row of omics
    for (i = 0; i < omicsNum; i++){
        sampleSizeCurr = sampleSize - NASignMarkCurr[i].size(); // only non-missing locus are involved in SD calculation
        omicsRowMeanCurr = omicsRowSum[i] / sampleSizeCurr;
        omicsRowSDCurr = omicsRowSD[i];
        for (j = 0; j < sampleSize; j++){
            if (find(NASignMarkCurr[i].begin(), NASignMarkCurr[i].end(), j) == NASignMarkCurr[i].end()) {
                omicsData[i][j] = (omicsData[i][j] - omicsRowMeanCurr) / omicsRowSDCurr;
            }
        }
    }
}

linearFitRlt linearFit(int currentOmics1, int currentOmics2, 
                     vector<vector<double> >& omics1Data, vector<vector<double> >& omics2Data,
                     vector<vector<int> >& NASignMark1, vector<vector<int> >& NASignMark2,
                     vector<double>& omics1RowSum, vector<double>& omics1RowSD, vector<double>& omics2RowSD, 
                     double rCriticalValue)
{
    int i;
    int sampleSizeCurr, df;
    vector<int> NASignMarkCurr;
    int status = 0;
    MKL_INT varNum = VARNUM, inc = 1;
    linearFitRlt rlt;
    rlt.currentOmics1=currentOmics1; rlt.currentOmics2=currentOmics2; // locus index start from 0

    // calc pearson correlation
    double corr = cblas_ddot ((MKL_INT) sampleSize, omics1Data[currentOmics1].data(), 1, omics2Data[currentOmics2].data(), 1);
    if (corr < rCriticalValue) {
        rlt.status = 4;
        return rlt;
    }

    gsl_set_error_handler_off(); // GSL error handler off

    // omit NA
    for (auto l : NASignMark1[currentOmics1]){
        NASignMarkCurr.push_back(l);
    }
    for (auto l : NASignMark2[currentOmics2]){
        NASignMarkCurr.push_back(l);
    }
    sort(NASignMarkCurr.begin(), NASignMarkCurr.end());
    NASignMarkCurr.erase(unique(NASignMarkCurr.begin(), NASignMarkCurr.end()), NASignMarkCurr.end());

    // degree of freedom
    sampleSizeCurr = sampleSize - NASignMarkCurr.size();
    df = sampleSizeCurr - 2;

    // QC by missing-rate threshold
    if (sampleSizeCurr <= (1 - missingRateThd) * sampleSize) {
      rlt.status = 1;
      return rlt;
    }
    
    // calc MAF considering NA 
    double omics1RowSumWithNA = omics1RowSum[currentOmics1];
    for (i = 0; i < NASignMarkCurr.size(); i++){
        omics1RowSumWithNA -= omics1Data[currentOmics1][i];
    }
    // QC by MAF threshold
    if (omics1RowSumWithNA * 0.5 / sampleSizeCurr < MAFThd) {
      rlt.status = 3;
      return rlt;
    }

    // test correlation significant
    double t, p_t;
    t = sqrt(df) * corr / sqrt(1 - pow(corr, 2)); rlt.t = (float)t;
    p_t = gsl_cdf_tdist_Q(abs(t), df) * 2; rlt.p = (float)-log10(p_t);
    // handle error of out of range
    if (p_t < 1e-308){
        p_t = 1e-308;
    }
    // other paraments
    rlt.b = (float)corr * omics2RowSD[currentOmics2] / omics1RowSD[currentOmics1];
    rlt.se = (float)rlt.b / t; 
    rlt.r2 = (float)pow(r, 2);
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
    vector<vector<int> > NASignMark1(omics1Num, vector<int>(0)); // NA mark for second omics, SNPnum * NAs
    vector<string> omics1Name(omics1Num, string); // locus name for first omics
    // input bfile
    getBfileSNPid(omics1FileName, omics1Num, omics1Name);
    char bfileName[strlen(omics1FileName)+10];
    strcpy(bfileName, omics1FileName); strcat(bfileName, ".bed");
    snplib::UnpackGeno(bfileName, omics1Data, NASignMark1, sampleSize, omics1Num);

    // calculate input file size
    calcInputSize(omics2FileName, omics2Num, sampleSize);
    vector<vector<double> > omics2Data(omics2Num, vector<double> (sampleSize)); // Traitnum * samples
    vector<vector<int> > NASignMark2(omics2Num, vector<int>(0)); // NA mark for second omics, Traitnum * NAs
    vector<string> omics2Name(omics2Num, string); // locus name for second omics
    //input omics2ylation data
    input2Ddouble(omics2Data, omics2FileName, NASignMark2, NASign, omics2Name, omics2Num, sampleSize);

    // create log file
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

    // processing
    // normalization
    vector<double> omics1RowSum(omics1Num); vector<double> omics1RowSD(omics1Num);
    vector<double> omics2RowSum(omics2Num); vector<double> omics2RowSD(omics2Num);
    preprocessing(omics1Data, omics1RowSum, omics1RowSD, sampleSize, NASignMark1);
    preprocessing(omics2Data, omics2RowSum, omics2RowSD, sampleSize, NASignMark2);
    //  critical value of t test
    auto tCriticalValue = gsl_cdf_tdist_Qinv(pFilter2Level/2, sampleSize - 2);
    auto rCriticalValue = sqrt(pow(tCriticalValue, 2) / (sampleSize - 2 + pow(tCriticalValue, 2)));

    // alloc space for current result on each thread
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

    // main process
    double time_start_whole = omp_get_wtime(), time_end_whole;
    double time_end;
    int progressFlag = max((int)omics1Num / 10,1);

#       pragma omp parallel for num_threads(thread_count) schedule(dynamic) \
        private(i, j, time_end, rltArr) \
        shared(omics1Data, omics2Data,  \
        omics1Num, omics2Num, sampleSize, \
        omics1RowSum, omics2RowSum, omics1SqrRowSum, \
        NASignMark1, NASignMark2, model, \
        progressFlag, \
        precision_config, output1LevelFileName, output2LevelFileName)
        for (i = 0; i < omics1Num; i++){
            if (i % progressFlag == 0 && i != 0){
                time_end=omp_get_wtime();
                cout << (float)i / omics1Num * 100 << "% finish, time use:" << time_end - time_start_whole << "s" << endl;
            }

            // number of items in current result
            int sigNum = 0;

            // alloc tmporary matrixs and vectors at once
            for (j = 0; j < omics2Num; j++){
                linearFitRlt rlt;
                rlt = linearFit(i, j, 
                      omics1Data, omics2Data, 
                      NASignMark1, NASignMark2, 
                      omics1RowSum, omics1RowSD, omics2RowSD);

                if (rlt.p >= pFilter2Level & rlt.status == 0){
                    rltArr[sigNum] = rlt;
                    sigNum++;
                }
            }

            // P< P1, maker BETA SE R2 T P
            // P1 < P < P2, maker BETA T P
            // P > P2, do not output
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
                        output1LevelFile << omics1Name[rltArr[j].currentOmics1] << "\t";
                        output1LevelFile << omics2Name[rltArr[j].currentOmics2] << "\t";
                        output1LevelFile << rltArr[j].b << "\t";
                        output1LevelFile << rltArr[j].se << "\t";
                        output1LevelFile << rltArr[j].r2 << "\t";
                        output1LevelFile << setiosflags(ios::fixed);
                        output1LevelFile << rltArr[j].t << "\t";
                        output1LevelFile << rltArr[j].p << "\t";
                        output1LevelFile << rltArr[j].df << "\n";
                        output1LevelFile << resetiosflags(ios::fixed);
                    } else {
                        output2LevelFile << omics1Name[rltArr[j].currentOmics1] << "\t";
                        output2LevelFile << omics2Name[rltArr[j].currentOmics2] << "\t";
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
        }

    // whole calc time stamp
    time_end_whole = omp_get_wtime();
    outputLogFile << "whole script time use:" << time_end_whole - time_start_whole << "s" << endl;
    outputLogFile << endl;
    outputLogFile.close();

    return 0;
}
