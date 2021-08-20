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
#include "mkl.h"
#include "omp.h"
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

void UnpackGeno(char *bfileName, float* geno_data, vector<vector<int> >& NA_data, size_t num_samples, size_t num_snps) {
  std::ifstream inputFile;
  float *snp_mask_d = new float[num_samples];

  inputFile.open(bfileName, std::ios::in | std::ios::binary);
  std::istreambuf_iterator<char> inputFile_it(inputFile);
  inputFile_it++; inputFile_it++; inputFile_it++;
  std::vector<uint8_t> geno( (inputFile_it),
                          std::istreambuf_iterator<char>() );

  const std::array<float, 4> geno_table{2.0, 0.0, 1.0, 0.0}; // Homozygote A1, missing, Heterozygote, Homozygote A2
  const std::array<float, 4> mask_table{0.0, 1.0, 0.0, 0.0};
  SNP snp(geno.data(), num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    float *snp_geno_d = &geno_data[i * num_samples];    
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

void getBfileSNPid(char *bfileNameRoot, int num_snps, 
                   vector<string>& omicsName, vector<int>& omicsCHR, vector<long>& omicsBP) {
    char bfileName[strlen(bfileNameRoot)+10];
    string s, one_item;
    ifstream inputFile;

    strcpy(bfileName, bfileNameRoot);
    strcat(bfileName, ".bim");
    inputFile.open(bfileName, ios::in);
    for (int i = 0; i < num_snps; i++) {
        getline(inputFile, s);
        istringstream is(s);
        
        is >> one_item; // CHR
        try {
            omicsCHR[i] = stoi(one_item);
        } catch (std::invalid_argument) {
            omicsCHR[i] = -1;
        }
        
        is >> one_item; // SNP
        omicsName[i] = one_item;
        
        is >> one_item; // CM
        
        is >> one_item; // BP
        try {
            omicsBP[i] = stol(one_item);
        } catch (std::invalid_argument) {
            omicsBP[i] = -1;
        }
    }
    inputFile.close();
}

void calcInputSize(char* omicsFileName, int& omicsNum, int& sampleSize) {
    int i, rowsCount, columnsCount;
    ifstream inputFile;
    string s, one_item;

    inputFile.open(omicsFileName);
    assert(inputFile.is_open());
    getline(inputFile, s);
    istringstream is(s);
    columnsCount = 0;
    while (is >> one_item){
        columnsCount++;
    }
    rowsCount = 1;
    while (getline(inputFile, s)){
        rowsCount++;
    }
    inputFile.close();
    sampleSize = columnsCount - 3; // first 3 cols are traits id
    omicsNum = rowsCount;
    s.clear();
}

//input 2 dimension omics data
void input2Dfloat(float* omicsData, char* fileName, vector<vector<int> >& NASignMark, char* NASign, 
                   vector<string>& omicsName, vector<int>& omicsCHR, vector<long>& omicsBP, 
                   int omicsNum, int sampleSize, 
                   string* dataArea, 
                   int thread_count) {
    ifstream inputFile;
    string one_line, one_item;
    string delimiter = " \t,";
    string::size_type pos,lastPos;

    inputFile.open(fileName);
    for (int i = 0; i < omicsNum; i++) {
        getline(inputFile, one_line);

        lastPos = one_line.find_first_not_of(delimiter, 0);
        pos = one_line.find_first_of(delimiter, lastPos);
        one_item = one_line.substr(lastPos, pos - lastPos); // SNP
        omicsName[i] = one_item;
        
        lastPos = one_line.find_first_not_of(delimiter, pos);
        pos = one_line.find_first_of(delimiter, lastPos);
        one_item = one_line.substr(lastPos, pos - lastPos); // CHR
        try {
            omicsCHR[i] = stoi(one_item);
        } catch (std::invalid_argument) {
            omicsCHR[i] = -1;
        }
        
        lastPos = one_line.find_first_not_of(delimiter, pos);
        pos = one_line.find_first_of(delimiter, lastPos);
        one_item = one_line.substr(lastPos, pos - lastPos); // BP
        try {
            omicsBP[i] = stol(one_item);
        } catch (std::invalid_argument) {
            omicsBP[i] = -1;
        }

        dataArea[i] = one_line.substr(pos, one_line.length() - pos);
    }
    inputFile.close();

    // input omics data
#   pragma omp parallel for schedule(dynamic) \
    num_threads(thread_count) \
    shared(omicsData, sampleSize, NASignMark, NASign, dataArea)
    for (int i = 0; i < omicsNum; i++) {
        int lineLength = dataArea[i].length();
        char s[lineLength]; strcpy(s, dataArea[i].c_str());
        int s_p = 0;

        for (int j = 0; j < sampleSize; j++) {
            while (s[s_p] == ' ' || s[s_p] == '\t' || s[s_p] == ',') s_p++; // skip space
            if (s_p >= lineLength) exit(1); // column number lack

            if (!strncmp(&s[s_p], NASign, strlen(NASign))) { // missing sign
                omicsData[i * sampleSize + j] = 0.0;  // transpose input matrix, and set NA value to 0.0
                NASignMark[i].push_back(j);
                s_p += strlen(NASign);
            } else if (s[s_p] == '-' || s[s_p] >= '0' && s[s_p] <= '9') {
                int sb = 1;
                if (s[s_p] == '-') {
                    sb = -1; s_p++;
                } // symbol

                double realData = 0;
                while (s[s_p] >= '0' && s[s_p] <= '9')
                    realData = realData * 10 + (s[s_p++] - '0'); // integer

                if (s[s_p] == '.') {
                    float k = 0.1;
                    s_p++;
                    while (s[s_p] >= '0' && s[s_p] <= '9') {
                        realData += (s[s_p++] - '0') * k;
                        k *= 0.1;
                    }
                } // decimal

                if (s[s_p] == 'e' || s[s_p] == 'E') {
                    s_p++;
                    int zd = 1;
                    if (s[s_p] == '+') {
                        zd = 1; s_p++;
                    } else if (s[s_p] == '-') {
                        zd = -1; s_p++;
                    } // zoom dirction

                    while (s[s_p] == '0') s_p++; // skip leading 0

                    int exponent = 0;
                    while (s[s_p] >= '0' && s[s_p] <= '9') {
                        exponent = exponent * 10 + (s[s_p++] - '0');
                    } // exponent
                    exponent = exponent * zd;

                    realData *= pow(10, exponent);
                } // scitific notation

                omicsData[i * sampleSize + j] = (float)sb * realData; // transpose input matrix
            } else {
                exit(2); // unknown symbol
            }
        }
        while (s[s_p] == ' ' || s[s_p] == '\t' || s[s_p] == ',') s_p++; if (s_p < lineLength) exit(1); // column number overflow
    }
}

// normalize each row in omics
void preprocessing(int omicsId, float* omicsData, float* omicsDataNorm, 
                   vector<double>& omicsRowSum, vector<double>& omicsRowSD, 
                   int sampleSize, vector<vector<int> >& NASignMarkCurr) {
    double omicsRowSumCurr, omicsRowSDCurr, omicsRowMeanCurr;
    int sampleSizeCurr;
    int j;

    // preprocess the row sum for omics
    omicsRowSumCurr = 0;
    for (j = 0; j < sampleSize; j++){
        omicsRowSumCurr += omicsData[omicsId * sampleSize + j];
    }
    omicsRowSum[omicsId] = omicsRowSumCurr;

    // preprocess the row SD for omics
    sampleSizeCurr = sampleSize - NASignMarkCurr[omicsId].size(); // only non-missing locus are involved in SD calculation
    omicsRowMeanCurr = omicsRowSum[omicsId] / sampleSizeCurr;
    omicsRowSDCurr = 0;
    for (j = 0; j < sampleSize; j++){
        if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), j) == NASignMarkCurr[omicsId].end()) {
            omicsRowSDCurr += pow(omicsData[omicsId * sampleSize + j] - omicsRowMeanCurr, 2);
        }
    }
    omicsRowSD[omicsId] = sqrt(omicsRowSDCurr / (sampleSizeCurr - 1));

    // normalize the row of omics
    sampleSizeCurr = sampleSize - NASignMarkCurr[omicsId].size(); // only non-missing locus are involved in SD calculation
    omicsRowMeanCurr = omicsRowSum[omicsId] / sampleSizeCurr;
    omicsRowSDCurr = omicsRowSD[omicsId];
    if (omicsRowSDCurr != 0) { // constant variants
        for (j = 0; j < sampleSize; j++){
            if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), j) == NASignMarkCurr[omicsId].end()) {
                omicsDataNorm[omicsId * sampleSize + j] = (omicsData[omicsId * sampleSize + j] - omicsRowMeanCurr) / omicsRowSDCurr;
            }
        }
    }
}

linearFitRlt linearFit(int currentOmics1, int currentOmics2, 
                     int sampleSize, float MAFThd, float missingRateThd,
                     float* omics1Data, float* omics2Data, 
                     float* omics1DataCurr, float* omics2DataCurr, 
                     vector<vector<int> >& NASignMark1, vector<vector<int> >& NASignMark2,
                     vector<double>& omics1RowSum, vector<double>& omics2RowSum) {
    int i;
    double corr;
    int sampleSizeCurr, df_t;
    vector<int> NASignMarkCurr;
    double A, B, C, D, E;
    double SXX, SYY, SXY;
    double beta0, beta1, MSE;
    double t, p_t;
    int status = 0;
    linearFitRlt rlt;
    rlt.currentOmics1=currentOmics1; rlt.currentOmics2=currentOmics2; // locus index start from 0

    // omit NA
    for (auto l : NASignMark1[currentOmics1]) {
        NASignMarkCurr.push_back(l);
    }
    for (auto l : NASignMark2[currentOmics2]) {
        NASignMarkCurr.push_back(l);
    }
    sort(NASignMarkCurr.begin(), NASignMarkCurr.end());
    NASignMarkCurr.erase(unique(NASignMarkCurr.begin(), NASignMarkCurr.end()), NASignMarkCurr.end());

    // degree of freedom
    sampleSizeCurr = sampleSize - NASignMarkCurr.size();
    df_t = sampleSizeCurr - 2;

    // QC by missing-rate threshold
    if (sampleSizeCurr <= (1 - missingRateThd) * sampleSize) {
      rlt.status = 1;
      return rlt;
    }

    // build current omics data with NA
    for (i = 0; i < sampleSize; i++) {
        omics1DataCurr[i] = omics1Data[i];
        omics2DataCurr[i] = omics2Data[i];
    }
    A = omics1RowSum[currentOmics1];
    C = omics2RowSum[currentOmics2];
    for (auto l : NASignMarkCurr) {
        omics1DataCurr[l] = 0;
        omics2DataCurr[l] = 0;
        A -= omics1Data[l];
        C -= omics2Data[l];
    }

    // QC by MAF threshold
    if (A * 0.5 < MAFThd * sampleSizeCurr) {
      rlt.status = 3;
      return rlt;
    }

    // calculate regression coefficients
    // A = sigma(X)
    // B = sigma(XX)
    // C = sigma(Y)
    // D = sigma(YY)
    // E = sigma(XY)
    B = (double)cblas_sdot ((MKL_INT) sampleSize, omics1DataCurr, 1, omics1DataCurr, 1);
    D = (double)cblas_sdot ((MKL_INT) sampleSize, omics2DataCurr, 1, omics2DataCurr, 1);
    E = (double)cblas_sdot ((MKL_INT) sampleSize, omics1DataCurr, 1, omics2DataCurr, 1);
    SXX = B - A * A / sampleSizeCurr;
    SYY = D - C * C / sampleSizeCurr;
    SXY = E - A * C / sampleSizeCurr;
    beta1 = SXY / SXX;
    MSE = (SYY - beta1 * SXY)/ df_t;
    t = beta1 * sqrt(SXX / MSE);

    // calculate pearson correlation
    corr = sqrt(pow(t, 2) / (df_t + pow(t, 2)));

    // test correlation significant
    rlt.t = (float)t;
    p_t = gsl_cdf_tdist_Q(abs(t), df_t) * 2; rlt.p = p_t;

    // handle error of out of range
    if (p_t < 1e-308){
        p_t = 1e-308;
    }
    // other paraments
    rlt.b = (float)beta1;
    rlt.se = (float)rlt.b / t; 
    rlt.r2 = (float)pow(corr, 2);
    rlt.nmiss = sampleSizeCurr;
    rlt.status = 0;

    return rlt;
}
}  // namespace meqtllib

int main(int argc, char **argv) {
    char *omics1FileName, *omics2FileName, *outputFileName;
    int thread_count;
    char* NASign;
    long precision_config;
    double cisP, transP;
    int omics1Num, omics2Num, sampleSize;
    long cisDist;
    float missingRateThd, MAFThd;

    omics1FileName = argv[1];
    omics2FileName = argv[2];
    outputFileName = argv[3];
    NASign = argv[4];
    missingRateThd = atof(argv[5]);
    MAFThd = atof(argv[6]);
    cisDist = stol(argv[7]);
    cisP = atof(argv[8]);
    transP = atof(argv[9]);
    thread_count = stoi(argv[10]);

    // global starting time stamp 
    double time_start_whole = omp_get_wtime(), time_end_whole;

    // initializing log file
    char outputLogFileName[strlen(outputFileName)+10];
    strcpy(outputLogFileName, outputFileName);
    strcat(outputLogFileName, ".log");
    ofstream outputLogFile;
    outputLogFile.open(outputLogFileName);

    // calculate bfile size
    calcBfileSize(omics1FileName, sampleSize, omics1Num);
    // calculate input file size
    calcInputSize(omics2FileName, omics2Num, sampleSize);
    // record data scale into log file
    outputLogFile << "omics 1 file : " << omics1FileName << endl;
    outputLogFile << "omics 2 file : " << omics2FileName << endl;
    outputLogFile << "omics 1 number : " << omics1Num << endl;
    outputLogFile << "omics 2 number : " << omics2Num << endl;
    outputLogFile << "sample number : " << sampleSize << endl;
    outputLogFile << endl;

    // critical value of t test
    double tCriticalValue, cisRCriticalValue, transRCriticalValue;
    tCriticalValue = gsl_cdf_tdist_Qinv(cisP/2, sampleSize - 2); // t critical value
    cisRCriticalValue = sqrt(pow(tCriticalValue, 2) / (sampleSize - 2 + pow(tCriticalValue, 2))) * (sampleSize - 1); // correlation critical value, times sampleSize
    outputLogFile << "cis-QTL pearson correlation critical value: " << cisRCriticalValue / (sampleSize - 1) << endl << endl;
    tCriticalValue = gsl_cdf_tdist_Qinv(transP/2, sampleSize - 2); // t critical value
    transRCriticalValue = sqrt(pow(tCriticalValue, 2) / (sampleSize - 2 + pow(tCriticalValue, 2))) * (sampleSize - 1); // correlation critical value, times sampleSize
    outputLogFile << "trans-QTL pearson correlation critical value: " << transRCriticalValue / (sampleSize - 1) << endl << endl;

    // header of result file for P < 1e-10
    char cisOutputFileName[strlen(outputFileName)+10];
    strcpy(cisOutputFileName, outputFileName);
    strcat(cisOutputFileName, ".cis.rlt");
    ofstream cisOutputFile;
    cisOutputFile.open(cisOutputFileName);
    cisOutputFile << "SNP.id\t" << "Trait.id\t" << "BETA\t" << "SE\t" << "R2\t" << "T\t" << "P\t" << "NMISS\n";
    cisOutputFile.close();
    // header of result file for P < 1e-2
    char transOutputFileName[strlen(outputFileName)+10];
    strcpy(transOutputFileName, outputFileName);
    strcat(transOutputFileName, ".trans.rlt");
    ofstream transOutputFile;
    transOutputFile.open(transOutputFileName);
    transOutputFile << "SNP.id\t" << "Trait.id\t" << "BETA\t" << "SE\t" << "R2\t" << "T\t" << "P\t" << "NMISS\n";
    transOutputFile.close();

    // input bfile
    float* omics1Data = (float*) mkl_malloc(sizeof(float) * omics1Num * sampleSize, 64);
    vector<vector<int> > NASignMark1(omics1Num, vector<int>(0)); // NA mark for second omics, SNPnum * NAs
    vector<string> omics1Name(omics1Num); // locus name for first omics
    vector<int> omics1CHR(omics1Num); // CHR number for first omics
    vector<long> omics1BP(omics1Num); // BP number for first omics
    getBfileSNPid(omics1FileName, omics1Num, omics1Name, omics1CHR, omics1BP);
    char bfileName[strlen(omics1FileName)+10];
    strcpy(bfileName, omics1FileName); strcat(bfileName, ".bed");
    snplib::UnpackGeno(bfileName, omics1Data, NASignMark1, sampleSize, omics1Num);
    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "SNP data input has finished, time used: " << time_end_whole - time_start_whole << " s" << endl << endl;

    // normalization of SNP data
    vector<double> omics1RowSum(omics1Num); vector<double> omics1RowSD(omics1Num);
    float* omics1DataNorm = (float*) mkl_malloc(sizeof(float) * omics1Num * sampleSize, 64);
#   pragma omp parallel for \
    num_threads(thread_count) \
    shared(omics1Num, \
           omics1Data, omics1DataNorm, \
           sampleSize, \
           omics1RowSum, omics1RowSD, NASignMark1)
    for (int i = 0; i < omics1Num; i++) {
        preprocessing(i, omics1Data, omics1DataNorm, omics1RowSum, omics1RowSD, sampleSize, NASignMark1);
    }
    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "SNP data normalization has finished, time used: " << time_end_whole - time_start_whole << " s" << endl << endl;  

    // input omics2ylation data
    float* omics2Data = (float*) mkl_malloc(sizeof(float) * omics2Num * sampleSize, 64);
    vector<vector<int> > NASignMark2(omics2Num, vector<int>(0)); // NA mark for second omics, Traitnum * NAs
    vector<string> omics2Name(omics2Num); // locus name for second omics
    vector<int> omics2CHR(omics2Num); // CHR number for second omics
    vector<long> omics2BP(omics2Num); // BP number for second omics
    string* dataArea = new string[omics2Num];
    input2Dfloat(omics2Data, omics2FileName, NASignMark2, NASign, omics2Name, omics2CHR, omics2BP, omics2Num, sampleSize, dataArea, thread_count);
    delete [] dataArea;
    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "trait data input has finished, time used: " << time_end_whole - time_start_whole << " s" << endl << endl;

    // normalization of trait data
    vector<double> omics2RowSum(omics2Num); vector<double> omics2RowSD(omics2Num);
    float* omics2DataNorm = (float*) mkl_malloc(sizeof(float) * omics2Num * sampleSize, 64);
#   pragma omp parallel for \
    num_threads(thread_count) \
    shared(omics2Num, \
           omics2Data, omics2DataNorm, \
           sampleSize, \
           omics2RowSum, omics2RowSD, NASignMark2)
    for (int i = 0; i < omics2Num; i++) {
        preprocessing(i, omics2Data, omics2DataNorm, omics2RowSum, omics2RowSD, sampleSize, NASignMark2);
    }
    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "trait data normalization has finished, time used: " << time_end_whole - time_start_whole << " s" << endl << endl;

    // split Block schedule
    int omics1BlockStride = min(omics1Num, blockSize);
    int omics2BlockStride = min(omics2Num, blockSize);
    int omics1BlockNum = (int)(omics1Num + omics1BlockStride - 1) / omics1BlockStride;
    int omics2BlockNum = (int)(omics2Num + omics2BlockStride - 1) / omics2BlockStride;
    int blockNum = omics1BlockNum * omics2BlockNum;

    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "preprocessing has finished, time used: " << time_end_whole - time_start_whole << " s" << endl << endl;

#   pragma omp parallel \
    num_threads(min(thread_count, blockNum)) \
    shared(omics1Data, omics2Data, \
           omics1Name, omics2Name, \
           omics1CHR, omics2CHR, omics1BP, omics2BP, \
           omics1DataNorm, omics2DataNorm, \
           omics1Num, omics2Num, sampleSize, MAFThd, missingRateThd, \
           omics1BlockStride, omics2BlockStride, \
           omics1RowSum, omics2RowSum, \
           omics1RowSD, omics2RowSD, \
           NASignMark1, NASignMark2, \
           cisDist, cisRCriticalValue, transRCriticalValue, \
           precision_config, cisOutputFileName, transOutputFileName, cisP, transP, \
           outputLogFile)
    {
        // mark of omics id in current block
        int omics1BlockHead, omics2BlockHead;

        // alloc space for omics data contains NA
        float* omics1DataCurr = (float*) mkl_malloc(sizeof(float) * sampleSize, 64);
        float* omics2DataCurr = (float*) mkl_malloc(sizeof(float) * sampleSize, 64);

        // alloc space for current result on each thread
        linearFitRlt* cisRltArr = new linearFitRlt[omics1BlockStride * omics2BlockStride];
        linearFitRlt* transRltArr = new linearFitRlt[omics1BlockStride * omics2BlockStride];
        float* corr = (float*) mkl_malloc(sizeof(float) * omics1BlockStride * omics2BlockStride, 64);

#       pragma omp for schedule(dynamic)
        for (int i = 0; i < blockNum; i++) {
            omics1BlockHead = (int) i % omics1BlockNum * omics1BlockStride;
            omics2BlockHead = (int) i / omics1BlockNum * omics2BlockStride;

            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                        min(omics1BlockStride, omics1Num - omics1BlockHead), min(omics2BlockStride, omics2Num - omics2BlockHead), sampleSize, 
                        1, &omics1DataNorm[omics1BlockHead * sampleSize], sampleSize, &omics2DataNorm[omics2BlockHead * sampleSize], sampleSize, 
                        0, corr, omics2BlockStride);

            // number of items in current result
            int cisSigNum = 0;
            int transSigNum = 0;
            int corrP = 0;
            for (int j = omics1BlockHead; j < min(omics1BlockHead + omics1BlockStride, omics1Num); j++) {
                for (int k = omics2BlockHead; k < min(omics2BlockHead + omics2BlockStride, omics2Num); k++) {
                    float corrCurr = abs(corr[corrP++]);
                    bool cisFilter = false, transFilter = false;
                    if (omics1CHR[j] == omics2CHR[k] && omics1CHR[j] > 0 && omics2CHR[k] > 0 &&
                        abs(omics1BP[j] - omics2BP[k]) <= cisDist && omics1BP[j] > 0 && omics2BP[k] > 0) {
                        if (corrCurr > cisRCriticalValue) {
                            cisFilter = true;
                        }
                    } else if (corrCurr > transRCriticalValue) {
                        transFilter = true;
                    }
  
                    if (cisFilter || transFilter) {
                        linearFitRlt rlt = linearFit(j, k, 
                              sampleSize, MAFThd, missingRateThd, 
                              &omics1Data[j * sampleSize], &omics2Data[k * sampleSize], 
                              omics1DataCurr, omics2DataCurr, 
                              NASignMark1, NASignMark2,
                              omics1RowSum, omics2RowSum);

                        if (cisFilter && rlt.p <= cisP && rlt.status == 0) {
                            cisRltArr[cisSigNum] = rlt;
                            cisSigNum++;
                        } else
                        if (transFilter && rlt.p <= transP && rlt.status == 0) {
                            transRltArr[transSigNum] = rlt;
                            transSigNum++;
                        }
                    }
                }
            }

            // P< P1, maker BETA SE R2 T P
            // P1 < P < P2, maker BETA T P
            // P > P2, do not output
#           pragma omp critical
            {
                ofstream cisOutputFile;
                cisOutputFile.open(cisOutputFileName, ios::app);
                cisOutputFile << setprecision(2);
                for (int j = 0; j < cisSigNum; j++){
                        cisOutputFile << omics1Name[cisRltArr[j].currentOmics1] << "\t";
                        cisOutputFile << omics2Name[cisRltArr[j].currentOmics2] << "\t";
                        cisOutputFile << cisRltArr[j].b << "\t";
                        cisOutputFile << cisRltArr[j].se << "\t";
                        cisOutputFile << cisRltArr[j].r2 << "\t";
                        cisOutputFile << cisRltArr[j].t << "\t";
                        cisOutputFile << cisRltArr[j].p << "\t";
                        cisOutputFile << cisRltArr[j].nmiss << "\n";
                }
                cisOutputFile.close();

                ofstream transOutputFile;
                transOutputFile.open(transOutputFileName, ios::app);
                transOutputFile << setprecision(2);
                for (int j = 0; j < transSigNum; j++){
                        transOutputFile << omics1Name[transRltArr[j].currentOmics1] << "\t";
                        transOutputFile << omics2Name[transRltArr[j].currentOmics2] << "\t";
                        transOutputFile << transRltArr[j].b << "\t";
                        transOutputFile << transRltArr[j].se << "\t";
                        transOutputFile << transRltArr[j].r2 << "\t";
                        transOutputFile << transRltArr[j].t << "\t";
                        transOutputFile << transRltArr[j].p << "\t";
                        transOutputFile << transRltArr[j].nmiss << "\n";
                }
                transOutputFile.close();
            }
            
            double time_end = omp_get_wtime();
            outputLogFile << i << " / " << blockNum << " block has finished, time used: " << time_end - time_start_whole << " s" << endl;
        }
        mkl_free(omics1DataCurr); mkl_free(omics2DataCurr);
    }

    // whole calc time stamp
    time_end_whole = omp_get_wtime();
    outputLogFile << "whole script time used: " << time_end_whole - time_start_whole << " s" << endl;
    outputLogFile << endl;
    outputLogFile.close();

    mkl_free(omics1Data); mkl_free(omics2Data);
    mkl_free(omics1DataNorm); mkl_free(omics2DataNorm);

    return 0;
}
