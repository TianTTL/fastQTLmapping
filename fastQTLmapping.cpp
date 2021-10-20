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
#include "clipp.h"

using namespace meqtllib;
using namespace clipp;
using std::cout;
using std::vector;
using std::array;
using std::ifstream;
using std::ios;

namespace {
const uint64_t kTransUG = 0x0303030303030303ull;  // 0303 = 0000001100000011
const uint64_t kTransSG = 0x0003000300030003ull;  // 0003 = 0000000000000011
uint64_t ConvertG64(const array<uint64_t, 4> &g64) {
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
uint64_t ConvertG64(const array<uint64_t, 5> &g64) {
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
                      array<uint64_t, 4> &g64) {
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
                      array<uint64_t, 5> &g64) {
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
  array<uint64_t, 4> g64;
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
  array<uint64_t, 5> g64;
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

void UnpackGeno(string bfileNameRoot, float* geno_data, vector<vector<uint32_t> >& NA_data, size_t num_samples, size_t num_snps, 
                vector<bool>& sampleFltSign, uint32_t covarNANum, string teststr) {
  ifstream inputFile;
  float *snp_geno_d = new float[num_samples + covarNANum];
  float *snp_mask_d = new float[num_samples + covarNANum];
  inputFile.open(bfileNameRoot + ".bed", ios::in | ios::binary);

  std::istreambuf_iterator<char> inputFile_it(inputFile);
  inputFile_it++; inputFile_it++; inputFile_it++;
  std::vector<uint8_t> geno( (inputFile_it),
                          std::istreambuf_iterator<char>() );

  const array<float, 4> geno_table{2.0, 0.0, 1.0, 0.0}; // Homozygote A1, missing, Heterozygote, Homozygote A2
  const array<float, 4> mask_table{0.0, 1.0, 0.0, 0.0};
  SNP snp(geno.data(), num_samples);
  for (size_t i = 0; i < num_snps; ++i) {
    snp.UnpackGeno(geno_table, mask_table, snp_geno_d, snp_mask_d);
    uint64_t rowHeadPos = i * num_samples;
    uint32_t sid = 0;
    for (size_t j = 0; j < num_samples + covarNANum; ++j) {
        if (!sampleFltSign[j]) {
            geno_data[rowHeadPos + sid] = snp_geno_d[j];
            if (snp_mask_d[j] > 0) {
                NA_data[i].push_back(sid);
            }
            sid++;
        }
    }
    snp += 1;
  }
  inputFile.close();
}
}  // namespace snplib


namespace meqtllib {
void calcBfileSize(string bfileNameRoot, uint32_t &num_samples, uint64_t &num_snps) {
    string s;
    ifstream inputFile;

    inputFile.open(bfileNameRoot + ".fam", ios::in);
    assert(inputFile.is_open());
    for (num_samples = 0; getline(inputFile, s); ++num_samples)
    ;
    inputFile.close();

    inputFile.open(bfileNameRoot + ".bim", ios::in);
    assert(inputFile.is_open());
    for (num_snps = 0; getline(inputFile, s); ++num_snps)
    ;
    inputFile.close();
}

void getBfileSNPid(string bfileNameRoot, uint64_t num_snps, 
                   vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP) {
    string s, one_item;
    ifstream inputFile;

    inputFile.open(bfileNameRoot + ".bim", ios::in);
    for (uint32_t i = 0; i < num_snps; i++) {
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

void calcInputSize(string omicsFileName, uint32_t &sampleSize, uint64_t& omicsNum) {
    uint32_t i, rowsCount = 0, colsCount = 0;
    ifstream inputFile;
    string s, one_item;

    inputFile.open(omicsFileName);
    assert(inputFile.is_open());
    
    getline(inputFile, s);
    istringstream is(s);
    while (is >> one_item) {
        colsCount++;
    }
    sampleSize = colsCount - 3; // first 3 columns are loci info
    rowsCount++;

    while (getline(inputFile, s)) {
        rowsCount++;
    }
    inputFile.close();
    omicsNum = rowsCount;
    s.clear();
}

void calcCovarSize(string covarFileName, string NASign, uint64_t sampleSize, uint64_t& covarNum, 
                   vector<bool>& sampleFltSign, uint32_t& covarNANum) {
    uint32_t i, rowsCount = 0;
    ifstream inputFile;
    string s, one_item;

    inputFile.open(covarFileName);
    assert(inputFile.is_open());
    while (getline(inputFile, s)){
        rowsCount++;
        istringstream is(s);
        for (i = 0; i < sampleSize; i++) {
            is >> one_item;
            if (one_item == NASign){
                sampleFltSign[i] = true;
            }
        }
        s.clear(); is.str(""); is.clear(); 
    }
    inputFile.close();
    covarNum = rowsCount;
    s.clear();

    for (auto j : sampleFltSign) {
        if (j) {
            covarNANum += 1;
        }
    }
}

//input 2 dimension omics data
void input2DfloatParse(float* omicsData, string fileName, vector<vector<uint32_t> >& NASignMark, string NASign, 
                  vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP, 
                  uint64_t omicsNum, uint32_t sampleSize, 
                  string* dataArea, 
                  uint32_t threadMaxN, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum) {
    ifstream inputFile;
    string one_line, one_item;
    string delimiter = " \t";
    string::size_type pos,lastPos;
    inputFile.open(fileName);
    assert(inputFile.is_open());
    for (uint64_t i = 0; i < omicsNum; i++) {
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
    num_threads(threadMaxN) \
    shared(omicsData, sampleSize, NASignMark, NASign, dataArea, sampleFltSign, covarNANum)
    for (uint64_t i = 0; i < omicsNum; i++) {
        uint32_t lineLength = dataArea[i].length();
        char s[lineLength]; strcpy(s, dataArea[i].c_str());
        uint32_t s_p = 0;
        uint64_t rowHeadPos = i * sampleSize;
        uint32_t sid = 0;

        for (uint32_t j = 0; j < sampleSize + covarNANum; j++) {
            while (s[s_p] == ' ' || s[s_p] == '\t') s_p++; // skip space
            if (s_p >= lineLength) {
                cout << "column number lack\n";
                exit(1); // column number lack
            }

            if (!strncmp(&s[s_p], NASign.data(), NASign.length())) { // missing sign
                if (!sampleFltSign[j]) {
                    omicsData[rowHeadPos + sid] = 0.0;  // and set NA value to 0.0
                    NASignMark[i].push_back(sid);
                    sid++;    
                }
                s_p += NASign.length();
            } else if (s[s_p] == '-' || s[s_p] >= '0' && s[s_p] <= '9') {
                int32_t sb = 1;
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
                    int32_t zd = 1;
                    if (s[s_p] == '+') {
                        zd = 1; s_p++;
                    } else if (s[s_p] == '-') {
                        zd = -1; s_p++;
                    } // zoom dirction

                    while (s[s_p] == '0') s_p++; // skip leading 0

                    int32_t exponent = 0;
                    while (s[s_p] >= '0' && s[s_p] <= '9') {
                        exponent = exponent * 10 + (s[s_p++] - '0');
                    } // exponent
                    exponent = exponent * zd;

                    realData *= pow(10, exponent);
                } // scitific notation

                if (!sampleFltSign[j]) {
                    omicsData[rowHeadPos + sid] = (float)sb * realData;
                    sid++;
                }
            } else {
                cout << "unknown symbol\n";
                exit(2); // unknown symbol
            }
        }
        while (s[s_p] == ' ' || s[s_p] == '\t' || s[s_p] == '\n' || s[s_p] == '\r') s_p++; if (s_p < lineLength) {
            cout << "column number overflow\n";
            exit(1); // column number overflow
        }
    }
}

//input 2 dimension covariates data
void inputCovar(float* covarData, string fileName, 
                uint32_t covarNum, uint32_t sampleSize, vector<bool>& sampleFltSign, uint32_t covarNANum) {
    uint32_t i, j;
    ifstream inputFile;
    string s, one_item;

    inputFile.open(fileName);
    for (i = 0; i < covarNum; i++) {
        getline(inputFile, s);
        istringstream is(s);
        uint32_t sid = 0;
        for (j = 0; j < sampleSize + covarNANum; j++) {
            is >> one_item;
            if (!sampleFltSign[j]) {
                covarData[i * sampleSize + sid] = stod(one_item);
                sid++;
            }
        }
        s.clear(); is.str(""); is.clear();
    }
    inputFile.close();
}

// centralize one row of a matrix
void cntrl(float* a, uint64_t omicsId, 
           uint32_t sampleSize, 
           vector<double>& rowSD,
           vector<vector<uint32_t> >& NASignMarkCurr) {
    double rowSumTmp, rowSDTmp, rowMeanTmp;
    uint32_t i;
    uint64_t rowHeadPos = omicsId * sampleSize, sampleSizeTmp;

    // fill the missing value
    for (auto i : NASignMarkCurr[omicsId]) {
        a[rowHeadPos + i] = 0;
    }
    sampleSizeTmp = sampleSize - NASignMarkCurr[omicsId].size();

    // calc the row mean
    rowSumTmp = 0;
    for (i = 0; i < sampleSize; i++) {
        rowSumTmp += a[rowHeadPos + i];
    }
    rowMeanTmp = rowSumTmp / sampleSizeTmp;

    // calc the row SD for omics
    rowSDTmp = 0;
    for (i = 0; i < sampleSize; i++){
        if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), i) == NASignMarkCurr[omicsId].end()) {
            rowSDTmp += pow(a[rowHeadPos + i] - rowMeanTmp, 2);
        }
    }
    rowSD[omicsId] = sqrt(rowSDTmp / (sampleSizeTmp - 1));

    // centralization
    if (rowSD[omicsId] != 0) { // constant variants
        for (i = 0; i < sampleSize; i++){
            if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), i) == NASignMarkCurr[omicsId].end()) {
                a[rowHeadPos + i] = (a[rowHeadPos + i] - rowMeanTmp) / rowSD[omicsId];
            }
        }
    } else {
        for (i = 0; i < sampleSize; i++){
            a[rowHeadPos + i] = 0;
        }
    }
}

// quantile based normalization
void cntrlQuant(float *omicsData, uint64_t omicsId, uint32_t sampleSize,
                vector<vector<uint32_t> >& NASignMarkCurr) {
    vector<float> v_temp(sampleSize);
    vector<float> v_rank(sampleSize);
    uint64_t rowHeadPos = omicsId * sampleSize;
    float maxV;

    uint32_t sampleSizeTmp = sampleSize - NASignMarkCurr[omicsId].size();

    // extract omics data
    maxV = omicsData[rowHeadPos + 0];
    for (uint32_t i = 0; i < sampleSize; i++) {
        v_temp[i] = omicsData[rowHeadPos + i];
        if (maxV < v_temp[i]) {
            maxV = v_temp[i];
        }
    }
    for (auto i : NASignMarkCurr[omicsId]) {
        v_temp[i] = maxV + 1;
    } // set missing value as max value to give them max rank

    // rank
    v_rank = rankSort(v_temp, sampleSize);
    
    // normalization
    for (uint32_t i = 0; i < sampleSize; i++) {
        omicsData[rowHeadPos + i] = float(gsl_cdf_ugaussian_Pinv((v_rank[i] + 0.5) / sampleSizeTmp)); // caution: v_rank start from 0!
    }
    
    // fill the missing value
    for (auto i : NASignMarkCurr[omicsId]) {
        omicsData[rowHeadPos + i] = 0;
    }

    return;
}

// sort and rank a vector
vector<float> rankSort(const vector<float>& v_temp, uint64_t sampleSize) {
    vector<std::pair<float, uint32_t> > v_sort(sampleSize);

    for (uint32_t i = 0; i < v_sort.size(); ++i) {
        v_sort[i] = std::make_pair(v_temp[i], i);
    }

    std::sort(v_sort.begin(), v_sort.end());

    vector<float> result(sampleSize);
    uint32_t currentRankP = 0, preRankP = -1;
    float currentValue, currentRank, currentQuantile;

    while (currentRankP < sampleSize) {
        currentValue = v_sort[currentRankP].first;
        currentRankP++;
        while (currentRankP < sampleSize) {
            if (v_sort[currentRankP].first != currentValue) {
                break;
            }
            currentRankP++;
        }
        currentRank = float(preRankP + currentRankP) / 2;
        for (uint32_t i = preRankP + 1; i < currentRankP; i++) {
            result[v_sort[i].second] = currentRank;
        }
        preRankP = currentRankP - 1;
    }

    return(result);
}

linearFitRlt linearFit(float corr, 
                       uint32_t omics1Id, uint32_t omics2Id, 
                       uint32_t sampleSize, uint32_t covarNum, 
                       float missingRateThd,
                       vector<vector<uint32_t> >& NASignMark1, vector<vector<uint32_t> >& NASignMark2, 
                       vector<double>& omics1NormSqrInv, 
                       vector<double>& omics1Scaling, vector<double>& omics2Scaling) {
    uint32_t df_r, df_t, sampleSizeTmp;
    vector<uint32_t> NASignMarkCurr;
    linearFitRlt rlt;
    rlt.omics1Id = omics1Id; rlt.omics2Id = omics2Id; // locus index start from 0

    // omit NA
    for (auto l : NASignMark1[omics1Id]) {
        NASignMarkCurr.push_back(l);
    }
    for (auto l : NASignMark2[omics2Id]) {
        NASignMarkCurr.push_back(l);
    }
    sort(NASignMarkCurr.begin(), NASignMarkCurr.end());
    NASignMarkCurr.erase(unique(NASignMarkCurr.begin(), NASignMarkCurr.end()), NASignMarkCurr.end());
    // QC by missing-rate threshold
    sampleSizeTmp = sampleSize - NASignMarkCurr.size();
    if (sampleSizeTmp <= (1 - missingRateThd) * sampleSize) {
      rlt.status = 1;
      return rlt;
    }

    // degree of freedom
    df_r = sampleSizeTmp - 1; // not concern with covariate when calc correlation 
    df_t = sampleSizeTmp - covarNum - 2;
    corr = corr / df_r;

    // calculate regression coefficients
    rlt.r2 = (float)pow(corr, 2);
    rlt.t = (float)sqrt(df_t) * corr / sqrt(1 - rlt.r2);
    rlt.p = gsl_cdf_tdist_Q(abs(rlt.t), df_t) * 2;
    rlt.b = (float)omics1NormSqrInv[omics1Id] * corr * df_r / omics1Scaling[omics1Id] * omics2Scaling[omics2Id];
    rlt.se = (float)rlt.b / rlt.t;
    rlt.nmiss = sampleSizeTmp;

    // handle error of out of range
    if (rlt.p < 1e-308){
        rlt.p = 1e-308;
    }

    // no error
    rlt.status = 0;

    return rlt;
}

void rCriticalValueCalc(double P, uint32_t sampleSize, double &rCriticalValue) {
    double tCriticalValue;
    tCriticalValue = gsl_cdf_tdist_Qinv(P/2, sampleSize - 2); // t critical value
    rCriticalValue = sqrt(pow(tCriticalValue, 2) / (sampleSize - 2 + pow(tCriticalValue, 2))) * (sampleSize - 1); // correlation critical value, times sampleSize
}

}  // namespace meqtllib

int main(int argc, char *argv[]) {
    string omics1FileName, omics2FileName, outputFileName;
    bool bfileFlag = false;
    double globalP = -1;
    string covarFileName;
    string NASign = "NA";
    float missingRateThd = 0.1;
    vector<double> distLv;
    vector<double> distLvP;
    int32_t threadMaxN = 1;
    int32_t omics1NormMod = 0, omics2NormMod = 0;

    auto cli = (
        required("--omics1") & value("omics1FileName", omics1FileName)        % "first omics file path",
        option("bfile").set(bfileFlag)         % "useing plink binary format represente first omics",
        required("--omics2") & value("omics2FileName", omics2FileName)       % "second omics file path",
        required("--out") & value("outputFileName", outputFileName)                % "output file path",
        option("-p") & number("globalP", globalP)                         % "global P-value threshold",
        option("--cov") & value("covarFileName", covarFileName)                % "covariates file path",
        option("--na") & value("NASign", NASign)                              % "sign of missing value",
        option("--missing-rate") & number("missingRateThd", missingRateThd) % "missing rate threshold",
        option("--dl") & numbers("distLv", distLv)     % "distance thresholds for each distance level",
        option("--dlp") & numbers("distLvP", distLvP)   % "P-value thresholds for each distance level",
        option("--threads") & integer("threadMaxN", threadMaxN)                       % "max. threads",
        option("--omics1norm")                                 % "normalization model for omics1 data"
            & (required("zscore").set(omics1NormMod, 1)
                                 | required("rank").set(omics1NormMod, 2)),
        option("--omics2norm")                                 % "normalization model for omics2 data"
            & (required("zscore").set(omics2NormMod, 1)
                                 | required("rank").set(omics2NormMod, 2))
    );

    auto fmt = doc_formatting{} .first_column(4) .doc_column(25) .last_column(80);

    if(!parse(argc, argv, cli)) {
        cout << make_man_page(cli, "fastQTLmapping", fmt)
        .prepend_section("DESCRIPTION", "    Fastest QTL mapping tool.")
        .append_section("LICENSE", "    GPL3") << '\n';
        return 0;
    }

    // Bonferroni correction
    if (globalP < 0) {
        globalP = 0.05 / omics1Num / omics2Num;
    }
    // fill distLvP with globalP
    if (distLvP.size() == 0) {
        distLvP.assign(distLv.size(), globalP);
    }
    // check equal length about distLv and distLvP
    if (distLv.size() != distLvP.size()) {
        cout << "Error: the length of distLv and distLvP is not equal.\n";
        return 1;
    }
    uint32_t distLvNum = distLvP.size();
    // check ascending about distLv
    for (uint32_t i = 1; i < distLv.size(); i++) {
        if (distLv[i] <= distLv[i - 1]) {
            cout << "Error: distLv is not ascending.\n";
            return 2;
        }
    }

    // global starting time stamp 
    double time_start_whole = omp_get_wtime(), time_end_whole;

    // initializing log file
    string outputLogFileName = outputFileName + ".log";
    ofstream outputLogFile;
    outputLogFile.open(outputLogFileName);

    // record paraments into log file
    outputLogFile << "omics 1 file : " << omics1FileName << endl;
    if (bfileFlag) {
        outputLogFile << "omics 1 is in PLINK binary format"<< endl;
    }
    outputLogFile << "omics 2 file : " << omics2FileName << endl;
    outputLogFile << "covar file : " << covarFileName << endl;
    outputLogFile << "missing value sign : " << NASign << endl;
    outputLogFile << "missing rate threshold : " << missingRateThd << endl;
    outputLogFile << "maximun parallel number : " << threadMaxN << endl;
    switch (omics1NormMod) {
        case 0 : 
            outputLogFile << "not normalize omics1 data" << endl;
            break;
        case 1 : 
            outputLogFile << "normalize omics1 data by z-score" << endl;
            break;
        case 2 : 
            outputLogFile << "normalize omics1 data by rank-based" << endl; 
            break;
        default : 
            break;
    }
    switch (omics2NormMod) {
        case 0 : 
            outputLogFile << "not normalize omics2 data" << endl;
            break;
        case 1 : 
            outputLogFile << "normalize omics2 data by z-score" << endl;
            break;
        case 2 : 
            outputLogFile << "normalize omics2 data by rank-based" << endl; 
            break;
        default : 
            break;
    }
    outputLogFile << endl;
    
    // calculate bfile size
    if (bfileFlag) { // input plink bfile as first omics
        calcBfileSize(omics1FileName, sampleSize, omics1Num);
    } else {
        calcInputSize(omics1FileName, sampleSize, omics1Num);
    }
    // calculate input file size
    calcInputSize(omics2FileName, sampleSize, omics2Num);
    // calculate input file size
    vector<bool> sampleFltSign(sampleSize, false);
    uint32_t covarNANum = 0;
    if (!covarFileName.empty()) {
        calcCovarSize(covarFileName, NASign, sampleSize, covarNum, sampleFltSign, covarNANum);
    }
    sampleSize -= covarNANum;

    // record data scale into log file
    outputLogFile << "omics 1 number : " << omics1Num << endl;
    outputLogFile << "omics 2 number : " << omics2Num << endl;
    outputLogFile << "covariates number : " << covarNum << endl;
    outputLogFile << "valid sample number : " << sampleSize << endl;
    outputLogFile << covarNANum << " samples are excluded because covariates missing" << endl;
    outputLogFile << endl;

    // critical value of t test
    vector<double> rCriticalValue(1 + distLvNum);
    for (uint32_t i = 0; i < distLvNum; i++) {
        rCriticalValueCalc(distLvP[i], sampleSize, rCriticalValue[i]);
        outputLogFile << "level " << i + 1 << " distance threshold: " << distLv[i] << endl; // level number start from 1
        outputLogFile << "level " << i + 1 << " pearson correlation significant threshold: " << distLvP[i] << endl;
        outputLogFile << "level " << i + 1 << " pearson correlation critical value: " << rCriticalValue[i] / (sampleSize - 1) << endl;
    }
    distLvP.push_back(globalP); // set last distLvP element as global P threshold
    rCriticalValueCalc(distLvP[distLvNum], sampleSize, rCriticalValue[distLvNum]); // set last rCriticalValue element as global r critical value
    outputLogFile << "global pearson correlation significant threshold: " << distLvP[distLvNum] << endl;
    outputLogFile << "global pearson correlation critical value: " << rCriticalValue[distLvNum] / (sampleSize - 1) << endl;
    outputLogFile << endl;

    // header of result file for output file
    ofstream outputFile;
    outputFile.open(outputFileName);
    outputFile << "omics1\t" << "omics2\t" << "BETA\t" << "SE\t" << "T\t" << "P\t" << "R2\t" << "NMISS\t" << "distance_level\n";
    outputFile.close();
    
    // input omics1 data
    float* omics1Data = (float*) mkl_malloc(sizeof(float) * omics1Num * sampleSize, 64);
    vector<vector<uint32_t> > NASignMark1(omics1Num, vector<uint32_t>(0)); // NA mark for first omics, N.O1 * NAs
    vector<string> omics1Name(omics1Num); // locus name for first omics
    vector<int32_t> omics1CHR(omics1Num); // CHR number for first omics
    vector<int64_t> omics1BP(omics1Num); // BP number for first omics
    if (bfileFlag) { // input plink bfile as first omics
        getBfileSNPid(omics1FileName, omics1Num, omics1Name, omics1CHR, omics1BP);
        snplib::UnpackGeno(omics1FileName, omics1Data, NASignMark1, sampleSize, omics1Num, sampleFltSign, covarNANum, omics2FileName);        
    } else {
        string* dataArea = new string[omics1Num];
        input2DfloatParse(omics1Data, omics1FileName, NASignMark1, NASign, 
                          omics1Name, omics1CHR, omics1BP, 
                          omics1Num, sampleSize, 
                          dataArea, threadMaxN, 
                          sampleFltSign, covarNANum);
        delete [] dataArea;
    }
    
    // time stamp for input omics1
    time_end_whole = omp_get_wtime();
    outputLogFile << "first omics data has been input, time used: " << time_end_whole - time_start_whole << " s" << endl;

    // input omics2 data
    float* omics2Data = (float*) mkl_malloc(sizeof(float) * omics2Num * sampleSize, 64);
    vector<vector<uint32_t> > NASignMark2(omics2Num, vector<uint32_t>(0)); // NA mark for second omics, N.O2 * NAs
    vector<string> omics2Name(omics2Num); // locus name for second omics
    vector<int32_t> omics2CHR(omics2Num); // CHR number for second omics
    vector<int64_t> omics2BP(omics2Num); // BP number for second omics
    string* dataArea = new string[omics2Num];
    input2DfloatParse(omics2Data, omics2FileName, NASignMark2, NASign, 
                      omics2Name, omics2CHR, omics2BP, 
                      omics2Num, sampleSize, 
                      dataArea, threadMaxN, 
                      sampleFltSign, covarNANum);
    delete [] dataArea;
    // time stamp for input omics2
    time_end_whole = omp_get_wtime();
    outputLogFile << "second omics data has been input, time used: " << time_end_whole - time_start_whole << " s" << endl;
    
    // input covariates data
    float* covarData = (float*) mkl_malloc(sizeof(float) * covarNum * sampleSize, 64);
    vector<vector<uint32_t> > NASignMarkC(covarNum, vector<uint32_t>(0)); // NA mark for first omics, SNPnum * NAs
    if (covarNum > 0) {
        inputCovar(covarData, covarFileName, 
                   covarNum, sampleSize, sampleFltSign, covarNANum);    
    }

    // orthogonal projection start

    // first centralization of SNP data
    vector<double> omics1RowSDCntrl1(omics1Num);
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics1Num, \
           omics1Data, \
           sampleSize, \
           omics1RowSDCntrl1, NASignMark1)
    for (uint64_t i = 0; i < omics1Num; i++) {
        cntrl(omics1Data, i, sampleSize, omics1RowSDCntrl1, NASignMark1);
    }
    
    // first centralization of trait data
    vector<double> omics2RowSDCntrl1(omics2Num);
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics2Num, \
           omics2Data, \
           sampleSize, \
           omics2RowSDCntrl1, NASignMark2)
    for (uint64_t i = 0; i < omics2Num; i++) {
        cntrl(omics2Data, i, sampleSize, omics2RowSDCntrl1, NASignMark2);
    }
    
    if (covarNum > 0) {
        // first centralization of covariates
        vector<double> covarRowSDCntrl1(covarNum);
    #   pragma omp parallel for \
        num_threads(threadMaxN) \
        shared(covarNum, \
               covarData, \
               sampleSize, \
               covarRowSDCntrl1)
        for (uint64_t i = 0; i < covarNum; i++) {
            cntrl(covarData, i, sampleSize, covarRowSDCntrl1, NASignMarkC);
        }

        // transpose covar
        float* covarDataT = (float*) mkl_malloc(sizeof(float) * sampleSize * covarNum, 64);
        for (uint64_t i = 0; i < covarNum; i++)
            for (uint64_t j = 0; j < sampleSize; j++)
                covarDataT[j * covarNum + i] = covarData[i * sampleSize + j];
        
        // QR decompose
        float *tau = new float[covarNum];
        LAPACKE_sgeqrf(LAPACK_ROW_MAJOR, sampleSize, covarNum, covarDataT, covarNum, tau);
        LAPACKE_sorgqr(LAPACK_ROW_MAJOR, sampleSize, covarNum, covarNum, covarDataT, covarNum, tau);
        
        // backup omics data
        float* omics1DataTmp = (float*) mkl_malloc(sizeof(float) * omics1Num * sampleSize, 64); 
        memcpy(omics1DataTmp, omics1Data, sizeof(float) * omics1Num * sampleSize);
        float* omics2DataTmp = (float*) mkl_malloc(sizeof(float) * omics2Num * sampleSize, 64); 
        memcpy(omics2DataTmp, omics2Data, sizeof(float) * omics2Num * sampleSize);

        // projection
        float* covarNormSqr = (float*) mkl_malloc(sizeof(float) * sampleSize * sampleSize, 64);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, sampleSize, sampleSize, covarNum, 1, covarDataT, covarNum, covarDataT, covarNum, 0, covarNormSqr, sampleSize);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, omics1Num, sampleSize, sampleSize, -1, omics1DataTmp, sampleSize, covarNormSqr, sampleSize, 1, omics1Data, sampleSize);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, omics2Num, sampleSize, sampleSize, -1, omics2DataTmp, sampleSize, covarNormSqr, sampleSize, 1, omics2Data, sampleSize);

        // free backup of omics data
        mkl_free(omics1DataTmp); mkl_free(omics2DataTmp);
    }

    // second centralization
    vector<double> omics1RowSDCntrl2(omics1Num);
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics1Num, \
           omics1Data, \
           sampleSize, \
           omics1RowSDCntrl2, NASignMark1)
    for (uint64_t i = 0; i < omics1Num; i++) {
        if (omics1NormMod != 2) {
            cntrl(omics1Data, i, sampleSize, omics1RowSDCntrl2, NASignMark1);
        } else {
            cntrlQuant(omics1Data, i, sampleSize, NASignMark1);
        }    
    }
    vector<double> omics2RowSDCntrl2(omics2Num);
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics2Num, \
           omics2Data, \
           sampleSize, \
           omics2RowSDCntrl2, NASignMark2)
    for (uint64_t i = 0; i < omics2Num; i++) {
        if (omics2NormMod != 2) {
            cntrl(omics2Data, i, sampleSize, omics2RowSDCntrl2, NASignMark2);
        } else {
            cntrlQuant(omics2Data, i, sampleSize, NASignMark2);
        }    
    }

    // Scalings for beta coefficents
    vector<double> omics1Scaling(omics1Num, 1);
    if (omics1NormMod == 0) {
        for (uint64_t i = 0; i < omics1Num; i++) {
            omics1Scaling[i] = omics1RowSDCntrl1[i] * omics1RowSDCntrl2[i];
        }
    }
    vector<double> omics2Scaling(omics2Num, 1);
    if (omics2NormMod == 0) {
        for (uint64_t i = 0; i < omics2Num; i++) {
            omics2Scaling[i] = omics2RowSDCntrl1[i] * omics2RowSDCntrl2[i];
        }
    }
    
    // inverse of norm 2 of omics1
    vector<double> omics1NormSqrInv(omics1Num);
    for (uint64_t i = 0; i < omics1Num; i++) {
        omics1NormSqrInv[i] = 1 / cblas_sdot(sampleSize, &omics1Data[i * sampleSize], 1, &omics1Data[i * sampleSize], 1);
    }

    // time stamp for orthogonal projection
    time_end_whole = omp_get_wtime();
    outputLogFile << "orthogonal projection has finished, time used: " << time_end_whole - time_start_whole << " s" << endl;

    // orthogonal projection finished

    // split Block schedule
    uint64_t omics1BlockStride = min(omics1Num, blockSize);
    uint64_t omics2BlockStride = min(omics2Num, blockSize);
    uint64_t omics1BlockNum = (int)(omics1Num + omics1BlockStride - 1) / omics1BlockStride;
    uint64_t omics2BlockNum = (int)(omics2Num + omics2BlockStride - 1) / omics2BlockStride;
    uint64_t blockNum = omics1BlockNum * omics2BlockNum;

    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "preprocessing has finished, time used: " << time_end_whole - time_start_whole << " s" << endl;

#   pragma omp parallel \
    num_threads(min((uint64_t)threadMaxN, blockNum)) \
    shared(omics1Name, omics2Name, \
           omics1CHR, omics2CHR, omics1BP, omics2BP, \
           omics1Num, omics2Num, covarNum, sampleSize, missingRateThd, \
           omics1Data, omics2Data, \
           omics1Scaling, omics2Scaling, \
           omics1BlockStride, omics2BlockStride, \
           NASignMark1, NASignMark2, \
           distLv, distLvP, distLvNum, \
           rCriticalValue, \
           outputFileName, \
           outputLogFile)
    {
        // mark of omics id in current block
        uint64_t omics1BlockHead, omics2BlockHead;

        // alloc space for current result on each thread
        linearFitRlt* rltArr = new linearFitRlt[omics1BlockStride * omics2BlockStride];
        float* corr = (float*) mkl_malloc(sizeof(float) * omics1BlockStride * omics2BlockStride, 64);

#       pragma omp for schedule(dynamic)
        for (uint32_t i = 0; i < blockNum; i++) {
            omics1BlockHead = (int) i % omics1BlockNum * omics1BlockStride;
            omics2BlockHead = (int) i / omics1BlockNum * omics2BlockStride;

            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                        min(omics1BlockStride, omics1Num - omics1BlockHead), min(omics2BlockStride, omics2Num - omics2BlockHead), sampleSize, 
                        1, &omics1Data[omics1BlockHead * sampleSize], sampleSize, &omics2Data[omics2BlockHead * sampleSize], sampleSize, 
                        0, corr, min(omics2BlockStride, omics2Num - omics2BlockHead));
            // number of items in current result
            uint64_t sigCnt = 0;
            uint64_t corrP = 0;
            for (uint64_t j = omics1BlockHead; j < min(omics1BlockHead + omics1BlockStride, omics1Num); j++) {
                for (uint64_t k = omics2BlockHead; k < min(omics2BlockHead + omics2BlockStride, omics2Num); k++) {
                    int32_t levelTmp = -1;
                    float corrCurr = corr[corrP++];
                    if (distLvNum >= 1 &&
                        omics1CHR[j] == omics2CHR[k] && 
                        omics1CHR[j] > 0 && omics2CHR[k] > 0 && 
                        abs(omics1BP[j] - omics2BP[k]) <= distLv[distLvNum - 1] && 
                        omics1BP[j] > 0 && omics2BP[k] > 0) { // two locus are localed in the broadest level distance
                        for (uint32_t l = 0; l < distLvNum; l++) {
                            if (abs(omics1BP[j] - omics2BP[k]) <= distLv[l] && 
                                abs(corrCurr) >= rCriticalValue[l]) {
                                levelTmp = l; break;
                            }
                        }
                    } else { // two locus are not localed in any level distance
                        if (abs(corrCurr) >= rCriticalValue[distLvNum]) {
                            levelTmp = distLvNum;
                        }
                    }
  
                    if (levelTmp >= 0) {
                        linearFitRlt rltTmp = 
                          linearFit(corrCurr, 
                                    j, k, 
                                    sampleSize, covarNum, 
                                    missingRateThd,
                                    NASignMark1, NASignMark2, 
                                    omics1NormSqrInv, 
                                    omics1Scaling, omics2Scaling);

                        if (rltTmp.p <= distLvP[levelTmp] && rltTmp.status == 0) {
                            rltArr[sigCnt] = rltTmp;
                            rltArr[sigCnt].level = levelTmp + 1; // level number start from 1
                            sigCnt++;
                        }
                    }
                }
            }

            // P< P1, maker BETA SE R2 T P
            // P1 < P < P2, maker BETA T P
            // P > P2, do not output
#           pragma omp critical
            {
                ofstream outputFile;
                outputFile.open(outputFileName, ios::app);
                outputFile << setprecision(2);
                for (uint32_t j = 0; j < sigCnt; j++){
                        outputFile << omics1Name[rltArr[j].omics1Id] << "\t";
                        outputFile << omics2Name[rltArr[j].omics2Id] << "\t";
                        outputFile << rltArr[j].b << "\t";
                        outputFile << rltArr[j].se << "\t";
                        outputFile << rltArr[j].t << "\t";
                        outputFile << rltArr[j].p << "\t";
                        outputFile << rltArr[j].r2 << "\t";
                        outputFile << rltArr[j].nmiss << "\t";
                        outputFile << rltArr[j].level << "\n";
                }
                outputFile.close();
            }

            double time_end = omp_get_wtime();
            outputLogFile << i + 1 << " / " << blockNum << " block has finished, time used: " << time_end - time_start_whole << " s" << endl;
        }
    }

    // whole calc time stamp
    time_end_whole = omp_get_wtime();
    outputLogFile << "whole script time used: " << time_end_whole - time_start_whole << " s" << endl;
    outputLogFile << endl;
    outputLogFile.close();

    mkl_free(omics1Data); mkl_free(omics2Data);

    return 0;
}
