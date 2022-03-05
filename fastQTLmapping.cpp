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
using std::pair;
using std::make_pair;
using std::sort;

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

void UnpackGeno(string bfileNameRoot, double* geno_data, vector<vector<uint32_t> >& NA_data, size_t num_samples, size_t num_snps, 
                vector<bool>& sampleFltSign, uint32_t covarNANum, string teststr) {
  ifstream inputFile;
  double *snp_geno_d = new double[num_samples + covarNANum];
  double *snp_mask_d = new double[num_samples + covarNANum];
  inputFile.open(bfileNameRoot + ".bed", ios::in | ios::binary);

  std::istreambuf_iterator<char> inputFile_it(inputFile);
  inputFile_it++; inputFile_it++; inputFile_it++;
  vector<uint8_t> geno( (inputFile_it),
                          std::istreambuf_iterator<char>() );

  const array<double, 4> geno_table{2.0, 0.0, 1.0, 0.0}; // Homozygote A1, missing, Heterozygote, Homozygote A2
  const array<double, 4> mask_table{0.0, 1.0, 0.0, 0.0};
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
    string s, oneItem;
    ifstream inputFile;

    inputFile.open(bfileNameRoot + ".bim", ios::in);
    for (uint32_t i = 0; i < num_snps; i++) {
        getline(inputFile, s);
        istringstream is(s);
        
        is >> oneItem; // CHR
        try {
            omicsCHR[i] = stoi(oneItem);
        } catch (std::invalid_argument) {
            omicsCHR[i] = -1;
        }
        
        is >> oneItem; // SNP
        omicsName[i] = oneItem;
        
        is >> oneItem; // CM
        
        is >> oneItem; // BP
        try {
            omicsBP[i] = stol(oneItem);
        } catch (std::invalid_argument) {
            omicsBP[i] = -1;
        }
    }
    inputFile.close();
}

void calcInputSize(string omicsFileName, uint32_t &sampleSize, uint64_t& omicsNum) {
    uint32_t i, rowsCount = 0, colsCount = 0;
    ifstream inputFile;
    string s, oneItem;

    inputFile.open(omicsFileName);
    assert(inputFile.is_open());
    
    getline(inputFile, s);
    istringstream is(s);
    while (is >> oneItem) {
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

//input 2 dimension omics data
void input2DfloatParse(double* omicsData, string fileName, vector<vector<uint32_t> >& NASignMark, string NASign, 
                  vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP, 
                  uint64_t omicsNum, uint32_t sampleSize, 
                  string* dataArea, 
                  uint32_t threadMaxN, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum) {
    ifstream inputFile;
    string one_line, oneItem;
    string delimiter = " \t";
    string::size_type pos,lastPos;
    inputFile.open(fileName);
    assert(inputFile.is_open());
    for (uint64_t i = 0; i < omicsNum; i++) {
        getline(inputFile, one_line);

        lastPos = one_line.find_first_not_of(delimiter, 0);
        pos = one_line.find_first_of(delimiter, lastPos);
        oneItem = one_line.substr(lastPos, pos - lastPos); // SNP
        omicsName[i] = oneItem;
        
        lastPos = one_line.find_first_not_of(delimiter, pos);
        pos = one_line.find_first_of(delimiter, lastPos);
        oneItem = one_line.substr(lastPos, pos - lastPos); // CHR
        try {
            omicsCHR[i] = stoi(oneItem);
        } catch (std::invalid_argument) {
            omicsCHR[i] = -1;
        }
        
        lastPos = one_line.find_first_not_of(delimiter, pos);
        pos = one_line.find_first_of(delimiter, lastPos);
        oneItem = one_line.substr(lastPos, pos - lastPos); // BP
        try {
            omicsBP[i] = stol(oneItem);
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
                    double k = 0.1;
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
                    omicsData[rowHeadPos + sid] = sb * realData;
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

void calcCovarSize(string covarFileName, string NASign, uint64_t sampleSize, uint32_t& covarNum, 
                   vector<bool>& sampleFltSign, uint32_t& covarNANum, 
                   vector<int32_t>& categFlag, uint32_t& covarCategNum) {
    uint32_t i, rowsCount = 0;
    covarCategNum = 0;
    ifstream inputFile;
    string s, oneItem;
    vector<int32_t> categFlagRef;

    inputFile.open(covarFileName);
    assert(inputFile.is_open());
    while (getline(inputFile, s)){
        rowsCount++;
        if (find(categFlag.begin(), categFlag.end(), rowsCount) != categFlag.end()) { // categFlag starts from 1
            covarCategNum++;
            categFlagRef.push_back(rowsCount - 1);
        }
        istringstream is(s);
        for (i = 0; i < sampleSize; i++) {
            is >> oneItem;
            if (oneItem == NASign){
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

    // extract legal rows id of categorical covariates
    // unique 
    // sort
    // minus 1 to make it start from 0
    categFlag.assign(categFlagRef.begin(), categFlagRef.end());
}

//input 2 dimension covariates data
double* inputCovar(string fileName, 
                  uint32_t& covarNum, uint32_t sampleSize, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum, 
                  vector<int32_t>& categFlag, uint32_t covarCategNum) {
    uint32_t i, j, k, sid;
    ifstream inputFile;
    string s, oneItem;

    // transcript categorical covariates into int vector
    // & calculate the number of levels of each categorical covariates
    int categVar[covarCategNum * sampleSize]; // categVar starts from 0
    vector<string> categVarTmp;
    vector<uint32_t> categLevN;
    uint32_t rowsCount = 0, categIdx = 0;
    uint32_t categLevSum = 0;

    inputFile.open(fileName);
    for (auto categFlagCurr : categFlag) {
        for (; rowsCount <= categFlagCurr; rowsCount++) {
            getline(inputFile, s);
        }
        istringstream is(s);
        sid = 0;
        for (j = 0; j < sampleSize + covarNANum; j++) {
            is >> oneItem;
            if (!sampleFltSign[j]) {
                auto p = find(categVarTmp.begin(), categVarTmp.end(), oneItem);
                if (p == categVarTmp.end()) {
                    categVarTmp.push_back(oneItem);
                    categVar[categIdx * sampleSize + sid] = categVarTmp.size() - 1;
                } else {
                    categVar[categIdx * sampleSize + sid] = p - categVarTmp.begin();
                }
                sid++;
            }
        }

        s.clear(); is.str(""); is.clear();
        categLevN.push_back(categVarTmp.size());
        categVarTmp.clear();
        categLevSum += *(categLevN.end() - 1);
        categIdx++;
    }
    inputFile.close();

    // input numeric and categorical covariates
    double* covarData = (double*) mkl_malloc(sizeof(double) * 
                                           ((covarNum - covarCategNum) + (categLevSum - covarCategNum)) * sampleSize, 
                                           64);
    uint32_t lastCovData = 0;
    categIdx = 0;
    inputFile.open(fileName);
    for (i = 0; i < covarNum; i++) {
        getline(inputFile, s);
        istringstream is(s);
        sid = 0;

        if (categIdx < covarCategNum && categFlag[categIdx] == i) {
            for (j = 0; j < sampleSize + covarNANum; j++) {
                if (!sampleFltSign[j]) {
                    for (k = 0; k < categLevN[categIdx] - 1; k++) { // only first categLevN[categIdx] - 1 covariates will be coded as dummy
                        covarData[(lastCovData + k) * sampleSize + sid] = 0;
                    }
                    auto categVarCur = categVar[categIdx * sampleSize + sid];
                    if (categVarCur < categLevN[categIdx] - 1) {
                        covarData[(lastCovData + categVarCur) * sampleSize + sid] = 1;
                    }
                    sid++;
                }
            }
            lastCovData += categLevN[categIdx] - 1;
            categIdx++;
        } else {
            for (j = 0; j < sampleSize + covarNANum; j++) {
                is >> oneItem;
                if (!sampleFltSign[j]) {
                    covarData[lastCovData * sampleSize + sid] = stod(oneItem);
                    sid++;
                }
            }
            lastCovData++;
        }
        
        s.clear(); is.str(""); is.clear();
    }
    inputFile.close();

    covarNum = (covarNum - covarCategNum) + (categLevSum - covarCategNum);

    return(covarData);
}

template <class ForwardIterator, class T>
uint64_t binarySearch(ForwardIterator head, ForwardIterator tail, const T& val) {
    ForwardIterator it, headBak = head;
    typename std::iterator_traits<ForwardIterator>::difference_type count, step;
    count = distance(head, tail);
    while (count > 0)
    {
        it = head; step=count/2; advance (it,step);
        if (*it < val) {
            head = ++it;
            count -= step+1;
        }
        else count = step;
    }
    if (!(head == tail) && !(val < *head)) {
        return(distance(headBak, head));
    } else {
        return(-1);
    }
}

void inputRplList(string rplFileName, vector<pair<uint64_t, uint64_t> >& rplList, vector<string>& omics1Name, vector<string>& omics2Name, uint32_t threadMaxN) {
    ifstream inputFile;
    string s;
    
    // input rplFile
    vector<pair<string, string> > rplFile;
    inputFile.open(rplFileName);
    assert(inputFile.is_open());
    while (getline(inputFile, s)) {
        istringstream is(s);
        pair<string, string> oneItem;
        is >> oneItem.first >> oneItem.second;
        rplFile.push_back(oneItem);
    }
    rplList.resize(rplFile.size());

    // parallel binary search
#   pragma omp parallel \
    num_threads(threadMaxN)
    {
        // build rank vector of omics name
        vector<pair<string, uint64_t> > omics1NameSortTmp(omics1Name.size()); // combine omics id and order
        for (uint64_t i = 0; i < omics1Name.size(); ++i) { omics1NameSortTmp[i] = make_pair(omics1Name[i], i); }
        sort(omics1NameSortTmp.begin(), omics1NameSortTmp.end());
        vector<string> omics1NameSort(omics1Name.size()); vector<uint64_t> omics1NameRank(omics1Name.size()); // split sorted omics id and rank
        for (uint64_t i = 0; i < omics1Name.size(); ++i) { omics1NameSort[i] = omics1NameSortTmp[i].first; }
        for (uint64_t i = 0; i < omics1Name.size(); ++i) { omics1NameRank[i] = omics1NameSortTmp[i].second; }

        vector<pair<string, uint64_t> > omics2NameSortTmp(omics2Name.size()); // combine omics id and order
        for (uint64_t i = 0; i < omics2Name.size(); ++i) { omics2NameSortTmp[i] = make_pair(omics2Name[i], i); }
        sort(omics2NameSortTmp.begin(), omics2NameSortTmp.end());
        vector<string> omics2NameSort(omics2Name.size()); vector<uint64_t> omics2NameRank(omics2Name.size()); // split sorted omics id and rank
        for (uint64_t i = 0; i < omics2Name.size(); ++i) { omics2NameSort[i] = omics2NameSortTmp[i].first; }
        for (uint64_t i = 0; i < omics2Name.size(); ++i) { omics2NameRank[i] = omics2NameSortTmp[i].second; }

#       pragma omp parallel for schedule(static)
        for (uint64_t i = 0; i < rplFile.size(); i++) {
            pair<uint64_t, uint64_t> oneItem;
            oneItem.first = binarySearch(omics1NameSort.begin(), omics1NameSort.end(), rplFile[i].first);
            oneItem.second = binarySearch(omics2NameSort.begin(), omics2NameSort.end(), rplFile[i].second);
            if (oneItem.first >= 0 && oneItem.second >= 0) {
                rplList[i] = make_pair(omics1NameRank[oneItem.first], omics2NameRank[oneItem.second]);
            } else {
                rplList[i] = make_pair(-1, -1);
            }
        }
    }
}

// centralize one row of a matrix
void cntrl(double* a, uint64_t omicsId, 
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
void cntrlQuant(double *omicsData, uint64_t omicsId, uint32_t sampleSize,
                vector<vector<uint32_t> >& NASignMarkCurr) {
    vector<double> v_temp(sampleSize);
    vector<double> v_rank(sampleSize);
    uint64_t rowHeadPos = omicsId * sampleSize;
    double maxV;

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
        omicsData[rowHeadPos + i] = double(gsl_cdf_ugaussian_Pinv((v_rank[i] + 0.5) / sampleSizeTmp)); // caution: v_rank start from 0!
    }
    
    // fill the missing value
    for (auto i : NASignMarkCurr[omicsId]) {
        omicsData[rowHeadPos + i] = 0;
    }

    return;
}

// sort and rank a vector
vector<double> rankSort(const vector<double>& v_temp, uint64_t sampleSize) {
    vector<pair<double, uint64_t> > v_sort(sampleSize);

    for (uint64_t i = 0; i < v_sort.size(); ++i) {
        v_sort[i] = make_pair(v_temp[i], i);
    }

    sort(v_sort.begin(), v_sort.end());

    vector<double> result(sampleSize);
    uint64_t currentRankP = 0, preRankP = -1;
    double currentValue, currentRank, currentQuantile;

    while (currentRankP < sampleSize) {
        currentValue = v_sort[currentRankP].first;
        currentRankP++;
        while (currentRankP < sampleSize) {
            if (v_sort[currentRankP].first != currentValue) {
                break;
            }
            currentRankP++;
        }
        currentRank = double(preRankP + currentRankP) / 2;
        for (uint64_t i = preRankP + 1; i < currentRankP; i++) {
            result[v_sort[i].second] = currentRank;
        }
        preRankP = currentRankP - 1;
    }

    return(result);
}

linearFitRlt linearFit(double corr, 
                       uint64_t omics1Id, uint64_t omics2Id, uint32_t level, 
                       uint32_t sampleSize, uint32_t covarNum, 
                       float missingRateThd,
                       double* omics1Data, double* omics2Data, double* covarData, 
                       vector<vector<uint32_t> >& NASignMark1, vector<vector<uint32_t> >& NASignMark2, 
                       vector<double>& omics1Sum, vector<double>& omics2Sum, vector<double>& covarSum, 
                       vector<double>& omics1Sqr, vector<double>& omics1OrtgSqrInv, 
                       vector<vector<double> >& omics1DotCov, vector<vector<double> >& omics2DotCov, vector<vector<double> >& CovarInter,
                       vector<double>& omics1Scaling, vector<double>& omics2Scaling, 
                       vector<double>& omics1RowSDCntrl1, vector<double>& omics2RowSDCntrl1) {
    uint32_t df_r, df_t, sampleSizeCurr;
    vector<uint32_t> NASignMarkCurr; NASignMarkCurr.reserve(sampleSize);
    linearFitRlt rlt;
    rlt.omics1Id = omics1Id; rlt.omics2Id = omics2Id; // locus index start from 0
    rlt.level = level + 1; // level number start from 1
    double b, se, t, p;

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
    sampleSizeCurr = sampleSize - NASignMarkCurr.size();
    if (NASignMarkCurr.size() > missingRateThd * sampleSize) {
      rlt.status = 1;
      return rlt;
    }
    rlt.nmiss = sampleSizeCurr;

    if (NASignMarkCurr.size() == 0) {
        // no missing-value, produce statistics from correlation
        // degree of freedom
        df_r = sampleSizeCurr - 1; // not concern with covariate when calc correlation 
        df_t = sampleSizeCurr - (1 + covarNum) - 1;
        corr = corr / df_r;

        // calculate regression coefficients
        double r2 = pow(corr, 2);
        t = sqrt(df_t) * corr / sqrt(1 - r2);
        p = gsl_cdf_tdist_Q(abs(t), df_t) * 2;
        b = omics1OrtgSqrInv[omics1Id] * corr * df_r / omics1Scaling[omics1Id] * omics2Scaling[omics2Id];
        se = b / t;
    } else {
        // Gauss-Jordan algorithm
        double* omics1DataCurr = (double*) mkl_malloc(sizeof(double) * sampleSize, 64);
        double* omics2DataCurr = (double*) mkl_malloc(sizeof(double) * sampleSize, 64);
        memcpy(omics1DataCurr, &omics1Data[omics1Id * sampleSize], sizeof(double) * sampleSize);
        memcpy(omics2DataCurr, &omics2Data[omics2Id * sampleSize], sizeof(double) * sampleSize);
        double* A = (double*) mkl_malloc(sizeof(double) * (covarNum + 2) * (covarNum + 2), 64); // column major; lower triangle
        double* B = (double*) mkl_malloc(sizeof(double) * (covarNum + 2), 64); // column major

        // generate matrix A and B
        A[0] = sampleSizeCurr;
        A[1] = omics1Sum[omics1Id];
        A[1 + (covarNum + 2)] = omics1Sqr[omics1Id];
        B[0] = omics2Sum[omics2Id];        
        B[1] = cblas_ddot(sampleSize, omics1DataCurr, 1, omics2DataCurr, 1); // do not need to remove NA
        for (uint32_t i = 0; i < covarNum; i++) {
            A[2 + i] = covarSum[i];
            A[(2 + i) + (2 + covarNum)] = omics1DotCov[omics1Id][i];
            B[2 + i] = omics2DotCov[omics2Id][i];
            for (uint32_t j = 0; j <= i; j++) {
                A[(2 + i) + (2 + j) * (covarNum + 2)] = CovarInter[i][j];
            }
        }

        // remove NA element which missing in only one omic
        for (auto l : NASignMarkCurr) {
            auto omics1NACurr = omics1DataCurr[l];
            auto omics2NACurr = omics2DataCurr[l];
            A[1] -= omics1NACurr;
            A[(covarNum + 2) + 1] -= omics1NACurr * omics1NACurr;
            B[0] -= omics2NACurr;
            for (uint32_t i = 0; i < covarNum; i++) {
                auto covarNACurr = covarData[i * sampleSize + l];
                A[2 + i] -= covarNACurr;
                A[(2 + i) + (2 + covarNum)] -= omics1NACurr * covarNACurr;
                B[2 + i] -= omics2NACurr * covarNACurr;
                for (uint32_t j = 0; j <= i; j++) {
                    A[(2 + i) + (2 + j) * (covarNum + 2)] -= covarNACurr * covarData[j * sampleSize + l];
                }
            }
        }

        // solve lm
        LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', 2 + covarNum, A, 2 + covarNum);
        LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', 2 + covarNum, 1, 
                       A , 2 + covarNum, 
                       B , 2 + covarNum);
        LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', 2 + covarNum, A, 2 + covarNum);

        // beta = beta.norm * sd(y) / sd(x1)
        b = B[1];
        df_t = sampleSizeCurr - (1 + covarNum) - 1;
        // calc y_het
        double* ivMt = (double*) mkl_malloc(sizeof(double) * sampleSize * (covarNum + 2), 64); // Independent Variable Matrix, column-major
        for (uint32_t i = 0; i < sampleSize; i++) { ivMt[i] = 1; }
        memcpy(&ivMt[0 + sampleSize], omics1DataCurr, sizeof(double) * sampleSize);
        memcpy(&ivMt[0 + 2*sampleSize], covarData, sizeof(double) * covarNum * sampleSize);
        double* resid = (double*) mkl_malloc(sizeof(double) * sampleSize, 64); // residuals vector
        memcpy(resid, omics2DataCurr, sizeof(double) * sampleSize);
        // calc residuals
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            sampleSize, 1, covarNum + 2, 
            1, ivMt, sampleSize, 
            B, covarNum + 2, 
            -1, resid, sampleSize);
        for (auto l : NASignMarkCurr) { resid[l] = 0; } // mask out residuals at NA sample
        double SSE = 0;
        for (uint32_t i = 0; i < sampleSize; i++) { SSE += pow(resid[i], 2); }
        double MSE = SSE / df_t;
        se = sqrt(MSE * A[1 + covarNum + 2]);
        t = b / se;
        p = gsl_cdf_tdist_Q(abs(t), df_t) * 2;
        // fix data scale
        b = b * omics2RowSDCntrl1[omics2Id] / omics1RowSDCntrl1[omics1Id];
        se = b / t;

        mkl_free(omics1DataCurr); mkl_free(omics2DataCurr);
        mkl_free(A); mkl_free(B);
        mkl_free(ivMt); mkl_free(resid);
    }

    // handle error of out of range
    if (p < 1e-308){ p = 1e-308; }

    // cummarize results
    rlt.b = (float)b; rlt.se = (float)se; rlt.t = (float)t; rlt.p = p;
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
    string omics1FileName, omics2FileName, outputFileName, rplFileName = "NA";
    bool bfileFlag = false;
    double globalP = -1;
    string covarFileName;
    string NASign = "NA";
    float missingRateThd = 0.1;
    vector<int32_t> categFlag;
    vector<double> distLv;
    vector<double> distLvP;
    int32_t threadMaxN = 1;
    int32_t omics1NormMod = 0, omics2NormMod = 0;
    double PLooseMarg = 100;

    auto cli = (
        required("--omics1") & value("omics1FileName", omics1FileName)       % "first omics file path",
        option("bfile").set(bfileFlag)           % "useing plink binary format represente first omics",
        required("--omics2") & value("omics2FileName", omics2FileName)      % "second omics file path",
        required("--out") & value("outputFileName", outputFileName)               % "output file path",
        option("-p") & number("globalP", globalP)                         % "global P-value threshold",
        option("--ploose") & number("PLooseMarg", PLooseMarg) %   "margin for loose P-value threshold",
        option("--cov") & value("covarFileName", covarFileName)               % "covariates file path",
        option("--categ") & integers("categFlag", categFlag)  % "flag indicate categorical covariates",
        option("--na") & value("NASign", NASign)                             % "sign of missing value",
        option("--missing-rate") & number("missingRateThd", missingRateThd) % "missing rate threshold",
        option("--dl") & numbers("distLv", distLv)     % "distance thresholds for each distance level",
        option("--dlp") & numbers("distLvP", distLvP)   % "P-value thresholds for each distance level",
        option("--threads") & integer("threadMaxN", threadMaxN)                       % "max. threads",
        option("--omics1norm")                                 % "normalization model for omics1 data"
            & (required("zscore").set(omics1NormMod, 1)
                                 | required("rank").set(omics1NormMod, 2)),
        option("--omics2norm")                                 % "normalization model for omics2 data"
            & (required("zscore").set(omics2NormMod, 1)
                                 | required("rank").set(omics2NormMod, 2)),
        option("--rpl") & value("rplFileName", rplFileName)             % "replication list file path"
    );

    auto fmt = doc_formatting{} .first_column(4) .doc_column(25) .last_column(80);

    if(!parse(argc, argv, cli)) {
        cout << make_man_page(cli, "fastQTLmapping", fmt)
        .prepend_section("DESCRIPTION", "    Fastest QTL mapping tool. Version 1.9.1")
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
    outputLogFile << "fastQTLmapping version 1.9.1 start" << endl << endl;

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

    // Specifies the global number of threads for MKL
    mkl_set_num_threads(1);
    
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
    uint32_t covarCategNum = 0;
    if (!covarFileName.empty()) {
        calcCovarSize(covarFileName, NASign, sampleSize, covarNum, sampleFltSign, covarNANum, categFlag, covarCategNum);
    }
    sampleSize -= covarNANum;

    // check the range of categFlag is legal
    for (uint32_t i = 1; i < distLv.size(); i++) {
        if (distLv[i] <= distLv[i - 1]) {
            cout << "Error: distLv is not ascending.\n";
            return 2;
        }
    }

    // record data scale into log file
    outputLogFile << "omics 1 number : " << omics1Num << endl;
    outputLogFile << "omics 2 number : " << omics2Num << endl;
    outputLogFile << "covariates number : " << covarNum << endl;
    outputLogFile << "numeric covariates number : " << covarNum - covarCategNum << endl;
    outputLogFile << "categorical covariates number : " << covarCategNum << endl;
    outputLogFile << "valid sample number : " << sampleSize << endl;
    outputLogFile << covarNANum << " samples are excluded because covariates missing" << endl;
    outputLogFile << endl;

    // critical value of t test
    vector<double> rCriticalValue(1 + distLvNum);
    for (uint32_t i = 0; i < distLvNum; i++) {
        rCriticalValueCalc(min(distLvP[i] * PLooseMarg, 1.0), sampleSize, rCriticalValue[i]);
        outputLogFile << "level " << i + 1 << " distance threshold: " << 
            distLv[i] << endl; // level number start from 1
        outputLogFile << "level " << i + 1 << " significant threshold: P-value <= " << 
            distLvP[i] << endl;
        outputLogFile << "level " << i + 1 << " loose significant threshold: P-value <= " << 
            min(distLvP[i] * PLooseMarg, 1.0) * PLooseMarg << endl;
        outputLogFile << "level " << i + 1 << " pearson correlation critical value under loose significant threshold : " << 
            rCriticalValue[i] / (sampleSize - 1) << endl << endl;
    }
    distLvP.push_back(globalP); // set last distLvP element as global P threshold
    rCriticalValueCalc(min(distLvP[distLvNum] * PLooseMarg, 1.0), sampleSize, rCriticalValue[distLvNum]); // set last rCriticalValue element as global r critical value
    outputLogFile << "global significant threshold: P-value <= " << 
        distLvP[distLvNum] << endl;
    outputLogFile << "global loose significant threshold: P-value <= " << 
        min(distLvP[distLvNum] * PLooseMarg, 1.0) << endl;
    outputLogFile << "global pearson correlation critical value under loose significant threshold : " << 
        rCriticalValue[distLvNum] / (sampleSize - 1) << endl;
    outputLogFile << endl;

    // header of result file for output file
    ofstream outputFile;
    outputFile.open(outputFileName);
    outputFile << "omics1\t" << "omics2\t" << "BETA\t" << "SE\t" << "T\t" << "P\t" << "NMISS\t" << "distance_level\n";
    outputFile.close();
    
    // input omics1 data
    double* omics1Data = (double*) mkl_malloc(sizeof(double) * omics1Num * sampleSize, 64);
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
    outputLogFile << "First omics data has been input, time used: " << time_end_whole - time_start_whole << " s" << endl;

    // input omics2 data
    double* omics2Data = (double*) mkl_malloc(sizeof(double) * omics2Num * sampleSize, 64);
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
    outputLogFile << "Second omics data has been input, time used: " << time_end_whole - time_start_whole << " s" << endl;
    
    // input covariates data
    double* covarData;
    if (covarNum > 0) {
        covarData = inputCovar(covarFileName, 
                               covarNum, sampleSize, 
                               sampleFltSign, covarNANum, 
                               categFlag, covarCategNum);
    }
    vector<vector<uint32_t> > NASignMarkC(covarNum, vector<uint32_t>(0)); // NA mark for covarates. in fact it is empty

    // input replication list
    vector<pair<uint64_t, uint64_t> > rplList;
    bool rplFlag = false;
    if (rplFileName != "NA") {
        inputRplList(rplFileName, rplList, omics1Name, omics2Name, threadMaxN);
        rplFlag = true;
        outputLogFile << rplList.size() << " pairs of xQTL will be tested." << endl << endl;
    } else {
        outputLogFile << omics1Num * omics2Num << " pairs of xQTL will be tested." << endl << endl;
    }

    // orthogonal projection start

    // first centralization of SNP data
    vector<double> omics1RowSDCntrl1(omics1Num, 1);
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
    vector<double> omics2RowSDCntrl1(omics2Num, 1);
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics2Num, \
           omics2Data, \
           sampleSize, \
           omics2RowSDCntrl1, NASignMark2)
    for (uint64_t i = 0; i < omics2Num; i++) {
        cntrl(omics2Data, i, sampleSize, omics2RowSDCntrl1, NASignMark2);
    }
    
    // generate orthogonal omics data and intermediate variables
    double* omics1DataOrtg = (double*) mkl_malloc(sizeof(double) * omics1Num * sampleSize, 64); 
    memcpy(omics1DataOrtg, omics1Data, sizeof(double) * omics1Num * sampleSize);
    double* omics2DataOrtg = (double*) mkl_malloc(sizeof(double) * omics2Num * sampleSize, 64); 
    memcpy(omics2DataOrtg, omics2Data, sizeof(double) * omics2Num * sampleSize);
    if (covarNum > 0) {
        // first centralization of covariates
        vector<double> covarRowSDCntrl1(covarNum);
#       pragma omp parallel for \
        num_threads(threadMaxN) \
        shared(covarNum, \
               covarData, \
               sampleSize, \
               covarRowSDCntrl1)
        for (uint64_t i = 0; i < covarNum; i++) {
            cntrl(covarData, i, sampleSize, covarRowSDCntrl1, NASignMarkC);
        }

        // transpose covar
        double* covarDataT = (double*) mkl_malloc(sizeof(double) * sampleSize * covarNum, 64);
        for (uint64_t i = 0; i < covarNum; i++)
            for (uint64_t j = 0; j < sampleSize; j++)
                covarDataT[j * covarNum + i] = covarData[i * sampleSize + j];
        
        // QR decompose
        double *tau = new double[covarNum];
        LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, sampleSize, covarNum, covarDataT, covarNum, tau);
        LAPACKE_dorgqr(LAPACK_ROW_MAJOR, sampleSize, covarNum, covarNum, covarDataT, covarNum, tau);
        
        // projection
        double* covarNormSqr = (double*) mkl_malloc(sizeof(double) * sampleSize * sampleSize, 64);
        mkl_set_num_threads_local(threadMaxN); // Specifies the number of threads for MKL
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, sampleSize, sampleSize, covarNum, 
            1, covarDataT, covarNum, covarDataT, covarNum, 0, covarNormSqr, sampleSize);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, omics1Num, sampleSize, sampleSize, 
            -1, omics1Data, sampleSize, covarNormSqr, sampleSize, 1, omics1DataOrtg, sampleSize);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, omics2Num, sampleSize, sampleSize, 
            -1, omics2Data, sampleSize, covarNormSqr, sampleSize, 1, omics2DataOrtg, sampleSize);
        mkl_set_num_threads_local(0); // reset the thread-local number to global number

        // free intermediate variables
        mkl_free(covarNormSqr);
    }

    // second centralization for origin omics data
    // NormMod = 0 : not process/1 : cntrl/2 : cntrlQuant
    vector<double> omics1RowSDCntrl2(omics1Num, 1);
    vector<double> omics2RowSDCntrl2(omics2Num, 1);
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics1Num, omics1Data, sampleSize, omics1RowSDCntrl2, NASignMark1)
    for (uint64_t i = 0; i < omics1Num; i++) {
        if (omics1NormMod == 1) {
            cntrl(omics1Data, i, sampleSize, omics1RowSDCntrl2, NASignMark1);
        } else if (omics1NormMod == 2) {
            cntrlQuant(omics1Data, i, sampleSize, NASignMark1);
        }
    }
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics2Num, omics2Data, sampleSize, omics2RowSDCntrl2, NASignMark2)
    for (uint64_t i = 0; i < omics2Num; i++) {
        if (omics2NormMod == 1) {
            cntrl(omics2Data, i, sampleSize, omics2RowSDCntrl2, NASignMark2);
        } else if (omics2NormMod == 2) {
            cntrlQuant(omics2Data, i, sampleSize, NASignMark2);
        }
    }

    // second centralization for orthognal omics data
    // NormMod = 0 : cntrl/1 : cntrl/2 : cntrlQuant
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics1Num, omics1DataOrtg, sampleSize, omics1RowSDCntrl2, NASignMark1)
    for (uint64_t i = 0; i < omics1Num; i++) {
        if (omics1NormMod != 2) {
            cntrl(omics1DataOrtg, i, sampleSize, omics1RowSDCntrl2, NASignMark1);
        } else {
            cntrlQuant(omics1DataOrtg, i, sampleSize, NASignMark1);
        }
    }
#   pragma omp parallel for \
    num_threads(threadMaxN) \
    shared(omics2Num, omics2DataOrtg, sampleSize, omics2RowSDCntrl2, NASignMark2)
    for (uint64_t i = 0; i < omics2Num; i++) {
        if (omics2NormMod != 2) {
            cntrl(omics2DataOrtg, i, sampleSize, omics2RowSDCntrl2, NASignMark2);
        } else {
            cntrlQuant(omics2DataOrtg, i, sampleSize, NASignMark2);
        }
    }
    // Scalings of beta coefficents for orthognal omics data
    // if user do not need to centralize the variates, we should restore the coeffecient to the original scale
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
    
    // preprocess for linear model solver
    vector<double> omics1Sum(omics1Num, 0);
    vector<double> omics1Sqr(omics1Num);
    vector<double> omics1OrtgSqrInv(omics1Num);
    for (uint64_t i = 0; i < omics1Num; i++) {
        for (uint64_t j = i * sampleSize; j < (i + 1) * sampleSize; j++) {
            omics1Sum[i] += omics1Data[j];
        }
        omics1Sqr[i] = cblas_ddot(sampleSize, &omics1Data[i * sampleSize], 1, &omics1Data[i * sampleSize], 1);
        omics1OrtgSqrInv[i] = 1 / cblas_ddot(sampleSize, &omics1DataOrtg[i * sampleSize], 1, &omics1DataOrtg[i * sampleSize], 1);
    }
    vector<double> omics2Sum(omics2Num, 0);
    for (uint64_t i = 0; i < omics2Num; i++) {
        for (uint64_t j = i * sampleSize; j < (i + 1) * sampleSize; j++) {
            omics2Sum[i] += omics2Data[j];
        }
    }

    vector<double> covarSum(covarNum, 0);
    vector<vector<double> > omics1DotCov;
    omics1DotCov.resize(omics1Num, vector<double>(covarNum));
    vector<vector<double> > omics2DotCov;
    omics2DotCov.resize(omics2Num, vector<double>(covarNum));
    vector<vector<double> > CovarInter;
    CovarInter.resize(covarNum, vector<double>(covarNum));
    for (uint32_t i = 0; i < covarNum; i++) { CovarInter[i].resize(i + 1); }
    for (uint32_t i = 0; i < covarNum; i++) {
        for (uint64_t j = i * sampleSize; j < (i + 1) * sampleSize; j++) {
            covarSum[i] += covarData[j];
        }
        for (uint64_t j = 0; j < omics1Num; j++) {
            omics1DotCov[j][i] = cblas_ddot(sampleSize, &omics1Data[j * sampleSize], 1, &covarData[i * sampleSize], 1);
        }
        for (uint64_t j = 0; j < omics2Num; j++) {
            omics2DotCov[j][i] = cblas_ddot(sampleSize, &omics2Data[j * sampleSize], 1, &covarData[i * sampleSize], 1);
        }
        for (uint32_t j = 0; j <= i; j++) {
            CovarInter[i][j] = cblas_ddot(sampleSize, &covarData[i * sampleSize], 1, &covarData[j * sampleSize], 1);
        }
    }

    // time stamp for orthogonal projection
    time_end_whole = omp_get_wtime();
    outputLogFile << "Orthogonal projection has finished, time used: " << time_end_whole - time_start_whole << " s" << endl;

    // orthogonal projection finished

    // count the number of test of each distance level
    vector<uint64_t> testCnt(distLvNum + 1, 0);

    // split Block schedule
    uint64_t omics1BlockStride, omics2BlockStride, omics1BlockNum, omics2BlockNum, blockNum, rplListLen;
    if (rplFlag) {
        rplListLen = rplList.size();
        blockNum = (rplListLen + blockSize - 1) / blockSize;
    } else {
        omics1BlockStride = min(omics1Num, blockSize);
        omics2BlockStride = min(omics2Num, blockSize);
        omics1BlockNum = (int)(omics1Num + omics1BlockStride - 1) / omics1BlockStride;
        omics2BlockNum = (int)(omics2Num + omics2BlockStride - 1) / omics2BlockStride;
        blockNum = omics1BlockNum * omics2BlockNum;
    }

    // time stamp for preprocessing
    time_end_whole = omp_get_wtime();
    outputLogFile << "Preprocessing has finished, time used: " << time_end_whole - time_start_whole << " s" << endl;

#   pragma omp parallel \
    num_threads(min((uint64_t)threadMaxN, blockNum)) \
    shared(omics1Name, omics2Name, \
           omics1CHR, omics2CHR, omics1BP, omics2BP, \
           omics1Num, omics2Num, covarNum, sampleSize, missingRateThd, \
           omics1DataOrtg, omics2DataOrtg, \
           omics1BlockStride, omics2BlockStride, \
           distLv, distLvP, distLvNum, \
           rCriticalValue, \
           outputFileName, \
           outputLogFile, \
           NASignMark1, NASignMark2, \
           omics1Sum, omics2Sum, covarSum, \
           omics1Sqr, omics1OrtgSqrInv, \
           omics1DotCov, omics2DotCov, CovarInter, \
           omics1Scaling, omics2Scaling, \
           omics1RowSDCntrl1, omics1RowSDCntrl2, \
           testCnt)
    {
        // mark of omics id in current block
        uint64_t omics1BlockHead, omics2BlockHead, omics1BlockStrideCurr, omics2BlockStrideCurr;
        uint64_t pairAmt;

        // alloc space for current result on each thread
        vector<linearFitRlt> rltArr; rltArr.reserve(omics1BlockStride * omics2BlockStride / 4); // statistics results of significant variates pairs
        double* corr; if (!rplFlag) { corr = (double*) mkl_malloc(sizeof(double) * omics1BlockStride * omics2BlockStride, 64); }
        vector<uint64_t> testCntCurr(distLvNum + 1);

#       pragma omp for schedule(dynamic)
        for (uint32_t i = 0; i < blockNum; i++) {
            if (rplFlag) {
                pairAmt = min(blockSize, rplListLen - i * blockSize);
            } else { 
                omics1BlockHead = (int) i / omics2BlockNum * omics1BlockStride;
                omics2BlockHead = (int) i % omics2BlockNum * omics2BlockStride;
                omics1BlockStrideCurr = min(omics1BlockStride, omics1Num - omics1BlockHead);
                omics2BlockStrideCurr = min(omics2BlockStride, omics2Num - omics2BlockHead);
                pairAmt = omics1BlockStrideCurr * omics2BlockStrideCurr;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                            omics1BlockStrideCurr, omics2BlockStrideCurr, sampleSize, 
                            1, &omics1DataOrtg[omics1BlockHead * sampleSize], sampleSize, &omics2DataOrtg[omics2BlockHead * sampleSize], sampleSize, 
                            0, corr, min(omics2BlockStride, omics2Num - omics2BlockHead));
            }

            for (uint32_t it = 0; it < distLvNum + 1; it++) { testCntCurr[it] = 0; }
            uint64_t corrP = 0; double corrCurr; uint64_t j, k;
            for (uint64_t it = 0; it < pairAmt; it ++) {
                if (rplFlag) {
                    j = rplList[i * blockSize + it].first;
                    k = rplList[i * blockSize + it].second;
                    if (j < 0 || k < 0) { continue; }
                    corrCurr = cblas_ddot(sampleSize, &omics1DataOrtg[j * sampleSize], 1, &omics2DataOrtg[k * sampleSize], 1);
                } else {
                    j = it / omics2BlockStrideCurr + omics1BlockHead;
                    k = it % omics2BlockStrideCurr + omics2BlockHead;
                    corrCurr = corr[corrP++];
                }
                int32_t levelCurr = -1;
                if (distLvNum >= 1 &&
                    omics1CHR[j] == omics2CHR[k] && 
                    omics1CHR[j] > 0 && omics2CHR[k] > 0 && 
                    abs(omics1BP[j] - omics2BP[k]) <= distLv[distLvNum - 1] && 
                    omics1BP[j] > 0 && omics2BP[k] > 0) { // two locus are localed in the broadest level distance
                    for (uint32_t l = 0; l < distLvNum; l++) {
                        if (abs(omics1BP[j] - omics2BP[k]) <= distLv[l]) {
                            testCntCurr[l] += 1;
                            if (abs(corrCurr) >= rCriticalValue[l]) { levelCurr = l; }
                            break;
                        }
                    }
                } else { // two locus are not localed in any level distance
                    testCntCurr[distLvNum] += 1;
                    if (abs(corrCurr) >= rCriticalValue[distLvNum]) {
                        levelCurr = distLvNum;
                    }
                }

                if (levelCurr >= 0) {
                    linearFitRlt rltTmp = 
                    linearFit(corrCurr, 
                              j, k, levelCurr, 
                              sampleSize, covarNum, 
                              missingRateThd,
                              omics1Data, omics2Data, covarData,
                              NASignMark1, NASignMark2, 
                              omics1Sum, omics2Sum, covarSum, 
                              omics1Sqr, omics1OrtgSqrInv, 
                              omics1DotCov, omics2DotCov, CovarInter,
                              omics1Scaling, omics2Scaling, 
                              omics1RowSDCntrl1, omics2RowSDCntrl1);
                    if (rltTmp.p <= distLvP[levelCurr] && rltTmp.status == 0) {
                        rltArr.push_back(rltTmp);
                    }
                }                
            }

            // P< P1, maker BETA SE T P
            // P1 < P < P2, maker BETA T P
            // P > P2, do not output
#           pragma omp critical
            {
                ofstream outputFile;
                outputFile.open(outputFileName, ios::app);
                outputFile << setprecision(2);
                for (auto rltArrIt : rltArr){
                        outputFile << omics1Name[rltArrIt.omics1Id] << "\t";
                        outputFile << omics2Name[rltArrIt.omics2Id] << "\t";
                        outputFile << rltArrIt.b << "\t";
                        outputFile << rltArrIt.se << "\t";
                        outputFile << rltArrIt.t << "\t";
                        outputFile << rltArrIt.p << "\t";
                        outputFile << rltArrIt.nmiss << "\t";
                        outputFile << rltArrIt.level << "\n";
                }
                outputFile.close();
                for (uint32_t it = 0; it < distLvNum + 1; it++) { testCnt[it] += testCntCurr[it]; }
            }

            double time_end = omp_get_wtime();
            outputLogFile << i + 1 << " / " << blockNum << " block has finished, time used: " << time_end - time_start_whole << " s" << endl;
            
            // clear elements in vectors
            rltArr.clear();
        }

        // free vectors
        vector<linearFitRlt>(rltArr).swap(rltArr);
    }

    // whole calc time stamp
    time_end_whole = omp_get_wtime();
    outputLogFile << "Whole script time used: " << time_end_whole - time_start_whole << " s" << endl << endl;
    for (uint32_t it = 0; it < distLvNum + 1; it++) { 
        outputLogFile << "Number of test of " << it + 1 << " distance level : " << testCnt[it] << endl;
    }
    outputLogFile << endl;
    outputLogFile.close();

    mkl_free(omics1Data); mkl_free(omics2Data);
    if (covarNum > 0) {
        mkl_free(omics1DataOrtg); mkl_free(omics2DataOrtg);
    }

    return 0;
}
