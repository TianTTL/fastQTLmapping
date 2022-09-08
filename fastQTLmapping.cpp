#include <iostream>
#include <iomanip>
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
#include <gsl/gsl_statistics_double.h>
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
using std::copy;

namespace meqtllib {
// calc PLINK BED file size
void calcBfileSize(string bfileNameRoot, uint32_t &sampleSize, uint32_t &omicsSize) {
    string s;
    ifstream inputFile;

    inputFile.open(bfileNameRoot + ".fam", ios::in);
    for (sampleSize = 0; getline(inputFile, s); ++sampleSize)
    ;
    inputFile.close();

    inputFile.open(bfileNameRoot + ".bim", ios::in);
    for (omicsSize = 0; getline(inputFile, s); ++omicsSize)
    ;
    inputFile.close();
}

//input loci information of BED file
void getBfileSNPid(string bfileNameRoot, uint32_t num_snps, 
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

// input BED file
void inputOmicsBed(istreambuf_iterator<char>& inputFile,
                   float *omicsData, 
                   uint32_t locusCount, uint32_t sampleSize, vector<bool>& sampleFltSign, uint32_t sampleFltNum, 
                   vector<vector<uint32_t> >& NASignMark) {
#   pragma omp single
    {
        std::string s, one_item;

        // flushing NASignMark vectors
        for (uint32_t i = 0; i < locusCount; i++) {
            NASignMark[i].clear();
        }

        // variates for unpackGeno
        float *snp_geno_d = new float[sampleSize + sampleFltNum];
        float *snp_mask_d = new float[sampleSize + sampleFltNum];
        const float geno_table[4] = {2.0, 0.0, 1.0, 0.0}; // Homozygote A1, missing, Heterozygote, Homozygote A2
        const float mask_table[4] = {0.0, 1.0, 0.0, 0.0};
        int num_full_bytes_ = sampleSize / 4;
        int num_samples_left_ = sampleSize % 4;
        uint8_t t;

        // total element count
        uint32_t elementCount = 0;

        for (int i = 0; i < locusCount; i++) {
            for (size_t j = 0; j < num_full_bytes_; ++j) {
                t = *inputFile++;
                snp_geno_d[4 * j] = geno_table[t & 3u];
                snp_mask_d[4 * j] = mask_table[t & 3u];
                t >>= 2;
                snp_geno_d[4 * j + 1] = geno_table[t & 3u];
                snp_mask_d[4 * j + 1] = mask_table[t & 3u];
                t >>= 2;
                snp_geno_d[4 * j + 2] = geno_table[t & 3u];
                snp_mask_d[4 * j + 2] = mask_table[t & 3u];
                t >>= 2;
                snp_geno_d[4 * j + 3] = geno_table[t & 3u];
                snp_mask_d[4 * j + 3] = mask_table[t & 3u];
            }
            if (num_samples_left_ > 0u) {
                t = *inputFile++;
                for (size_t j = 0; j < num_samples_left_; ++j) {
                    snp_geno_d[4 * num_full_bytes_ + j] = geno_table[t & 3u];
                    snp_mask_d[4 * num_full_bytes_ + j] = mask_table[t & 3u];
                    t >>= 2;
                }
            }

            uint32_t sid = 0;
            for (int j = 0; j < sampleSize + sampleFltNum; j++) {
                if (!sampleFltSign[j]) {
                    if (snp_mask_d[j] > 0) {
                        omicsData[elementCount] = 0.0; // fill NA value to 0
                        NASignMark[i].push_back(sid);
                    } else {
                        omicsData[elementCount] = snp_geno_d[j];
                    }
                    elementCount++;
                    sid++;
                }
            }
        }
    }
}

// calc 2 dimension omics file size
void calcInputSize(string omicsFileName, uint32_t &sampleSize, uint32_t& omicsNum) {
    uint32_t i, rowsCount = 0, colsCount = 0;
    ifstream inputFile;
    string s, oneItem;

    inputFile.open(omicsFileName);

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

// input loci information of 2 dimension omics file
void get2DfloatId(string fileName, uint32_t omicsNum, 
                  vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP) {
    ifstream inputFile;
    string one_line, oneItem;
    string delimiter = " \t";
    string::size_type pos,lastPos;
    inputFile.open(fileName);
    for (uint32_t i = 0; i < omicsNum; i++) {
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
    }
    inputFile.close();
}

// input 2 dimension omics file
void input2DfloatParse(std::ifstream& inputFile, 
                       float* omicsData, 
                       uint32_t locusCount, uint32_t sampleSize, 
                       string* dataArea, 
                       uint32_t threadMaxN, 
                       vector<bool>& sampleFltSign, uint32_t sampleFltNum, 
                       vector<vector<uint32_t> >& NASignMark, string NASign) {
#   pragma omp single
    {
        string one_line, oneItem;
        string delimiter = " \t";
        string::size_type pos,lastPos;

        // flushing NASignMark vectors
        for (uint32_t i = 0; i < locusCount; i++) {
            NASignMark[i].clear();
        }

        for (uint32_t i = 0; i < locusCount; i++) {
            getline(inputFile, one_line);

            lastPos = one_line.find_first_not_of(delimiter, 0);
            pos = one_line.find_first_of(delimiter, lastPos); // SNP
            
            lastPos = one_line.find_first_not_of(delimiter, pos);
            pos = one_line.find_first_of(delimiter, lastPos); // CHR
            
            lastPos = one_line.find_first_not_of(delimiter, pos);
            pos = one_line.find_first_of(delimiter, lastPos); // BP
            
            dataArea[i] = one_line.substr(pos, one_line.length() - pos);
        }
    }

    // input omics data
#   pragma omp for schedule(dynamic) 
    for (uint32_t i = 0; i < locusCount; i++) {
        uint32_t lineLength = dataArea[i].length();
        char s[lineLength]; strcpy(s, dataArea[i].c_str());
        uint32_t s_p = 0;
        uint64_t rowHeadPos = i * sampleSize;
        uint32_t sid = 0;

        for (uint32_t j = 0; j < sampleSize + sampleFltNum; j++) {
            while (s[s_p] == ' ' || s[s_p] == '\t') s_p++; // skip space
            if (s_p >= lineLength) {
                cout << "column number lack\n";
                exit(1); // column number lack
            }

            if (!strncmp(&s[s_p], NASign.data(), NASign.length())) { // missing sign
                if (!sampleFltSign[j]) {
                    omicsData[rowHeadPos + sid] = 0.0;  // fill NA value to 0.0
                    NASignMark[i].push_back(sid);
                    sid++;
                }
                s_p += NASign.length();
            } else if (s[s_p] == '-' || s[s_p] >= '0' && s[s_p] <= '9') {
                int32_t sb = 1;
                if (s[s_p] == '-') {
                    sb = -1; s_p++;
                } // symbol

                float realData = 0;
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

void calcCovarSize(string covarFileName, string NASign, uint32_t sampleSize, uint32_t& covarNum, 
                   vector<bool>& sampleFltSign, uint32_t& covarNANum, 
                   vector<uint32_t>& categFlag, uint32_t& covarCategNum) {
    uint32_t i, rowsCount = 0;
    covarCategNum = 0;
    ifstream inputFile;
    string s, oneItem;
    vector<uint32_t> categFlagRef;

    inputFile.open(covarFileName);
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
float* inputCovar(string fileName, 
                  uint32_t& covarNum, uint32_t sampleSize, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum, 
                  vector<uint32_t>& categFlag, uint32_t covarCategNum) {
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
    float* covarData = (float*) mkl_malloc(sizeof(float) * 
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

// binary search for an ascending iterator
// Output the position of the first element greater than val
template <class ForwardIterator, class T>
int64_t binarySearch(ForwardIterator head, ForwardIterator tail, const T& val) {
    ForwardIterator it, headBak = head;
    typename std::iterator_traits<ForwardIterator>::difference_type count, step;
    count = distance(head, tail);
    while (count > 0)
    {
        it = head; step = count/2; advance(it,step);
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

// binary search for an decending iterator
// Output the position of the last element more than val
template <class ForwardIterator, class T>
int64_t binarySearchDec(ForwardIterator head, ForwardIterator tail, const T& val) {
    ForwardIterator it, headBak = head;
    typename std::iterator_traits<ForwardIterator>::difference_type count, step;
    count = distance(head, tail);
    while (count > 0)
    {
        it = head; step = count/2; advance(it,step);
        if (*it > val) {
            head = ++it;
            count -= step+1;
        }
        else count = step;
    }
    if (!(head == tail) && !(val < *head)) {
        return(distance(headBak, head) - 1);
    } else {
        return(-1);
    }
}

void inputRplList(string rplFileName, vector<pair<int64_t, int64_t> >& rplList, 
                  vector<string>& omics1Name, vector<string>& omics2Name, uint32_t threadMaxN) {
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
        vector<pair<string, uint32_t> > omics1NameSortTmp(omics1Name.size()); // combine omics id and order
        for (uint32_t i = 0; i < omics1Name.size(); ++i) { omics1NameSortTmp[i] = make_pair(omics1Name[i], i); }
        sort(omics1NameSortTmp.begin(), omics1NameSortTmp.end());
        vector<string> omics1NameSort(omics1Name.size()); vector<uint32_t> omics1NameRank(omics1Name.size()); // split sorted omics id and rank
        for (uint32_t i = 0; i < omics1Name.size(); ++i) { omics1NameSort[i] = omics1NameSortTmp[i].first; }
        for (uint32_t i = 0; i < omics1Name.size(); ++i) { omics1NameRank[i] = omics1NameSortTmp[i].second; }

        vector<pair<string, uint32_t> > omics2NameSortTmp(omics2Name.size()); // combine omics id and order
        for (uint32_t i = 0; i < omics2Name.size(); ++i) { omics2NameSortTmp[i] = make_pair(omics2Name[i], i); }
        sort(omics2NameSortTmp.begin(), omics2NameSortTmp.end());
        vector<string> omics2NameSort(omics2Name.size()); vector<uint32_t> omics2NameRank(omics2Name.size()); // split sorted omics id and rank
        for (uint32_t i = 0; i < omics2Name.size(); ++i) { omics2NameSort[i] = omics2NameSortTmp[i].first; }
        for (uint32_t i = 0; i < omics2Name.size(); ++i) { omics2NameRank[i] = omics2NameSortTmp[i].second; }

        // build replication list by binary searching
#       pragma omp parallel for schedule(static)
        for (uint32_t i = 0; i < rplFile.size(); i++) {
            pair<int64_t, int64_t> oneItem;
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
template <typename T>
void cntrl(T* omicsData, uint32_t omicsId, 
           uint32_t sampleSize, 
           vector<float>& rowSD,
           float sdThd, 
           vector<vector<uint32_t> >& NASignMarkCurr) {
    float rowSumTmp, rowSDTmp, rowMeanTmp;
    uint64_t rowHeadPos = omicsId * sampleSize;

    uint32_t sampleSizeCurr = sampleSize - NASignMarkCurr[omicsId].size();

    // calc the row mean
    rowSumTmp = 0;
    for (uint32_t i = 0; i < sampleSize; i++) {
        rowSumTmp += omicsData[rowHeadPos + i];
    }
    rowMeanTmp = rowSumTmp / sampleSizeCurr;

    // calc the row SD for omics
    rowSDTmp = 0;
    for (uint32_t i = 0; i < sampleSize; i++){
        if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), i) == NASignMarkCurr[omicsId].end()) {
            rowSDTmp += pow(omicsData[rowHeadPos + i] - rowMeanTmp, 2);
        }
    }
    rowSD[omicsId] = sqrt(rowSDTmp / (sampleSizeCurr - 1));

    // centralization
    if (rowSD[omicsId] <= sdThd) { // constant variants
        for (uint32_t i = 0; i < sampleSize; i++){
            omicsData[rowHeadPos + i] = 0;
        }
    } else {
        for (uint32_t i = 0; i < sampleSize; i++){
            if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), i) == NASignMarkCurr[omicsId].end()) {
                omicsData[rowHeadPos + i] = (omicsData[rowHeadPos + i] - rowMeanTmp) / rowSD[omicsId];
            }
        }
    }
}

// quantile based normalization
template <typename T>
void cntrlQuant(T *omicsData, uint32_t omicsId, 
                uint32_t sampleSize,
                vector<float>& rowSD, 
                float sdThd, 
                vector<vector<uint32_t> >& NASignMarkCurr) {
    float rowSumTmp, rowSDTmp, rowMeanTmp;
    vector<T> v_temp(sampleSize);
    vector<float> v_rank(sampleSize);
    uint64_t rowHeadPos = omicsId * sampleSize;
    T maxV;

    uint32_t sampleSizeCurr = sampleSize - NASignMarkCurr[omicsId].size();

    // calc the row mean
    rowSumTmp = 0;
    for (uint32_t i = 0; i < sampleSize; i++) {
        rowSumTmp += omicsData[rowHeadPos + i];
    }
    rowMeanTmp = rowSumTmp / sampleSizeCurr;

    // calc the row SD for omics
    rowSDTmp = 0;
    for (uint32_t i = 0; i < sampleSize; i++){
        if (find(NASignMarkCurr[omicsId].begin(), NASignMarkCurr[omicsId].end(), i) == NASignMarkCurr[omicsId].end()) {
            rowSDTmp += pow(omicsData[rowHeadPos + i] - rowMeanTmp, 2);
        }
    }
    rowSD[omicsId] = sqrt(rowSDTmp / (sampleSizeCurr - 1));

    if (rowSD[omicsId] <= sdThd) { // constant variants
        for (uint32_t i = 0; i < sampleSize; i++){
            omicsData[rowHeadPos + i] = 0;
        }
        return;
    }

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
        omicsData[rowHeadPos + i] = gsl_cdf_ugaussian_Pinv((v_rank[i] + 0.5) / sampleSizeCurr); // caution: v_rank start from 0!
    }

    // fill the missing value
    for (auto i : NASignMarkCurr[omicsId]) {
        omicsData[rowHeadPos + i] = 0;
    }

    return;
}

// sort and rank a vector
template <typename T>
vector<float> rankSort(const vector<T>& v_temp, uint32_t sampleSize) {
    vector<pair<T, uint32_t> > v_sort(sampleSize);

    for (uint32_t i = 0; i < v_sort.size(); ++i) {
        v_sort[i] = make_pair(v_temp[i], i);
    }

    sort(v_sort.begin(), v_sort.end());

    vector<float> result(sampleSize);
    uint32_t currentRankP = 0, preRankP = -1;
    T currentValue;
    float currentRank;

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

// comparator of results of linear fitting
bool operator<(const fitRltPart &A, const fitRltPart &B) {
    if (A.level != B.level) {
        return(A.level < B.level);
    } else {
        return(A.p < B.p);
    }
}

fitRlt linearFit(float corr, 
                 uint32_t omics1Id, uint32_t omics2Id, 
                 uint32_t omics1GlobalId, uint32_t omics2GlobalId, uint8_t level, 
                 uint32_t sampleSize, uint32_t covarNum, 
                 float msRtThd,
                 float* omics1Data, float* omics2Data, float* covarData, 
                 vector<vector<uint32_t> >& NASignMark1, vector<vector<uint32_t> >& NASignMark2, 
                 vector<float>& omics1Sum, vector<float>& omics2Sum, vector<float>& covarSum, 
                 vector<float>& omics1Sqr, vector<float>& omics1OrtgSqrInv, 
                 vector<vector<float> >& omics1DotCov, vector<vector<float> >& omics2DotCov, vector<vector<float> >& CovarInter,
                 int32_t omics1NormMod, int32_t omics2NormMod, 
                 vector<float>& omics1Scaling, vector<float>& omics2Scaling, 
                 vector<float>& omics1RowSDCntrl1, vector<float>& omics2RowSDCntrl1) {
    uint32_t df_r, df_t, sampleSizeCurr;
    vector<uint32_t> NASignMarkCurr; NASignMarkCurr.reserve(sampleSize);
    fitRlt rlt;
    rlt.omics1Id = omics1GlobalId; rlt.omics2Id = omics2GlobalId; // locus index start from 0
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
    if (NASignMarkCurr.size() > msRtThd * sampleSize) {
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
        b = omics1OrtgSqrInv[omics1Id] * corr * df_r;
        // fix data scale
        if (omics1NormMod == 0) {
            b = b / omics1Scaling[omics1Id];
        }
        if (omics2NormMod == 0) {
            b = b * omics2Scaling[omics2Id];
        }
        se = b / t;
    } else {
        // Gauss-Jordan algorithm
        double* omics1DataCurr = (double*) mkl_malloc(sizeof(double) * sampleSize, 64);
        double* omics2DataCurr = (double*) mkl_malloc(sizeof(double) * sampleSize, 64);
        copy(&omics1Data[omics1Id * sampleSize], &omics1Data[omics1Id * sampleSize + sampleSize], omics1DataCurr);
        copy(&omics2Data[omics2Id * sampleSize], &omics2Data[omics2Id * sampleSize + sampleSize], omics2DataCurr);
        
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
        copy(omics1DataCurr, omics1DataCurr + sampleSize, &ivMt[0 + sampleSize]);
        copy(covarData, covarData + covarNum * sampleSize, &ivMt[0 + 2*sampleSize]);
        double* resid = (double*) mkl_malloc(sizeof(double) * sampleSize, 64); // residuals vector
        copy(omics2DataCurr, omics2DataCurr + sampleSize, resid);
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
        if (omics1NormMod == 0) {
            b = b / omics1RowSDCntrl1[omics1Id];
        }
        if (omics2NormMod == 0) {
            b = b * omics2RowSDCntrl1[omics2Id];
        }
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

void rCriticalValueCalc(double P, uint32_t sampleSize, uint32_t covarNum, double &rCriticalValue) {
    double tCriticalValue;
    tCriticalValue = gsl_cdf_tdist_Qinv(P/2, sampleSize - (1 + covarNum) - 1); // t critical value
    rCriticalValue = sqrt(pow(tCriticalValue, 2) / (sampleSize - (1 + covarNum) - 1 + pow(tCriticalValue, 2))) * (sampleSize - 1); // correlation critical value, times sampleSize
}

// global variables
extern string omics1FileName, omics2FileName, outputFileName, rplFileName;
extern string covarFileName;
extern bool bfileFlag1, bfileFlag2;
extern string NASign;
extern vector<uint32_t> categFlag;

extern double globalP;
extern float FWER;
extern float msRtThd;

extern vector<float> distLv;
extern vector<double> distLvP;
extern uint8_t distLvNum;

extern uint32_t threadMaxN;
extern int32_t omics1NormMod, omics2NormMod;
extern float PLooseMarg;
extern uint32_t outPcs;
extern int32_t helpFlag;
extern int32_t modeFlag;
extern uint32_t chunkSize;

extern uint32_t omics1Num, omics2Num, covarNum, sampleSize;
extern ofstream outputLogFile;
extern ostringstream oss;

}  // namespace meqtllib

void dualOutput(std::ostringstream &oss, std::ostream &os1, std::ostream &os2) {
    os1 << oss.str();
    os2 << oss.str();
    oss.str(""); oss.clear(); 
}

int cntModeProc();
int discModeProc();

int main(int argc, char *argv[]) {
    // clipp
    auto firstOpt = "mode-independent options:" % (
        required("--omics1") & value("omics1FileName", omics1FileName)       % "first omics file path"
            & option("bfile").set(bfileFlag1)     % "using plink binary format represente first omics",
        required("--omics2") & value("omics2FileName", omics2FileName)      % "second omics file path"
            & option("bfile").set(bfileFlag2)    % "using plink binary format represente second omics",
        required("--out") & value("outputFileName", outputFileName)               % "output file path",
        option("--outPcs") & integer("outPcs", outPcs)            % "the number of significant digits",
        option("--threads") & integer("threadMaxN", threadMaxN)                       % "max. threads",
        option("--chunk") & value("chunkSize", chunkSize)             % "dimension of splitting chunk",
        option("-h", "--help").set(helpFlag, 1)                                          % "show help"
    );

    auto cntMode = "loci-pairs counting mode:" % (
        command("count").set(modeFlag, 1)                                             % "mode command",
        option("--dl") & numbers("distLv", distLv)     % "distance thresholds for each distance level", 
        option("--FWER") & value("FWER", FWER)                             % "Family-wise Error Error"
    );

    auto discMode = "discovery mode:" % (
        command("discovery").set(modeFlag, 2)                                         % "mode command",
        option("-p") & number("globalP", globalP)                         % "global P-value threshold",
        option("--ploose") & number("PLooseMarg", PLooseMarg) %   "margin for loose P-value threshold",
        option("--cov") & value("covarFileName", covarFileName)               % "covariates file path",
        option("--categ") & integers("categFlag", categFlag)  % "flag indicate categorical covariates",
        option("--na") & value("NASign", NASign)                             % "sign of missing value",
        option("--MR") & number("msRtThd", msRtThd)               % "missing rate threshold",
        option("--SD") & number("sdThd", sdThd)       % "standard deviation threshold",
        option("--dl") & numbers("distLv", distLv)     % "distance thresholds for each distance level",
        option("--dlp") & numbers("distLvP", distLvP)   % "P-value thresholds for each distance level",
        option("--omics1norm")                                 % "normalization model for omics1 data"
            & (required("zscore").set(omics1NormMod, 1)
                                 | required("rank").set(omics1NormMod, 2)),
        option("--omics2norm")                                 % "normalization model for omics2 data"
            & (required("zscore").set(omics2NormMod, 1)
                                 | required("rank").set(omics2NormMod, 2))
    );

    auto rplMode = "replication mode:" % (
        command("replication").set(modeFlag, 3)                                       % "mode command",
        option("--cov") & value("covarFileName", covarFileName)               % "covariates file path",
        option("--categ") & integers("categFlag", categFlag)  % "flag indicate categorical covariates",
        option("--na") & value("NASign", NASign)                             % "sign of missing value",
        option("--MR") & number("msRtThd", msRtThd)               % "missing rate threshold",
        option("--SD") & number("sdThd", sdThd)       % "standard deviation threshold",
        option("--dl") & numbers("distLv", distLv)     % "distance thresholds for each distance level",
        option("--omics1norm")                                 % "normalization model for omics1 data"
            & (required("zscore").set(omics1NormMod, 1)
                                 | required("rank").set(omics1NormMod, 2)),
        option("--omics2norm")                                 % "normalization model for omics2 data"
            & (required("zscore").set(omics2NormMod, 1)
                                 | required("rank").set(omics2NormMod, 2)),
         option("--rpl") & value("rplFileName", rplFileName)             % "replication list file path"
    );

    auto cli = (
        firstOpt,
        cntMode | discMode //| rplMode
    );

    // output man
    auto fmt = doc_formatting{} .first_column(4) .doc_column(25) .last_column(80);
    if(!parse(argc, argv, cli) || helpFlag == 1) {
        if (helpFlag != 1) {
            cout << endl << "Error: Unknown option." << endl << endl;    
        }
        cout << make_man_page(cli, "fastQTLmapping", fmt)
        .prepend_section("DESCRIPTION", "    Fastest QTL mapping tool. Version " + VERSION)
        .append_section("LICENSE", "    GPL3") << '\n';
        return 0;
    }

    // initializing log file
    string outputLogFileName = outputFileName + ".log";
    outputLogFile.open(outputLogFileName);
    outputLogFile << setprecision(outPcs);
    oss << "\n######################################\n" 
        << "  fastQTLmapping version " << VERSION << " start  \n" 
        << "######################################\n\n"; 
    dualOutput(oss, outputLogFile, std::cout);

    // check input file exist
    ifstream inputFileCk;
    if (bfileFlag1) {
        inputFileCk.open(omics1FileName + ".bed", ios::in | ios::binary);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics1 bed file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
        inputFileCk.open(omics1FileName + ".bim", ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics1 bim file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
        inputFileCk.open(omics1FileName + ".fam", ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics1 fam file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
    } else {
        inputFileCk.open(omics1FileName, ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics1 data file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
    }
    if (bfileFlag2) {
        inputFileCk.open(omics2FileName + ".bed", ios::in | ios::binary);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics1 bed file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
        inputFileCk.open(omics2FileName + ".bim", ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics1 bim file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
        inputFileCk.open(omics2FileName + ".fam", ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics2 fam file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
    } else {
        inputFileCk.open(omics2FileName, ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Omics2 data file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
    }
    if (!covarFileName.empty()) {
        inputFileCk.open(covarFileName, ios::in);
        if (!inputFileCk.is_open()) {
            oss << "Error: Covariates file not exist.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        inputFileCk.close();
    }

    // check global P
    if (globalP <= 0 || globalP > 1) {
        oss << "Error: P-value threshold should in (0, 1].\n"; dualOutput(oss, outputLogFile, std::cout);
        return 1;
    }

    // fill distLvP with globalP
    if (distLvP.size() == 0) {
        distLvP.assign(distLv.size(), globalP);
    }

    // check equal length about distLv and distLvP
    if (distLv.size() != distLvP.size()) {
        oss << "Error: Length of distLv and distLvP is not equal.\n"; dualOutput(oss, outputLogFile, std::cout);
        return 2;
    }
    distLvNum = distLv.size();

    // check ascending about distLv
    for (uint32_t i = 1; i < distLv.size(); i++) {
        if (distLv[i] <= distLv[i - 1]) {
            oss << "Error: DistLv is not ascending.\n"; dualOutput(oss, outputLogFile, std::cout);
            return 3;
        }
    }

    int FQMstat;
    if (modeFlag == 1) {
        FQMstat = cntModeProc();
    } else if (modeFlag == 2) {
        FQMstat = discModeProc();
    } else if (modeFlag == 3) {
        // FQMstat = rplModeProc();
    }

    // close log file
    outputLogFile.close();

    return(FQMstat);
}

int cntModeProc() {
    // global starting time stamp 
    double time_start_whole = omp_get_wtime(), time_end;

    // record omics file into log file
    oss << "omics 1 file : " << omics1FileName << endl;
    if (bfileFlag1) {
        oss << "  omics 1 is in PLINK binary format\n";
    }
    oss << "omics 2 file : " << omics2FileName << endl;
    if (bfileFlag2) {
        oss << "  omics 2 is in PLINK binary format\n";
    }
    oss << "output file : " << outputFileName << endl;
    oss << "maximun parallel number : " << threadMaxN << endl; dualOutput(oss, outputLogFile, std::cout);

    // calculate input file size
    if (bfileFlag1) { // input plink bfile as first omics
        calcBfileSize(omics1FileName, sampleSize, omics1Num);
    } else {
        calcInputSize(omics1FileName, sampleSize, omics1Num);
    }
    if (bfileFlag2) { // input plink bfile as second omics
        calcBfileSize(omics2FileName, sampleSize, omics2Num);
    } else {
        calcInputSize(omics2FileName, sampleSize, omics2Num);
    }
    
    // record data scale into log file
    oss << "omics 1 number : " << omics1Num << endl; dualOutput(oss, outputLogFile, std::cout);
    oss << "omics 2 number : " << omics2Num << endl; dualOutput(oss, outputLogFile, std::cout);
    oss << "sample number : " << sampleSize << endl << endl; dualOutput(oss, outputLogFile, std::cout);
    
    // input variates informatrion of omics1 data
    vector<string> omics1Name(omics1Num); // locus name for first omics
    vector<int32_t> omics1CHR(omics1Num); // CHR number for first omics
    vector<int64_t> omics1BP(omics1Num); // BP number for first omics
    if (bfileFlag1) { // input plink bfile as first omics
        getBfileSNPid(omics1FileName, omics1Num, omics1Name, omics1CHR, omics1BP);
    } else {
        get2DfloatId(omics1FileName, omics1Num, omics1Name, omics1CHR, omics1BP);
    }

    // input variates informatrion of omics2 data
    vector<string> omics2Name(omics2Num); // locus name for second omics
    vector<int32_t> omics2CHR(omics2Num); // CHR number for second omics
    vector<int64_t> omics2BP(omics2Num); // BP number for second omics
    if (bfileFlag2) { // input plink bfile as second omics
        getBfileSNPid(omics2FileName, omics2Num, omics2Name, omics2CHR, omics2BP);
    } else {
        get2DfloatId(omics2FileName, omics2Num, omics2Name, omics2CHR, omics2BP);
    }

    // count the number of test of each distance level
    vector<uint64_t> testCnt(distLvNum + 1, 0);
    int *testCntPrivate;
#   pragma omp parallel num_threads(threadMaxN)
    {
        const int nthreads = omp_get_num_threads();
        const int ithread = omp_get_thread_num();

#       pragma omp single
        {
            testCntPrivate = new int[(distLvNum + 1)*nthreads]();
        }

#       pragma omp for schedule(dynamic)
        for (uint32_t i = 0 ; i < omics1Num ; ++i)
        {
            uint32_t offset = (distLvNum + 1) * ithread;
            for (uint32_t j = 0 ; j < omics2Num ; ++j) {
                if (distLvNum >= 1 &&
                    omics1CHR[i] == omics2CHR[j] && 
                    omics1CHR[i] > 0 && omics2CHR[j] > 0 && 
                    abs(omics1BP[i] - omics2BP[j]) <= distLv[distLvNum - 1] && 
                    omics1BP[i] > 0 && omics2BP[j] > 0) { // two locus are localed in the broadest level distance
                    for (uint32_t k = 0; k < distLvNum; k++) {
                        if (abs(omics1BP[i] - omics2BP[j]) <= distLv[k]) {
                            testCntPrivate[offset + k] += 1;
                            break;
                        }
                    }
                } else {
                    testCntPrivate[offset + distLvNum] += 1;
                }
            }
        }
        
#       pragma omp for schedule(dynamic)
        for(uint8_t i = 0; i < (distLvNum + 1); i++) {
            uint32_t offset = i;
            for(uint32_t j = 0; j < nthreads; j++) {
                testCnt[i] += testCntPrivate[offset];
                offset += (distLvNum + 1);
            }
        }
    }
    delete[] testCntPrivate;

    // output counting results
    ofstream outputFile;
    outputFile.open(outputFileName + ".cnt");
    outputFile << setprecision(outPcs);
    for (uint8_t l = 0; l < distLvNum + 1; l++) { 
        outputFile << "Number of test of distance level " << l + 1 << " : \n\t" << testCnt[l] << endl;
        outputFile << "Bonferroni threshold of distance level " << l + 1 << " under FWER " << FWER << " : \n\t" << FWER / testCnt[l] << endl;
    }

    // global ending time stamp
    time_end = omp_get_wtime();
    oss << "Whole procedure time used : " << time_end - time_start_whole << " s" << endl << endl; dualOutput(oss, outputLogFile, std::cout);
}

int discModeProc() {
    // check PLooseMarg
    if (PLooseMarg < 1) {
        PLooseMarg = 1;
    }

    // global starting time stamp 
    double time_start_whole = omp_get_wtime(), time_end;

    // Specifies the global number of threads for MKL
    mkl_set_num_threads(1);

    // calculate input file size
    if (bfileFlag1) { // input plink bfile as first omics
        calcBfileSize(omics1FileName, sampleSize, omics1Num);
    } else {
        calcInputSize(omics1FileName, sampleSize, omics1Num);
    }
    if (bfileFlag2) { // input plink bfile as second omics
        calcBfileSize(omics2FileName, sampleSize, omics2Num);
    } else {
        calcInputSize(omics2FileName, sampleSize, omics2Num);
    }
    vector<bool> sampleFltSign(sampleSize, false);
    uint32_t covarNANum = 0;
    uint32_t covarCategNum = 0;
    if (!covarFileName.empty()) {
        calcCovarSize(covarFileName, NASign, sampleSize, covarNum, sampleFltSign, covarNANum, categFlag, covarCategNum);
    }
    sampleSize -= covarNANum;

    // input variates informatrion of omics1 data
    vector<string> omics1Name(omics1Num); // locus name for first omics
    vector<int32_t> omics1CHR(omics1Num); // CHR number for first omics
    vector<int64_t> omics1BP(omics1Num); // BP number for first omics
    if (bfileFlag1) { // input plink bfile as first omics
        getBfileSNPid(omics1FileName, omics1Num, omics1Name, omics1CHR, omics1BP);
    } else {
        get2DfloatId(omics1FileName, omics1Num, omics1Name, omics1CHR, omics1BP);
    }

    // input variates informatrion of omics2 data
    vector<string> omics2Name(omics2Num); // locus name for second omics
    vector<int32_t> omics2CHR(omics2Num); // CHR number for second omics
    vector<int64_t> omics2BP(omics2Num); // BP number for second omics
    if (bfileFlag2) { // input plink bfile as second omics
        getBfileSNPid(omics2FileName, omics2Num, omics2Name, omics2CHR, omics2BP);
    } else {
        get2DfloatId(omics2FileName, omics2Num, omics2Name, omics2CHR, omics2BP);
    }

    // record some parameters and data scale into log file 
    oss << "omics 1 file : " << omics1FileName << endl;
    if (bfileFlag1) {
        oss << "  omics 1 is in PLINK binary format\n";
    }
    oss << "  omics 1 variate number : " << omics1Num << endl;
    switch (omics1NormMod) {
        case 0 : 
            oss << "  not normalizing omics1 data\n";
            break;
        case 1 : 
            oss << "  normalizing omics1 data by z-score\n";
            break;
        case 2 : 
            oss << "  normalizing omics1 data by rank-based\n";
            break;
        default : 
            break;
    }
    oss << endl; dualOutput(oss, outputLogFile, std::cout);
    oss << "omics 2 file : " << omics2FileName << endl;
    if (bfileFlag2) {
        oss << "  omics 2 is in PLINK binary format\n";
    }
    oss << "  omics 2 variate number : " << omics2Num << endl;
    switch (omics2NormMod) {
        case 0 : 
            oss << "  not normalizing omics2 data\n";
            break;
        case 1 : 
            oss << "  normalizing omics2 data by z-score\n";
            break;
        case 2 : 
            oss << "  normalizing omics2 data by rank-based\n";
            break;
        default : 
            break;
    }
    oss << endl; dualOutput(oss, outputLogFile, std::cout);

    if (!covarFileName.empty()) {
        oss << "covar file : " << covarFileName << endl;
        oss << "covariates number : " << covarNum << endl;
        oss << "  numeric covariates number : " << covarNum - covarCategNum << endl;
        oss << "  categorical covariates number : " << covarCategNum << endl << endl;
        oss << "valid sample number : " << sampleSize << endl;
        oss << "  " << covarNANum << " samples are excluded because covariates missing" << endl;    
    } else {
        oss << "sample number : " << sampleSize << endl;
    }
    oss << endl; dualOutput(oss, outputLogFile, std::cout);

    oss << "output file : " << outputFileName << endl << endl;
    oss << "missing value sign : " << NASign << endl;
    oss << "missing rate threshold : " << msRtThd << endl << endl;
    oss << "number of OpenMP threads : " << threadMaxN << endl;
    oss << "dimension of splitting chunk : " << chunkSize << endl;
    oss << endl; dualOutput(oss, outputLogFile, std::cout);
    
    // critical value of t test
    vector<double> rCriticalValue(1 + distLvNum);
    for (uint8_t i = 0; i < distLvNum; i++) {
        rCriticalValueCalc(min(distLvP[i] * PLooseMarg, 1.0), sampleSize, covarNum, rCriticalValue[i]);
        oss << "level " << i + 1 << " distance threshold : " << 
            distLv[i] << endl;  dualOutput(oss, outputLogFile, std::cout); // level number start from 1
        oss << "level " << i + 1 << " significant threshold : P-value <= " << 
            distLvP[i] << endl; dualOutput(oss, outputLogFile, std::cout);
        oss << "  loose significant threshold : P-value <= " << 
            min(distLvP[i] * PLooseMarg, 1.0) << endl; dualOutput(oss, outputLogFile, std::cout);
        oss << "  pearson correlation critical value under loose significant threshold : " << 
            rCriticalValue[i] / (sampleSize - 1) << endl << endl; dualOutput(oss, outputLogFile, std::cout);
    }
    distLvP.push_back(globalP); // set last distLvP element as global P threshold
    rCriticalValueCalc(min(distLvP[distLvNum] * PLooseMarg, 1.0), sampleSize, covarNum, rCriticalValue[distLvNum]); // set last rCriticalValue element as global r critical value
    oss << "global significant threshold : P-value <= " << 
        distLvP[distLvNum] << endl; dualOutput(oss, outputLogFile, std::cout);
    oss << "  loose significant threshold : P-value <= " << 
        min(distLvP[distLvNum] * PLooseMarg, 1.0) << endl; dualOutput(oss, outputLogFile, std::cout);
    oss << "  pearson correlation critical value under loose significant threshold : " << 
        rCriticalValue[distLvNum] / (sampleSize - 1) << endl << endl; dualOutput(oss, outputLogFile, std::cout);

    // preprocess of ST FDR
    // Critical value of bins for each gradient significant level
    // seq(0.05, 0.95, 0.05)
    uint8_t st_ll = 19; // length of lambda
    vector<float> st_lambda(st_ll);
    // vector<double> st_lambdaRCV(st_ll + 2); // the first element is set to upper bound
    //                                         // the second element is set to lower bound 
    // ascending lambda
    generate(st_lambda.begin(), st_lambda.end(), [n = 0.00]() mutable { return n+=0.05; } );
    // decending lambda correlation critical value
    // for (uint32_t i = 0; i < st_ll; i++) {
    //     rCriticalValueCalc(st_lambda[i], sampleSize, covarNum, st_lambdaRCV[i + 1]);
    // }
    // st_lambdaRCV[0] = sampleSize; // upper boundary
    // st_lambdaRCV[st_ll+1] = 0; // lower boundary

    // header of result file for output file
    ofstream outputBinFile(outputFileName + ".bin", ios::out | ios::binary);
    
    // input covariates data
    float* covarData;
    if (covarNum > 0) {
        covarData = inputCovar(covarFileName, 
                               covarNum, sampleSize, 
                               sampleFltSign, covarNANum, 
                               categFlag, covarCategNum);
    }
    vector<vector<uint32_t> > NASignMarkC(covarNum, vector<uint32_t>(0)); // NA mark for covarates. in fact it is empty
    
    // Chunk spliting schedule
    uint32_t omics1ChunkStrideAllc, omics2ChunkStrideAllc, omics1ChunkNum, omics2ChunkNum, chunkNum;
    omics1ChunkStrideAllc = min(omics1Num, chunkSize * threadMaxN);
    omics2ChunkStrideAllc = min(omics2Num, chunkSize);
    omics1ChunkNum = (int)(omics1Num + omics1ChunkStrideAllc - 1) / omics1ChunkStrideAllc;
    omics2ChunkNum = (int)(omics2Num + omics2ChunkStrideAllc - 1) / omics2ChunkStrideAllc;
    chunkNum = omics1ChunkNum * omics2ChunkNum;

    // estimating peak memory
    uint64_t memCons = 
    max(
        (omics1ChunkStrideAllc + omics2ChunkStrideAllc) * sampleSize * (sizeof(float) * 2) + // omics data and orthgnal data
        omics1ChunkStrideAllc * omics2ChunkStrideAllc * sizeof(float) + // gemm
        omics1ChunkStrideAllc * omics2ChunkStrideAllc * globalP * PLooseMarg * 9.25 * sizeof(uint32_t), // tmp results
        omics1Num * omics2Num * globalP * 4.25 * sizeof(uint32_t) // all result of P and Q and level
    );
    if (memCons < 1024 * 1024 * 1024) {
        oss << "fastQTLmapping requires less than 1 GB RAM to complete the current task.\n\n"; dualOutput(oss, outputLogFile, std::cout);
    } else {
        oss << "fastQTLmapping requires about " << ceil(memCons * 1.0 / 1024 / 1024 / 1024) << " GB RAM to complete the current task.\n\n"; dualOutput(oss, outputLogFile, std::cout);
    }

    // clocking
    time_end = omp_get_wtime();
    oss << "Preprocessing has been completed, time used : " << time_end - time_start_whole << " s" << endl; dualOutput(oss, outputLogFile, std::cout);

    // open input data stream
    ifstream o1InputFile, o2InputFile;
    istreambuf_iterator<char> o1InputFile_it;
    istreambuf_iterator<char> o2InputFile_it;
    if (bfileFlag1) {
        o1InputFile.open(omics1FileName + ".bed", std::ios::in | std::ios::binary);
        o1InputFile_it = o1InputFile;
        o1InputFile_it++; o1InputFile_it++; o1InputFile_it++;
    } else {
        o1InputFile.open(omics1FileName);
    }

    if (bfileFlag2) {
        o2InputFile.open(omics2FileName + ".bed", std::ios::in | std::ios::binary);
        o2InputFile_it = o2InputFile;
        o2InputFile_it++; o2InputFile_it++; o2InputFile_it++;
    } else {
        o2InputFile.open(omics2FileName);
    }

    // assign input omics data
    float* omics1Data = (float*) mkl_malloc(sizeof(float) * omics1ChunkStrideAllc * sampleSize, 64);
    vector<vector<uint32_t> > NASignMark1(omics1ChunkStrideAllc, vector<uint32_t>(0)); // NA mark for first omics, N.O1 * NAs
    float* omics2Data = (float*) mkl_malloc(sizeof(float) * omics2ChunkStrideAllc * sampleSize, 64);
    vector<vector<uint32_t> > NASignMark2(omics2ChunkStrideAllc, vector<uint32_t>(0)); // NA mark for second omics, N.O2 * NAs
    string* dataArea = new string[max(omics1ChunkStrideAllc, omics2ChunkStrideAllc)];

    // assign precompute
    vector<float> omics1RowSDCntrl1(omics1ChunkStrideAllc, 1);
    vector<float> omics2RowSDCntrl1(omics2ChunkStrideAllc, 1);
    float* omics1DataOrtg = (float*) mkl_malloc(sizeof(float) * omics1ChunkStrideAllc * sampleSize, 64);
    float* omics2DataOrtg = (float*) mkl_malloc(sizeof(float) * omics2ChunkStrideAllc * sampleSize, 64);
    vector<float> covarRowSDCntrl1(covarNum);
    float* covarDataT = (float*) mkl_malloc(sizeof(float) * sampleSize * covarNum, 64);
    float *tau = new float[covarNum];
    float* covarNormSqr = (float*) mkl_malloc(sizeof(float) * sampleSize * sampleSize, 64);
    vector<float> omics1RowSDCntrl2(omics1ChunkStrideAllc, 1);
    vector<float> omics2RowSDCntrl2(omics2ChunkStrideAllc, 1);
    vector<float> omics1Scaling(omics1ChunkStrideAllc, 1);
    vector<float> omics2Scaling(omics2ChunkStrideAllc, 1);
    vector<float> omics1Sum(omics1ChunkStrideAllc, 0);
    vector<float> omics1Sqr(omics1ChunkStrideAllc);
    vector<float> omics1OrtgSqrInv(omics1ChunkStrideAllc);
    vector<float> omics2Sum(omics2ChunkStrideAllc, 0);
    vector<float> covarSum(covarNum, 0);
    vector<vector<float> > omics1DotCov;
    omics1DotCov.resize(omics1ChunkStrideAllc, vector<float>(covarNum));
    vector<vector<float> > omics2DotCov;
    omics2DotCov.resize(omics2ChunkStrideAllc, vector<float>(covarNum));
    vector<vector<float> > CovarInter;
    CovarInter.resize(covarNum, vector<float>(covarNum));
    
    // assign gemm results
    // malloc float type to save space
    float* corrMat = (float*) mkl_malloc(sizeof(float) * (uint32_t) omics1ChunkStrideAllc * omics2ChunkStrideAllc, 64);

    // assign counts of bins of each gradient significant levels conjunction with each distance levels
    uint64_t* st_W = new uint64_t[threadMaxN * (st_ll + 1) * (distLvNum + 1)]();
    vector<vector<uint64_t> > st_Wcdf(distLvNum + 1, vector<uint64_t>(st_ll, 0));

    // assign the number of test of each distance levels
    vector<uint64_t> testCnt(distLvNum + 1, 0);

    // assign the number of significant result
    uint64_t rltN = 0;

#   pragma omp parallel \
    num_threads(threadMaxN) default(shared)
    {
        // get thread id and threads number
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        // assign significant result on each thread
        vector<fitRlt> rltArr; 
        rltArr.reserve((uint64_t)omics1Num * omics2Num * globalP / nthreads);

        // set offset of Inter-variable ST
        uint32_t testCntOffset = tid * (distLvNum + 1);
        uint32_t st_WOffset = tid * (st_ll + 1) * (distLvNum + 1);
        
        // precompute covar data
        if (covarNum > 0) {
            // first centralization of covariates
#           pragma omp for
            for (uint32_t i = 0; i < covarNum; i++) {
                cntrl(covarData, i, sampleSize, covarRowSDCntrl1, sdThd, NASignMarkC);
            }

            // transpose covar
            for (uint32_t i = 0; i < covarNum; i++)
                for (uint32_t j = 0; j < sampleSize; j++)
                    covarDataT[j * covarNum + i] = covarData[i * sampleSize + j];
            // QR decompose
            LAPACKE_sgeqrf(LAPACK_ROW_MAJOR, sampleSize, covarNum, covarDataT, covarNum, tau);
            LAPACKE_sorgqr(LAPACK_ROW_MAJOR, sampleSize, covarNum, covarNum, covarDataT, covarNum, tau);

            // projection
            mkl_set_num_threads_local(threadMaxN); // Specifies the number of threads for MKL
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, sampleSize, sampleSize, covarNum, 
                1, covarDataT, covarNum, covarDataT, covarNum, 0, covarNormSqr, sampleSize);
            mkl_set_num_threads_local(0); // reset the thread-local number to global number
        }
        // calculate inner product between omics data and covariates
        for (uint32_t i = 0; i < covarNum; i++) { CovarInter[i].resize(i + 1); }
        for (uint32_t i = 0; i < covarNum; i++) {
            for (uint32_t j = i * sampleSize; j < (i + 1) * sampleSize; j++) {
                covarSum[i] += covarData[j];
            }
            for (uint32_t j = 0; j <= i; j++) {
                CovarInter[i][j] = cblas_sdot(sampleSize, &covarData[i * sampleSize], 1, &covarData[j * sampleSize], 1);
            }
        }

        for (uint32_t o1Chunk = 0; o1Chunk < omics1ChunkNum; o1Chunk++) {
            for (uint32_t o2Chunk = 0; o2Chunk < omics2ChunkNum; o2Chunk++) {
                // mapping chunk info
                uint32_t omics1ChunkHead, omics2ChunkHead, omics1ChunkStride, omics2ChunkStride;
                uint32_t omics1BlockStrideAllc, omics2BlockStrideAllc;
                omics1ChunkHead = omics1ChunkStrideAllc * o1Chunk;
                omics2ChunkHead = omics2ChunkStrideAllc * o2Chunk;
                omics1ChunkStride = min(omics1ChunkStrideAllc, omics1Num - omics1ChunkHead);
                omics2ChunkStride = min(omics2ChunkStrideAllc, omics2Num - omics2ChunkHead);
                omics1BlockStrideAllc = (uint32_t) (omics1ChunkStride + nthreads - 1) / nthreads;
                omics2BlockStrideAllc = omics2ChunkStride;
    
                // input omics data start
                // input omics1 data when scanned through omics 2
                if (o2Chunk == 0) {
                    if (bfileFlag1) { // input plink bfile as first omics
                        inputOmicsBed(o1InputFile_it, omics1Data, 
                                      omics1ChunkStride, sampleSize, sampleFltSign, covarNANum, 
                                      NASignMark1);
                    } else {
                        input2DfloatParse(o1InputFile, 
                                          omics1Data,
                                          omics1ChunkStride, sampleSize, 
                                          dataArea, 
                                          threadMaxN, 
                                          sampleFltSign, covarNANum, 
                                          NASignMark1, NASign);
                    }
#                   pragma omp single
                    {
                        // move omics2 file point to header
                        if (bfileFlag2) {
                            o2InputFile.seekg(0L, ios::beg);
                            o2InputFile_it = o2InputFile;
                            o2InputFile_it++; o2InputFile_it++; o2InputFile_it++;
                        } else {
                            o2InputFile.seekg(0L, ios::beg);
                        }
                    }
                }

                // input omics2 data
                if (bfileFlag2) { // input plink bfile as second omics
                    inputOmicsBed(o2InputFile_it, omics2Data, 
                                  omics2ChunkStride, sampleSize, sampleFltSign, covarNANum, 
                                  NASignMark2);
                } else {
                    input2DfloatParse(o2InputFile, 
                                      omics2Data,
                                      omics2ChunkStride, sampleSize, 
                                      dataArea, 
                                      threadMaxN, 
                                      sampleFltSign, covarNANum, 
                                      NASignMark2, NASign);
                }
                // input omics data end

                // clocking
                // if (tid == 0) {
                //     time_end = omp_get_wtime();
                //     oss << "Input for chunk " 
                //         << o1Chunk * omics2ChunkNum + o2Chunk + 1 << "/" << chunkNum 
                //         << " has been completed, time used : " << time_end - time_start_whole << " s" << endl;
                //     dualOutput(oss, outputLogFile, std::cout);
                // }

                // orthogonal projection start
                if (o2Chunk == 0) {
                    // first centralization of first omics
#                   pragma omp for 
                    for (uint32_t i = 0; i < omics1ChunkStride; i++) {
                        cntrl(omics1Data, i, sampleSize, omics1RowSDCntrl1, sdThd, NASignMark1);
                    }

                    // generate orthogonal omics data and intermediate variables
#                   pragma omp single
                    {
                        copy(omics1Data, omics1Data + omics1ChunkStride * sampleSize, omics1DataOrtg);
                    }                
                    if (covarNum > 0) {
#                       pragma omp single
                        {
                            // projection
                            mkl_set_num_threads_local(threadMaxN); // Specifies the number of threads for MKL
                            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, omics1ChunkStride, sampleSize, sampleSize, 
                                -1, omics1Data, sampleSize, covarNormSqr, sampleSize, 1, omics1DataOrtg, sampleSize);
                            mkl_set_num_threads_local(0); // reset the thread-local number to global number
                        }
                    }
                    
                    // second centralization for origin omics data
                    // NormMod = 0 : not process/1 : cntrl/2 : cntrlQuant
                    // omics1RowSDCntrl2 is only used as a placeholder here, and will be overwritten later
#                   pragma omp for
                    for (uint32_t i = 0; i < omics1ChunkStride; i++) {
                        if (omics1NormMod == 1) {
                            cntrl(omics1Data, i, sampleSize, omics1RowSDCntrl2, sdThd, NASignMark1);
                        } else if (omics1NormMod == 2) {
                            cntrlQuant(omics1Data, i, sampleSize, omics1RowSDCntrl2, sdThd, NASignMark1);
                        }
                    }

                    // second centralization for orthognal omics data
                    // NormMod = 0 : cntrl/1 : cntrl/2 : cntrlQuant
#                   pragma omp for
                    for (uint32_t i = 0; i < omics1ChunkStride; i++) {
                        if (omics1NormMod != 2) {
                            cntrl(omics1DataOrtg, i, sampleSize, omics1RowSDCntrl2, sdThd, NASignMark1);
                        } else {
                            cntrlQuant(omics1DataOrtg, i, sampleSize, omics1RowSDCntrl2, sdThd, NASignMark1);
                        }
                    }
                    // Scalings of beta coefficents for orthognal omics data
                    // if user do not need to centralize the variates, we should restore the coeffecient to the original scale
#                   pragma omp single 
                    {
                        if (omics1NormMod == 0) {
                            for (uint32_t i = 0; i < omics1ChunkStride; i++) {
                                omics1Scaling[i] = omics1RowSDCntrl1[i] * omics1RowSDCntrl2[i];
                            }
                        }
                    }
                    
                    // preprocess for linear model solver
#                   pragma omp single 
                    {
                        for (uint32_t i = 0; i < omics1ChunkStride; i++) {
                            for (uint32_t j = i * sampleSize; j < (i + 1) * sampleSize; j++) {
                                omics1Sum[i] += omics1Data[j];
                            }
                            omics1Sqr[i] = cblas_sdot(sampleSize, &omics1Data[i * sampleSize], 1, &omics1Data[i * sampleSize], 1);
                            omics1OrtgSqrInv[i] = 1 / cblas_sdot(sampleSize, &omics1DataOrtg[i * sampleSize], 1, &omics1DataOrtg[i * sampleSize], 1);
                        }
                    }

                    // calculate inner product between omics data and covariates
#                   pragma omp single 
                    {
                        for (uint32_t i = 0; i < covarNum; i++) {
                            for (uint32_t j = 0; j < omics1ChunkStride; j++) {
                                omics1DotCov[j][i] = cblas_sdot(sampleSize, &omics1Data[j * sampleSize], 1, &covarData[i * sampleSize], 1);
                            }
                        }
                    }
                }

                // first centralization of second omics
#               pragma omp for
                for (uint32_t i = 0; i < omics2ChunkStride; i++) {
                    cntrl(omics2Data, i, sampleSize, omics2RowSDCntrl1, sdThd, NASignMark2);
                }

                // generate orthogonal omics data and intermediate variables
#               pragma omp single
                {
                    copy(omics2Data, omics2Data + omics2ChunkStride * sampleSize, omics2DataOrtg);
                }          
                if (covarNum > 0) {
#                   pragma omp single
                    {
                        // projection
                        mkl_set_num_threads_local(threadMaxN); // Specifies the number of threads for MKL
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, omics2ChunkStride, sampleSize, sampleSize, 
                            -1, omics2Data, sampleSize, covarNormSqr, sampleSize, 1, omics2DataOrtg, sampleSize);
                        mkl_set_num_threads_local(0); // reset the thread-local number to global number
                    }
                }
                
                // second centralization for origin omics data
                // NormMod = 0 : not process/1 : cntrl/2 : cntrlQuant
                // omics2RowSDCntrl2 is only used as a placeholder here, and will be overwritten later
#               pragma omp for
                for (uint32_t i = 0; i < omics2ChunkStride; i++) {
                    if (omics2NormMod == 1) {
                        cntrl(omics2Data, i, sampleSize, omics2RowSDCntrl2, sdThd, NASignMark2);
                    } else if (omics2NormMod == 2) {
                        cntrlQuant(omics2Data, i, sampleSize, omics2RowSDCntrl2, sdThd, NASignMark2);
                    }
                }

                // second centralization for orthognal omics data
                // NormMod = 0 : cntrl/1 : cntrl/2 : cntrlQuant
#               pragma omp for
                for (uint32_t i = 0; i < omics2ChunkStride; i++) {
                    if (omics2NormMod != 2) {
                        cntrl(omics2DataOrtg, i, sampleSize, omics2RowSDCntrl2, sdThd, NASignMark2);
                    } else {
                        cntrlQuant(omics2DataOrtg, i, sampleSize, omics2RowSDCntrl2, sdThd, NASignMark2);
                    }
                }
                // Scalings of beta coefficents for orthognal omics data
                // if user do not need to centralize the variates, we should restore the coeffecient to the original scale
#               pragma omp single 
                {
                    if (omics2NormMod == 0) {
                        for (uint32_t i = 0; i < omics2ChunkStride; i++) {
                            omics2Scaling[i] = omics2RowSDCntrl1[i] * omics2RowSDCntrl2[i];
                        }
                    }
                }

                // preprocess for linear model solver
#               pragma omp single 
                {
                    for (uint32_t i = 0; i < omics2ChunkStride; i++) {
                        for (uint32_t j = i * sampleSize; j < (i + 1) * sampleSize; j++) {
                            omics2Sum[i] += omics2Data[j];
                        }
                    }
                }

                // calculate inner product between omics data and covariates
#               pragma omp single 
                {
                    for (uint32_t i = 0; i < covarNum; i++) {
                        for (uint32_t j = 0; j < omics2ChunkStride; j++) {
                            omics2DotCov[j][i] = cblas_sdot(sampleSize, &omics2Data[j * sampleSize], 1, &covarData[i * sampleSize], 1);
                        }
                    }
                }
                // orthogonal projection end

                // clocking
                // if (tid == 0) {
                //     time_end = omp_get_wtime();
                //     oss << "Orthogonal projection for chunk " 
                //         << o1Chunk * omics2ChunkNum + o2Chunk + 1 << "/" << chunkNum 
                //         << " has been completed, time used : " << time_end - time_start_whole << " s" << endl;
                //     dualOutput(oss, outputLogFile, std::cout);    
                // }
                
                // mapping block info of calling thread
                uint32_t omics1BlockStride, omics2BlockStride, omics1BlockHead, omics2BlockHead;
                uint64_t pairAmt;
                omics1BlockStride = min(omics1BlockStrideAllc, omics1ChunkStride - omics1BlockStrideAllc * tid);
                omics2BlockStride = omics2BlockStrideAllc;
                omics1BlockHead = tid * omics1BlockStrideAllc;
                omics2BlockHead = 0;
                pairAmt = omics1BlockStride * omics2BlockStride;

                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                            omics1BlockStride, omics2BlockStride, sampleSize, 
                            1, &omics1DataOrtg[omics1BlockHead * sampleSize], sampleSize, omics2DataOrtg, sampleSize, 
                            0, &corrMat[tid * omics1BlockStride * omics2BlockStride], omics2BlockStride);

                float *corrCurr = &corrMat[tid * omics1BlockStride * omics2BlockStride];
                for (uint32_t i = 0; i < omics1BlockStride; i++) {
                    for (uint32_t j = 0; j < omics2BlockStride; j++) {
                        uint32_t omics1VarInChunk = i + omics1BlockHead;
                        uint32_t omics2VarInChunk = j + omics2BlockHead;
                        uint32_t omics1VarGlobal = omics1VarInChunk + omics1ChunkHead;
                        uint32_t omics2VarGlobal = omics2VarInChunk + omics2ChunkHead;
                        float corrCurrAbs = abs(*corrCurr);

                        // find corresponding RCV bin for corrCurr
                        // need to recheck
                        uint32_t st_WOffset_sec;
                        if (corrCurrAbs < 21.3322) { // st_lambdaRCV[10]
                            if (corrCurrAbs < 10.0785) { // st_lambdaRCV[15]
                                if (corrCurrAbs < 3.97474) { // st_lambdaRCV[18]
                                    if (corrCurrAbs < 1.98346) {
                                        st_WOffset_sec = st_WOffset + 19 * (distLvNum + 1);
                                    } else {
                                        st_WOffset_sec = st_WOffset + 18 * (distLvNum + 1);
                                    }
                                } else {
                                    if (corrCurrAbs < 5.9819) { // st_lambdaRCV[17]
                                        st_WOffset_sec = st_WOffset + 17 * (distLvNum + 1);
                                    } else {
                                        if (corrCurrAbs < 8.01342) { // st_lambdaRCV[16]
                                            st_WOffset_sec = st_WOffset + 16 * (distLvNum + 1);
                                        } else {
                                            st_WOffset_sec = st_WOffset + 15 * (distLvNum + 1);
                                        }
                                    }
                                }
                            } else {
                                if (corrCurrAbs < 14.3521) { // st_lambdaRCV[13]
                                    if (corrCurrAbs < 12.1875) { // st_lambdaRCV[14]
                                        st_WOffset_sec = st_WOffset + 14 * (distLvNum + 1);
                                    } else {
                                        st_WOffset_sec = st_WOffset + 13 * (distLvNum + 1);
                                    }
                                } else {
                                    if (corrCurrAbs < 16.586) { // st_lambdaRCV[12]
                                        st_WOffset_sec = st_WOffset + 12 * (distLvNum + 1);
                                    } else {
                                        if (corrCurrAbs < 18.9059) { // st_lambdaRCV[11]
                                            st_WOffset_sec = st_WOffset + 11 * (distLvNum + 1);
                                        } else {
                                            st_WOffset_sec = st_WOffset + 10 * (distLvNum + 1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if (corrCurrAbs < 36.3743) { // st_lambdaRCV[5]
                                if (corrCurrAbs < 26.6164) { // st_lambdaRCV[8]
                                    if (corrCurrAbs < 23.8909) { // st_lambdaRCV[9]
                                        st_WOffset_sec = st_WOffset + 9 * (distLvNum + 1);
                                    } else {
                                        st_WOffset_sec = st_WOffset + 8 * (distLvNum + 1);
                                    }
                                } else {
                                    if (corrCurrAbs < 29.5553) { // st_lambdaRCV[7]
                                        st_WOffset_sec = st_WOffset + 7 * (distLvNum + 1);
                                    } else {
                                        if (corrCurrAbs < 32.7743) { // st_lambdaRCV[6]
                                            st_WOffset_sec = st_WOffset + 6 * (distLvNum + 1);
                                        } else {
                                            st_WOffset_sec = st_WOffset + 5 * (distLvNum + 1);
                                        }
                                    }
                                }
                            } else {
                                if (corrCurrAbs < 51.9926) { // st_lambdaRCV[2]
                                    if (corrCurrAbs < 40.5197) { // st_lambdaRCV[4]
                                        st_WOffset_sec = st_WOffset + 4 * (distLvNum + 1);
                                    } else {
                                        if (corrCurrAbs < 45.5098) { // st_lambdaRCV[3]
                                            st_WOffset_sec = st_WOffset + 3 * (distLvNum + 1);
                                        } else {
                                            st_WOffset_sec = st_WOffset + 2 * (distLvNum + 1);
                                        }
                                    }
                                } else {
                                    if (corrCurrAbs < 61.9354) { // st_lambdaRCV[1]
                                        st_WOffset_sec = st_WOffset + (distLvNum + 1);
                                    } else {
                                        st_WOffset_sec = st_WOffset + 0;
                                    }
                                }
                            }
                        }

                        int8_t levelCurr = -1;
                        if (distLvNum >= 1 &&
                            omics1CHR[omics1VarGlobal] == omics2CHR[omics2VarGlobal] && 
                            omics1CHR[omics1VarGlobal] > 0 && omics2CHR[omics2VarGlobal] > 0 && 
                            abs(omics1BP[omics1VarGlobal] - omics2BP[omics2VarGlobal]) <= distLv[distLvNum - 1] && 
                            omics1BP[omics1VarGlobal] > 0 && omics2BP[omics2VarGlobal] > 0) { // two locus are localed in the broadest level distance
                            for (uint8_t l = 0; l < distLvNum; l++) {
                                if (abs(omics1BP[omics1VarGlobal] - omics2BP[omics2VarGlobal]) <= distLv[l]) {
                                    st_W[st_WOffset_sec + l] += 1;
                                    if (corrCurrAbs >= rCriticalValue[l]) { levelCurr = l; }
                                    break;
                                }
                            }
                        } else { // two locus are not localed in any level distance
                            st_W[st_WOffset_sec + distLvNum] += 1;
                            if (corrCurrAbs >= rCriticalValue[distLvNum]) {
                                levelCurr = distLvNum;
                            }
                        }

                        if (levelCurr >= 0) {
                            fitRlt rltTmp = 
                            linearFit(*corrCurr, 
                                      omics1VarInChunk, omics2VarInChunk, 
                                      omics1VarGlobal, omics2VarGlobal, 
                                      levelCurr, 
                                      sampleSize, covarNum, 
                                      msRtThd,
                                      omics1Data, omics2Data, covarData,
                                      NASignMark1, NASignMark2, 
                                      omics1Sum, omics2Sum, covarSum, 
                                      omics1Sqr, omics1OrtgSqrInv, 
                                      omics1DotCov, omics2DotCov, CovarInter,
                                      omics1NormMod, omics2NormMod, 
                                      omics1Scaling, omics2Scaling, 
                                      omics1RowSDCntrl1, omics2RowSDCntrl1);
                            if (rltTmp.p <= distLvP[levelCurr] && rltTmp.status == 0) {
                                rltArr.push_back(rltTmp);
                            }
                        }

                        // next correlation
                        corrCurr++;
                    }
                }

#               pragma omp critical
                {
                    // output results
                    outputBinFile.write((char*)rltArr.data(), rltArr.size() * sizeof(fitRlt));

                    // wipe data, keep capacity
                    rltN += rltArr.size(); 
                    rltArr.clear();
                }

                // clocking
                if (tid == 0) {
                    time_end = omp_get_wtime();
                    oss << "Main process for chunk " 
                        << o1Chunk * omics2ChunkNum + o2Chunk + 1 << "/" << chunkNum 
                        << " has been completed, time used : " << time_end - time_start_whole << " s" << endl;
                    dualOutput(oss, outputLogFile, std::cout);
                }
            }
        }

        // free rltArr
        rltArr.clear(); // clear elements in vectors
        vector<fitRlt>().swap(rltArr); // free vectors

        // synchronize threads
#       pragma omp barrier

        // reduce bin count array
#       pragma omp single
        {
            for (uint8_t l = 0; l < (distLvNum + 1); l++) {
                for (uint32_t i = 0; i < nthreads; i++) {
                    st_Wcdf[l][st_ll - 1] += st_W[i * (st_ll + 1) * (distLvNum + 1) + st_ll * (distLvNum + 1) + l];
                }
                for (int32_t j = st_ll - 2; j >= 0; j--) { // 
                    st_Wcdf[l][j] = st_Wcdf[l][j + 1];
                    for (uint32_t i = 0; i < nthreads; i++) {
                        st_Wcdf[l][j] += st_W[i * (st_ll + 1) * (distLvNum + 1) + (j + 1) * (distLvNum + 1) + l];
                    }
                }
            }

            for(uint8_t l = 0; l < (distLvNum + 1); l++) {
                for(uint32_t i = 0; i < nthreads; i++) {
                    for (int32_t j = 0; j < st_ll + 1; j++) {
                        testCnt[l] += st_W[i * (st_ll + 1) * (distLvNum + 1) + j * (distLvNum + 1) + l];
                    }
                }
            }
        }
    }

    // close outputBinFile
    outputBinFile.close();

    // memory free
    mkl_free(omics1Data); mkl_free(omics2Data);
    delete[] st_W;

    // ST FDR
    // estimate pi0
    vector<double> st_pDens(st_ll), st_pDensSort(st_ll);
    vector<double > st_mse(st_ll);
    vector<double> st_pi0(1 + distLvNum, 1);

    for (uint8_t l = 0; l < 1 + distLvNum; l++) {
        for (uint32_t i = 0; i < st_ll; i++) {
            st_pDens[i] = st_Wcdf[l][i] / (testCnt[l] * (1 - st_lambda[i]));
        }

        st_pDensSort = st_pDens;
        sort(st_pDensSort.begin(), st_pDensSort.end());
        double minpi0 = gsl_stats_quantile_from_sorted_data(st_pDensSort.data(), 1, st_ll, 0.1);
        
        for (uint32_t i = 0; i < st_ll; i++) {
            if (st_Wcdf[l][i] == 0) {
                st_mse[i] = 999L; // extract elements which Wcdf == 0
            } else {
                st_mse[i] = (st_Wcdf[l][i] / (pow(testCnt[l], 2) * pow(1 - st_lambda[i], 2))) * 
                            (1 - 1.0*st_Wcdf[l][i]/testCnt[l]) + 
                            pow(st_pDens[i] - minpi0, 2);
            }
        }

        double st_mseMin = *min_element(st_mse.begin(), st_mse.end());
        for (uint32_t i = 0; i < st_ll; i++) {
            if (fabs(st_mse[i] - st_mseMin) < 1e-6 & st_pi0[l] > st_pDens[i]) {
                st_pi0[l] = st_pDens[i];
            }
        }
    }
    
    oss << endl;
    for (uint8_t l = 0; l < 1 + distLvNum; l++) {
        if (st_pi0[l] < 0) {
            oss << "Warning: The estimated PI0 <= 0 in distance level " << l << " . To keep the program running, pi0 is set to 1. \n" << endl; dualOutput(oss, outputLogFile, std::cout);
            st_pi0[l] = 1;
        }
        oss << "PI0 of distance level " << l + 1 << " : " << st_pi0[l] << endl;
    }
    oss << endl; dualOutput(oss, outputLogFile, std::cout);

    // assign P vector for all significant results
    vector<fitRltPart> rltArrP(rltN);
    // input level and P
    ifstream inputBinFile(outputFileName + ".bin", ios::in | ios::binary);
    inputBinFile.seekg(0); // set get point of inputBinFile into header
    for (uint64_t i = 0; i < rltN; i++) {
        fitRlt rltTmp;
        if (!inputBinFile.read((char*)&rltTmp, sizeof(fitRlt))) {
            oss << "Error: Output file is incomplete .\n"; dualOutput(oss, outputLogFile, std::cout);
            return 1;
        }
        rltArrP[i].level = rltTmp.level;
        rltArrP[i].p = rltTmp.p;
        rltArrP[i].rank = i;
    }

    // calc Q-value
    if (rltN > 0) {
        sort(rltArrP.begin(), rltArrP.end());
#       pragma omp parallel \
        num_threads(threadMaxN) default(shared) 
        {
            int64_t rltHead = 0, rltTail;
            for (uint8_t l = 0; l < 1 + distLvNum; l++) {
                for (rltTail = rltHead; rltArrP[rltTail].level == l + 1; rltTail++);
#               pragma omp for schedule(dynamic) 
                for (int64_t i = rltHead; i < rltTail; i++) {
                    // use p to storage q value
                    rltArrP[i].p = rltArrP[i].p * testCnt[l]/((i - rltHead + 1) * (1 - pow(1 - rltArrP[i].p, testCnt[l])));
                }
#               pragma omp for schedule(dynamic) 
                for (int64_t i = rltTail - 2; i >= rltHead; i--) {
                    if (rltArrP[i].p > rltArrP[i+1].p) {
                        rltArrP[i].p = rltArrP[i+1].p;
                    }            
                }
#               pragma omp for schedule(dynamic) 
                for (int64_t i = rltHead; i < rltTail; i++) {
                    rltArrP[i].p = min(1.0, rltArrP[i].p) * st_pi0[l];
                }
                rltHead = rltTail;
            }
        }
    }

    // translate rank to order of rltArrP
    // decoupling .p and .rank element of rltArrP
    vector<uint8_t> rltArrSign(rltN, 0);
    for (uint64_t i = 0; i < rltN; i++) {
        if (rltArrSign[i] == 0) {
            uint64_t k = i;
            uint64_t j = rltArrP[i].rank;
            while (rltArrSign[j] == 0) {
                rltArrSign[j] = 1;
                uint64_t r_1 = rltArrP[j].rank;
                rltArrP[j].rank = k;
                k = j;
                j = r_1;
            }
        }
    }

    // append Q-value to result file
    inputBinFile.seekg(0); // set get point of inputBinFile into header
    ofstream outputFile(outputFileName);
    outputFile << setprecision(outPcs);
    outputFile << "omics1\t" << "omics2\t" << "distance_level\t" << "NMISS\t" << "BETA\t" << "SE\t" << "T\t" << "P-value\t" << "Q-value\n";
    for (uint64_t i = 0; i < rltN; i++) {
        fitRlt rltTmp;
        inputBinFile.read((char*)&rltTmp, sizeof(fitRlt));

        outputFile << omics1Name[rltTmp.omics1Id] << "\t";
        outputFile << omics2Name[rltTmp.omics2Id] << "\t";
        outputFile << (uint32_t)rltTmp.level << "\t";
        outputFile << rltTmp.nmiss << "\t";
        outputFile << rltTmp.b << "\t";
        outputFile << rltTmp.se << "\t";
        outputFile << rltTmp.t << "\t";
        outputFile << rltTmp.p << "\t";
        outputFile << rltArrP[rltArrP[i].rank].p << "\n";
    }
    inputBinFile.close(); remove((outputFileName + ".bin").c_str());
    outputFile.close();
    
    // free rltArrP
    rltArrP.clear();
    vector<fitRltPart>().swap(rltArrP);

    // output count of test
    for (uint8_t l = 0; l < distLvNum + 1; l++) { 
        oss << "Number of test of distance level " << l + 1 << " : " << testCnt[l] << endl; dualOutput(oss, outputLogFile, std::cout);
    }

    // global ending time stamp
    time_end = omp_get_wtime();
    oss << endl << "Whole procedure time used : " << time_end - time_start_whole << " s" << endl << endl; dualOutput(oss, outputLogFile, std::cout);

    return 0;
}
