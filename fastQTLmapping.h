#ifndef _SNPLIB_SRC_SNP_H_
#define _SNPLIB_SRC_SNP_H_

#include <array>
#include <cstring>
#include <vector>

using namespace std;

namespace meqtllib {
// structure of linear regression result
struct fitRlt {
    uint32_t omics1Id;
    uint32_t omics2Id;
    uint8_t level;
    uint32_t nmiss;
    float b;
    float t;
    double p;
    float se; // beta / t
    uint32_t status; // 0: good; 1: missing rate error
};

// structure of partial of linear regression result
struct fitRltPart {
    uint8_t level;
    double pq;
    int64_t rank;
};

// global variables
string VERSION="0.9.9_alpha";
string omics1FileName, omics2FileName, outputFileName, rplFileName = "NA";
string covarFileName;
bool bfileFlag1 = false, bfileFlag2 = false;
string NASign = "NA";
vector<uint32_t> categFlag;

double globalP = 1.0;
float FWER = 0.05;
float msRtThd = 0.1;
double sdThd = 1e-6;

vector<int64_t> distLv;
vector<double> distLvP;
uint8_t distLvNum;

uint32_t threadMaxN = 1;
int32_t omics1NormMod = 0, omics2NormMod = 0;
float PLooseMarg = 100;
uint32_t outPcs = 4;
int32_t helpFlag = 0;
int32_t modeFlag = 0; // 1 pre-analysis, 2 disc, 3 rpl
uint32_t chunkSize = 5000;

uint32_t omics1Num, omics2Num, covarNum, sampleSize, sampleSizeAlt;
ofstream outputLogFile;
ostringstream oss;

void calcBfileSize(string bfileNameRoot, uint32_t &sampleSize, uint32_t &omicsSize);
void getBfileSNPid(string bfileNameRoot, uint32_t num_snps, 
                   vector<string>& omicsName, vector<int32_t>& omicsCHR, 
                   vector<int64_t>& omicsBPST, vector<int64_t>& omicsBPEN);
void inputOmicsBed(istreambuf_iterator<char>& inputFile,
                   float *omicsData, 
                   uint32_t locusCount, uint32_t sampleSize, vector<bool>& sampleFltSign, uint32_t sampleFltNum, 
                   vector<vector<uint32_t> >& NASignMark);
void calcInputSize(string omicsFileName, uint32_t &sampleSize, uint32_t& omicsNum);
uint32_t get2DfloatId(string fileName, uint32_t omicsNum, 
                      vector<string>& omicsName, vector<int32_t>& omicsCHR, 
                      vector<int64_t>& omicsBPST, vector<int64_t>& omicsBPEN);
void input2DfloatParse(std::ifstream& inputFile, 
                       float* omicsData, 
                       uint32_t locusCount, uint32_t sampleSize, 
                       string* dataArea, 
                       uint32_t threadMaxN, 
                       vector<bool>& sampleFltSign, uint32_t sampleFltNum, 
                       vector<vector<uint32_t> >& NASignMark, string NASign);
void calcCovarSize(string covarFileName, string NASign, uint32_t &sampleSize, uint32_t& covarNum, 
                   vector<bool>& sampleFltSign, uint32_t& covarNANum, 
                   vector<uint32_t>& categFlag, uint32_t& covarCategNum);
float* inputCovar(string fileName, 
                  uint32_t& covarNum, uint32_t sampleSize, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum, 
                  vector<uint32_t>& categFlag, uint32_t covarCategNum);
template <class ForwardIterator, class T>
int64_t binarySearch(ForwardIterator head, ForwardIterator tail, const T& val);
template <class ForwardIterator, class T>
int64_t binarySearchDec(ForwardIterator head, ForwardIterator tail, const T& val);
void inputRplList(string rplFileName, vector<pair<int64_t, int64_t> >& rplList, 
                  vector<string>& omics1Name, vector<string>& omics2Name, uint32_t threadMaxN);
template <typename T>
void norm(T* omicsData, uint32_t omicsId, 
           uint32_t sampleSize, 
           vector<double>& rowSD,
           double sdThd, 
           vector<vector<uint32_t> >& NASignMarkCurr);
template <typename T>
void normQuant(T *omicsData, uint32_t omicsId, 
                uint32_t sampleSize,
                vector<double>& rowSD, 
                double sdThd, 
                vector<vector<uint32_t> >& NASignMarkCurr);
template <typename T>
vector<float> rankSort(const vector<T>& v_temp, uint32_t sampleSize);
void rCriticalValueCalc(double P, uint32_t sampleSize, uint32_t covarNum, double &rCriticalValue);
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
                 vector<float>& omics1RowSDNorm1, vector<float>& omics2RowSDNorm1, 
                 vector<float>& omics1RowSDNorm2, vector<float>& omics2RowSDNorm2);
void rCriticalValueCalc(double P, uint32_t sampleSize, double &rCriticalValue);
}  // namespace meqtllib

#endif  //_SNPLIB_SRC_SNP_H_
