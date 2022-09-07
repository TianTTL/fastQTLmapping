#ifndef _SNPLIB_SRC_SNP_H_
#define _SNPLIB_SRC_SNP_H_

#include <array>
#include <cstring>
#include <vector>

using namespace std;

namespace meqtllib {
// structure of linear regression result
struct linearFitRlt {
    uint32_t omics1Id;
    uint32_t omics2Id;
    float b;
    float t;
    double p;
    double q;
    float se; // beta / t
    uint32_t nmiss;
    uint32_t status; // 0: good; 1: missing rate error
    uint32_t level;
};

// global variables
string VERSION="0.9.3";
string omics1FileName, omics2FileName, outputFileName, rplFileName = "NA";
string covarFileName;
bool bfileFlag1 = false, bfileFlag2 = false;
string NASign = "NA";
vector<uint32_t> categFlag;

double globalP = 1;
double FWER = 0.05;
float missingRateThd = 0.1;

vector<double> distLv;
vector<double> distLvP;
uint32_t distLvNum;

uint32_t threadMaxN = 1;
int32_t omics1NormMod = 0, omics2NormMod = 0;
double PLooseMarg = 100;
uint32_t outPcs = 4;
int32_t helpFlag = 0;
int32_t modeFlag = 0; // 1 cnt, 2 disc, 3 rpl
uint32_t chunkSize = 5000;

uint32_t omics1Num, omics2Num, covarNum, sampleSize;
ofstream outputLogFile;
ostringstream oss;

void calcBfileSize(string bfileNameRoot, uint32_t &sampleSize, uint32_t &omicsSize);
void getBfileSNPid(string bfileNameRoot, uint32_t num_snps, 
                   vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP);
void inputOmicsBed(istreambuf_iterator<uint8_t>& inputFile,
                       double *omicsData, 
                       int locusCount, int sampleSize, vector<bool>& sampleFltSign, uint32_t sampleFltNum, 
                       vector<vector<uint32_t> >& NA_data);
void calcInputSize(string omicsFileName, uint32_t &sampleSize, uint32_t& omicsNum);
void get2DfloatId(string fileName, uint32_t omicsNum, 
                  vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP);
void input2DfloatParse(std::ifstream& inputFile, 
                       double* omicsData, 
                       uint32_t locusCount, uint32_t sampleSize, 
                       string* dataArea, 
                       uint32_t threadMaxN, 
                       vector<bool>& sampleFltSign, uint32_t sampleFltNum, 
                       vector<vector<uint32_t> >& NASignMark, string NASign);
void calcCovarSize(string covarFileName, string NASign, uint32_t sampleSize, uint32_t& covarNum, 
                   vector<bool>& sampleFltSign, uint32_t& covarNANum, 
                   vector<uint32_t>& categFlag, uint32_t& covarCategNum);
double* inputCovar(string fileName, 
                  uint32_t& covarNum, uint32_t sampleSize, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum, 
                  vector<uint32_t>& categFlag, uint32_t covarCategNum);
template <class ForwardIterator, class T>
int64_t binarySearch(ForwardIterator head, ForwardIterator tail, const T& val);
template <class ForwardIterator, class T>
int64_t binarySearchDec(ForwardIterator head, ForwardIterator tail, const T& val);
void inputRplList(string rplFileName, vector<pair<int64_t, int64_t> >& rplList, 
                  vector<string>& omics1Name, vector<string>& omics2Name, uint32_t threadMaxN);
void cntrl(double* a, uint32_t omicsId, 
           uint32_t sampleSize, 
           vector<double>& rowSD,
           vector<vector<uint32_t> >& NASignMarkCurr);
void cntrlQuant(double *omicsData, uint32_t omicsId, uint32_t sampleSize,
                vector<vector<uint32_t> >& NASignMarkCurr);
vector<double> rankSort(const vector<double>& v_temp, uint32_t sampleSize);
bool cmpLinearFitRlt(struct linearFitRlt &A, struct linearFitRlt &B);
void rCriticalValueCalc(double P, uint32_t sampleSize, double &rCriticalValue);
linearFitRlt linearFit(double corr, 
                       uint32_t omics1Id, uint32_t omics2Id, 
                       uint32_t omics1GlobalId, uint32_t omics2GlobalId, uint32_t level, 
                       uint32_t sampleSize, uint32_t covarNum, 
                       float missingRateThd,
                       double* omics1Data, double* omics2Data, double* covarData, 
                       vector<vector<uint32_t> >& NASignMark1, vector<vector<uint32_t> >& NASignMark2, 
                       vector<double>& omics1Sum, vector<double>& omics2Sum, vector<double>& covarSum, 
                       vector<double>& omics1Sqr, vector<double>& omics1OrtgSqrInv, 
                       vector<vector<double> >& omics1DotCov, vector<vector<double> >& omics2DotCov, vector<vector<double> >& CovarInter,
                       int32_t omics1NormMod, int32_t omics2NormMod, 
                       vector<double>& omics1Scaling, vector<double>& omics2Scaling, 
                       vector<double>& omics1RowSDCntrl1, vector<double>& omics2RowSDCntrl1);
void rCriticalValueCalc(double P, uint32_t sampleSize, double &rCriticalValue);
}  // namespace meqtllib

#endif  //_SNPLIB_SRC_SNP_H_
