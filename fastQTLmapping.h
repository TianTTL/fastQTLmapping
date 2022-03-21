#ifndef _SNPLIB_SRC_SNP_H_
#define _SNPLIB_SRC_SNP_H_

#include <array>
#include <cstring>
#include <vector>

using namespace std;

namespace snplib {
class SNP {
 private:
  const uint8_t *geno_;
  const size_t num_samples_;
  const size_t num_full_bytes_;
  const size_t num_samples_left_;
  const size_t num_bytes_;

  void ConvertGeno(size_t num_snps, size_t idx, std::array<uint64_t, 4> &g64);
  void ConvertGeno(size_t num_snps, size_t idx, std::array<uint64_t, 5> &g64);

 public:
  SNP(const uint8_t *geno, size_t num_samples);
  SNP(const SNP &) = default;
  SNP(SNP &&) = default;
  SNP &operator=(const SNP &) = delete;
  SNP &operator=(SNP &&) = delete;
  ~SNP() = default;
  uint8_t operator[](size_t idx) const { return geno_[idx]; }
  uint8_t operator()(size_t idx) const {
    auto i = idx / 4;
    auto s = idx % 4;
    return ((*this)[i] >> (2 * s)) & 3u;
  }
  SNP &operator+=(size_t idx) {
    geno_ += idx * num_bytes_;
    return *this;
  }
  template <class T>
  void Copy(T *dest) const {
    memcpy((void *)dest, (const void *)geno_, sizeof(uint8_t) * num_bytes_);
  }
  template <class T>
  void UnpackGeno(const std::array<T, 4> &geno_table, T *geno) {
    for (size_t i = 0; i < num_full_bytes_; ++i) {
      auto t = geno_[i];
      geno[4 * i] = geno_table[t & 3u];
      t >>= 2;
      geno[4 * i + 1] = geno_table[t & 3u];
      t >>= 2;
      geno[4 * i + 2] = geno_table[t & 3u];
      t >>= 2;
      geno[4 * i + 3] = geno_table[t & 3u];
    }
    if (num_samples_left_ > 0u) {
      auto t = geno_[num_full_bytes_];
      for (size_t i = 0; i < num_samples_left_; ++i) {
        geno[4 * num_full_bytes_ + i] = geno_table[t & 3u];
        t >>= 2;
      }
    }
  }
  template <class T>
  void UnpackGeno(const std::array<T, 4> &geno_table,
                  const std::array<T, 4> &mask_table, T *geno, T *mask) {
    for (size_t i = 0; i < num_full_bytes_; ++i) {
      auto t = geno_[i];
      geno[4 * i] = geno_table[t & 3u];
      mask[4 * i] = mask_table[t & 3u];
      t >>= 2;
      geno[4 * i + 1] = geno_table[t & 3u];
      mask[4 * i + 1] = mask_table[t & 3u];
      t >>= 2;
      geno[4 * i + 2] = geno_table[t & 3u];
      mask[4 * i + 2] = mask_table[t & 3u];
      t >>= 2;
      geno[4 * i + 3] = geno_table[t & 3u];
      mask[4 * i + 3] = mask_table[t & 3u];
    }
    if (num_samples_left_ > 0u) {
      auto t = geno_[num_full_bytes_];
      for (size_t i = 0; i < num_samples_left_; ++i) {
        geno[4 * num_full_bytes_ + i] = geno_table[t & 3u];
        mask[4 * num_full_bytes_ + i] = mask_table[t & 3u];
        t >>= 2;
      }
    }
  }
  void TransposeGeno(size_t num_snps, size_t idx, uint64_t *geno64);
  void TransposeGeno(size_t num_snps, uint64_t *geno64);
};
}  // namespace snplib

namespace meqtllib {
/* Parameters */
#define VARNUM 2

struct linearFitRlt {
    uint64_t omics1Id;
    uint64_t omics2Id;
    float b;
    float t;
    double p;
    float se; // beta / t
    uint32_t nmiss;
    uint32_t status; // 0: good; 1: missing rate error
    uint32_t level;
};

uint64_t omics1Num;
uint64_t omics2Num;
uint32_t covarNum;
uint32_t sampleSize;
float missingRateThd, MAFThd;
const uint64_t blockSize = 10000;
const uint32_t precision_config = 2;

void calcBfileSize(string bfileNameRoot, uint32_t &num_samples, uint64_t &num_snps);
void getBfileSNPid(string bfileNameRoot, uint64_t num_snps, 
                   vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP);
void calcInputSize(string omicsFileName, uint64_t& omicsNum);
void input2DfloatParse(double* omicsData, string fileName, vector<vector<uint32_t> >& NASignMark, string NASign, 
                  vector<string>& omicsName, vector<int32_t>& omicsCHR, vector<int64_t>& omicsBP, 
                  uint64_t omicsNum, uint32_t sampleSize, 
                  string* dataArea, 
                  uint32_t threadMaxN, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum);
void calcCovarSize(string covarFileName, string NASign, uint64_t sampleSize, uint32_t& covarNum, 
                   vector<bool>& sampleFltSign, uint32_t& covarNANum, 
                   vector<int32_t>& categFlag, uint32_t& covarCategNum);
double* inputCovar(string fileName, 
                  uint32_t& covarNum, uint32_t sampleSize, 
                  vector<bool>& sampleFltSign, uint32_t covarNANum, 
                  vector<int32_t>& categFlag, uint32_t covarCategNum);
void inputRplList(string rplFileName, vector<pair<int64_t, int64_t> >& rplList, 
                  vector<string>& omics1Name, vector<string>& omics2Name, uint32_t threadMaxN);
void cntrl(double* a, uint64_t omicsId, 
           uint32_t sampleSize, 
           vector<double>& rowSD,
           vector<vector<uint32_t> >& NASignMarkCurr);
void cntrlQuant(double *omicsData, uint64_t omicsId, uint32_t sampleSize,
                vector<vector<uint32_t> >& NASignMarkCurr);
vector<double> rankSort(const vector<double>& v_temp, uint64_t sampleSize);
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
                       vector<double>& omics1RowSDCntrl1, vector<double>& omics2RowSDCntrl1);
void rCriticalValueCalc(double P, uint32_t sampleSize, double &rCriticalValue);
}  // namespace meqtllib

#endif  //_SNPLIB_SRC_SNP_H_
