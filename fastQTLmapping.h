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

struct linearFitRlt{
    int currentOmics1;
    int currentOmics2;
    float b;
    float t;
    double p;
    float se; // beta / t
    float r2; 
    int nmiss;
    int status; // 0: good; 1: multicollinearity; 2: constant value; 3: matrix inversion error
};

long omics1Num;
long omics2Num;
int sampleSize;
float missingRateThd, MAFThd;
const long blockSize = 10000;

void calcBfileSize(char *bfileNameRoot, int &num_samples, long &num_snps);
void getBfileSNPid(char *bfileNameRoot, int num_snps, vector<string>& omicsName);
void calcInputSize(char* omicsFileName, long& omicsNum, int& sampleSize);
void input2Dfloat(float* omicsData, char* fileName, vector<vector<int> >& NASignMark, char* NASign, 
                   vector<string>& omicsName, int omicsNum, int sampleSize, 
                   string* dataArea, 
                   int thread_count);
void preprocessing(float* omicsData, float* omicsDataNorm, 
                   vector<double>& omicsRowSum, vector<double>& omicsRowSD, 
                   int omicsNum, int sampleSize, vector<vector<int> >& NASignMarkCurr);
linearFitRlt linearFit(int currentOmics1, int currentOmics2, 
                     int sampleSize, float MAFThd, float missingRateThd,
                     float* omics1Data, float* omics2Data, 
                     float* omics1DataCurr, float* omics2DataCurr, 
                     vector<vector<int> >& NASignMark1, vector<vector<int> >& NASignMark2,
                     vector<double>& omics1RowSum, vector<double>& omics2RowSum);
void rCriticalValueCalc(double P, int sampleSize, double &rCriticalvalue);
}  // namespace meqtllib

#endif  //_SNPLIB_SRC_SNP_H_
