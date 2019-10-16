#ifndef RNG_H_
#define RNG_H_

#include <random>
#include <vector>

namespace rng {
class RandomNumberGenerator {
 public:
  RandomNumberGenerator();
  virtual ~RandomNumberGenerator();
  virtual int GenerateUniformInt(int min, int max);
  virtual double_t GenerateUniformDouble(double_t min, double_t max);
  virtual std::vector<size_t> SampleRandomIndices(size_t start_idx,
                                                  size_t end_idx,
                                                  size_t num_samples);
  void seed(uint_fast32_t seed);
  uint_fast32_t seed();

 protected:
  uint_fast32_t seed_;
  std::mt19937 gen_;
};
}  // namespace rng

#endif
