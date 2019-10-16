#include "rng.h"

#include <algorithm>

namespace rng {

RandomNumberGenerator::RandomNumberGenerator() {
  std::random_device rd;
  gen_ = std::mt19937(rd());
}

RandomNumberGenerator::~RandomNumberGenerator() {}

int RandomNumberGenerator::GenerateUniformInt(int min, int max) {
  std::uniform_int_distribution<> dis(min, max);
  return dis(gen_);
}

double_t RandomNumberGenerator::GenerateUniformDouble(double_t min,
                                                      double_t max) {
  std::uniform_real_distribution<> dis(min, max);
  return dis(gen_);
}

std::vector<size_t> RandomNumberGenerator::SampleRandomIndices(
    size_t start_idx, size_t end_idx, size_t num_samples) {
  size_t num_indices = end_idx - start_idx + 1;
  std::vector<size_t> indices(num_indices);
  std::iota(indices.begin(), indices.end(), start_idx);
  std::vector<size_t> ret;
  if (num_samples < num_indices) {
    std::shuffle(indices.begin(), indices.end(), gen_);
    ret = std::vector<size_t>(indices.begin(), indices.begin() + num_samples);
  } else {
    ret = indices;
  }
  return ret;
}

void RandomNumberGenerator::seed(uint_fast32_t seed) {
  seed_ = seed;
  gen_.seed(seed_);
}

uint_fast32_t RandomNumberGenerator::seed() { return seed_; }

}  // namespace rng
