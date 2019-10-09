#include "rng.h"

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

void RandomNumberGenerator::seed(uint_fast32_t seed) {
  seed_ = seed;
  gen_.seed(seed_);
}

uint_fast32_t RandomNumberGenerator::seed() { return seed_; }

}  // namespace rng
