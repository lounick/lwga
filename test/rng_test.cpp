#include "rng.h"
#include "gtest/gtest.h"

TEST(RNGTest, setSeedTest) {
  rng::RandomNumberGenerator gen;
  uint_fast32_t seed = 42;
  gen.seed(seed);
  EXPECT_EQ(gen.seed(), seed);
}

TEST(RNGTest, generateUniformIntegerTest) {
  rng::RandomNumberGenerator gen;
  int generated_number = gen.GenerateUniformInt(0, 10);
  EXPECT_GE(generated_number, 0);
  EXPECT_LE(generated_number, 10);
}

TEST(RNGTest, generateUniformDoubleTest) {
  rng::RandomNumberGenerator gen;
  int generated_number = gen.GenerateUniformDouble(0.0, 1.0);
  EXPECT_GE(generated_number, 0.0);
  EXPECT_LE(generated_number, 1.0);
}

TEST(RNGTest, sampleRandomIndicesTest) {
  rng::RandomNumberGenerator gen;
  size_t start_idx = 0;
  size_t end_idx = 9;
  size_t num_samples = 3;
  std::vector<size_t> idxs =
      gen.SampleRandomIndices(start_idx, end_idx, num_samples);
  EXPECT_EQ(idxs.size(), num_samples);
  for (size_t idx : idxs) {
    EXPECT_GE(idx, start_idx);
    EXPECT_LE(idx, end_idx);
  }
}

TEST(RNGTest, sampleRandomIndicesiTooManySamplesTest) {
  rng::RandomNumberGenerator gen;
  size_t start_idx = 0;
  size_t end_idx = 9;
  size_t num_samples = 20;
  std::vector<size_t> idxs =
      gen.SampleRandomIndices(start_idx, end_idx, num_samples);
  EXPECT_EQ(idxs.size(), end_idx - start_idx + 1);
  for (size_t idx : idxs) {
    EXPECT_GE(idx, start_idx);
    EXPECT_LE(idx, end_idx);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
