#include "rng.h"
#include "gtest/gtest.h"

TEST(RNGTest, setSeedTest) {
  rng::RandomNumberGenerator gen;
  uint_fast32_t seed = 42;
  gen.seed(seed);
  EXPECT_EQ(gen.seed(), seed);
}

TEST(RNGTest, generateUniformIntegerTest){
  rng::RandomNumberGenerator gen;
  int generated_number = gen.GenerateUniformInt(0, 10);
  EXPECT_GE(generated_number, 0);
  EXPECT_LE(generated_number, 10);
}

TEST(RNGTest, generateUniformDoubleTest){
  rng::RandomNumberGenerator gen;
  int generated_number = gen.GenerateUniformDouble(0.0, 1.0);
  EXPECT_GE(generated_number, 0.0);
  EXPECT_LE(generated_number, 1.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
