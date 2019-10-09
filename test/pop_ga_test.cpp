#include "pop_ga.h"
#include "rng.h"

#include "gmock/gmock.h"
//#include "gtest/gtest.h"

class MockRNG : public rng::RandomNumberGenerator {
 public:
  MOCK_METHOD(int, GenerateUniformInt, (int min, int max), (override));
  MOCK_METHOD(double_t, GenerateUniformDouble, (double_t min, double_t max),
              (override));
};

TEST(POPGATest, totalRewardTest) {
  Path p{0, 1, 5};
  Vector<double_t> r{0.0, 10.0, 10.0, 10.0, 10.0, 0.0};
  const double_t expected_total_reward = 10.0;
  EXPECT_EQ(pop_ga::CalculateTotalReward(p, r), expected_total_reward);
}

TEST(POPGATest, actualRewardTest) {
  Path p{0, 1, 2, 5};
  Vector<double_t> r{0.0, 10.0, 10.0, 10.0, 10.0, 0.0};
  Vector<double_t> s{1.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const double_t expected_actual_reward = 10.0;
  EXPECT_EQ(pop_ga::CalculateActualReward(p, r, s), expected_actual_reward);
}

TEST(POPGATest, expectedRewardTest) {
  Path p{0, 1, 2, 5};
  Vector<double_t> r{0.0, 10.0, 10.0, 10.0, 10.0, 0.0};
  Vector<double_t> probs{1.0, 0.5, 0.5, 0.5, 0.5, 1.0};
  const double_t expected_total_reward = 10.0;
  EXPECT_EQ(pop_ga::CalculateExpectedReward(p, r, probs),
            expected_total_reward);
}

TEST(POPGATest, totalCostTest) {
  Path p{0, 1, 5};
  Matrix<double_t> c{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};
  const double_t expected_total_cost = 4.0;
  EXPECT_EQ(pop_ga::CalculateTotalCost(p, c), expected_total_cost);
}

TEST(POPGATest, actualCostTest) {
  Path p{0, 1, 2, 5};
  Matrix<double_t> c{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};
  Vector<double_t> s{1.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const double_t expected_actual_cost = 4.0;
  EXPECT_EQ(pop_ga::CalculateActualCost(p, c, s), expected_actual_cost);
}

TEST(POPGATest, expectedCostTest) {
  Path p{0, 1, 2, 5};
  Matrix<double_t> c{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};
  Vector<double_t> probs{1.0, 0.5, 0.5, 0.5, 0.5, 1.0};
  const double_t expected_total_cost = 5.0;
  EXPECT_EQ(pop_ga::CalculateExpectedCost(p, c, probs), expected_total_cost);
}

TEST(POPGATest, tmaxObjectiveTest) {
  Path p{0, 1, 3, 5};
  Vector<double_t> r{0.0, 10.0, 10.0, 10.0, 10.0, 0.0};
  Vector<double_t> probs{1.0, 0.5, 0.5, 0.5, 0.5, 1.0};
  Matrix<double_t> c{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};
  double_t cost_per_time_unit = 1;
  double_t expected_objective = -4.0;
  EXPECT_EQ(pop_ga::CalculateTMaxObjective(p, r, probs, c, cost_per_time_unit),
            expected_objective);
}

TEST(POPGATest, expectedTimeObjectiveTest) {
  Path p{0, 1, 2, 5};
  Vector<double_t> r{0.0, 10.0, 10.0, 10.0, 10.0, 0.0};
  Vector<double_t> probs{1.0, 0.5, 0.5, 0.5, 0.5, 1.0};
  Matrix<double_t> c{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};
  double_t cost_per_time_unit = 1;
  double_t expected_objective = 5.0;
  EXPECT_EQ(pop_ga::CalculateExpectedTimeObjective(p, r, probs, c,
                                                   cost_per_time_unit),
            expected_objective);
}

TEST(POPGATest, generateRandomChromosomeCostTest) {
  pop_ga::Properties props;
  props.generation_method = pop_ga::GenerationMethod::RANDOM;
  props.start_id = 0;
  props.end_id = 5;
  props.maximum_cost = 7.0;
  props.cost_per_time_unit = 1.0;

  Matrix<double_t> costs{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};

  // TODO: REPLACE WITH MOCK
  // rng::RandomNumberGenerator g;
  // g.seed(5);
  MockRNG g;
  EXPECT_CALL(g, GenerateUniformInt(::testing::_, ::testing::_))
      .Times(2)
      .WillOnce(::testing::Return(0))
      .WillOnce(::testing::ReturnArg<1>());
  pop_ga::Chromosome c = pop_ga::GenerateChromosome(props, costs, g);
  double_t path_cost = 4.0;
  EXPECT_DOUBLE_EQ(c.cost, path_cost);
  EXPECT_EQ(c.p.front(), props.start_id);
  EXPECT_EQ(c.p[1], 1);
  EXPECT_EQ(c.p.back(), props.end_id);
}

TEST(POPGATest, generateRandomChromosomeCostLETest) {
  pop_ga::Properties props;
  props.generation_method = pop_ga::GenerationMethod::RANDOM;
  props.start_id = 0;
  props.end_id = 5;
  props.maximum_cost = 7.0;
  props.cost_per_time_unit = 1.0;

  Matrix<double_t> costs{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};

  rng::RandomNumberGenerator g;
  size_t num_executions = 100000;
  for (size_t i = 0; i < num_executions; ++i) {
    pop_ga::Chromosome c = pop_ga::GenerateChromosome(props, costs, g);
    EXPECT_LE(c.cost, props.maximum_cost);
    EXPECT_EQ(c.p.front(), props.start_id);
    EXPECT_EQ(c.p.back(), props.end_id);
  }
}

TEST(POPGATest, generateGRASPChromosomeTest) {
  pop_ga::Properties props;
  props.generation_method = pop_ga::GenerationMethod::GRASP;
  props.start_id = 0;
  props.end_id = 5;
  props.maximum_cost = 7.0;
  props.cost_per_time_unit = 1.0;

  Matrix<double_t> costs{
      {0.0, 2.0, 3.0, 7.0, 9.0, 0.0},   {2.0, 0.0, 5.0, 5.0, 7.0, 2.0},
      {3.0, 5.0, 0.0, 10.0, 12.0, 3.0}, {7.0, 5.0, 10.0, 0.0, 2.0, 7.0},
      {9.0, 7.0, 12.0, 2.0, 0.0, 9.0},  {0.0, 2.0, 3.0, 7.0, 9.0, 0.0}};

  // TODO: REPLACE WITH MOCK
  rng::RandomNumberGenerator g;
  pop_ga::Chromosome c = pop_ga::GenerateChromosome(props, costs, g);
  EXPECT_LE(c.cost, props.maximum_cost);
  EXPECT_GT(c.cost, 0);
  //  EXPECT_TRUE((c.cost > 0.0) && (c.cost <= props.maximum_cost));
  //using ::testing::AllOf;
  //using ::testing::Gt;
  //using ::testing::Le;
  //EXPECT_THAT(c.cost, AllOf(Gt(0), Le(props.maximum_cost)));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleMock(&argc, argv);
  return RUN_ALL_TESTS();
}
