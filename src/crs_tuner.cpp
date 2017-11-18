//
// Created by nick on 17/11/17.
//

#include <cmath>
#include <cstdint>
#include <vector>
#include <utility>
#include <ctop_ga.h>
#include <tuning_instances.h>

template<typename T>
struct limits {
  T min;
  T max;
};

class Configuration {
 public:
  Configuration(uint_fast16_t pop_size,
                uint_fast8_t num_gen,
                uint_fast8_t tour_size,
                double_t cx_rate,
                double_t mut_rate,
                double_t elitist_rate)
      : pop_size_(pop_size),
        num_gen_(num_gen),
        tour_size_(tour_size),
        cx_rate_(cx_rate),
        mut_rate_(mut_rate),
        elitist_rate_(elitist_rate) {}
  Configuration() {}
  double_t rating() { return rating_; }
  void rating(double_t rating) { rating_ = rating; }
  double_t rating_deviation() { return rating_deviation_; }
  void rating_deviation(double_t rating_deviation) {
    rating_deviation_ = rating_deviation;
  }
  double_t rating_interval() { return rating_interval_; }
  void rating_interval(double_t rating_interval) {
    rating_interval_ = rating_interval;
  }
  uint_fast16_t pop_size() { return pop_size_; }
  void pop_size(uint_fast16_t pop_size) { pop_size_ = pop_size; }
  uint_fast8_t num_gen() { return num_gen_; }
  void num_gen(uint_fast8_t num_gen) { num_gen_ = num_gen; }
  uint_fast8_t tour_size() { return tour_size_; }
  void tour_size(uint_fast8_t tour_size) { tour_size_ = tour_size; }
  double_t cx_rate() { return cx_rate_; }
  void cx_rate(double_t cx_rate) { cx_rate_ = cx_rate; }
  double_t mut_rate() { return mut_rate_; }
  void mut_rate(double_t mut_rate) { mut_rate_ = mut_rate; }
  double_t elitist_rate() { return elitist_rate_; }
  void elitist_rate(double_t elitist_rate) { elitist_rate_ = elitist_rate; }
 private:
  double_t rating_;
  double_t rating_deviation_;
  double_t rating_interval_;
  uint_fast16_t pop_size_;
  uint_fast8_t num_gen_;
  uint_fast8_t tour_size_;
  double_t cx_rate_;
  double_t mut_rate_;
  double_t elitist_rate_;
};

using Configurations = std::vector<Configuration>;

struct Problem {
  Matrix<double_t> cost_mat;
  std::vector<double> rewards;
  std::vector<double> max_cost_v;
};

using Problems = std::vector<Problem>;

Configurations generatePopulation(uint_fast8_t pop_size,
                                   limits<uint_fast16_t> pop_size_lim,
                                   limits<uint_fast8_t> num_generations_lim,
                                   limits<uint_fast8_t> tour_size_lim,
                                   limits<double_t> cx_rate_lim,
                                   limits<double_t> mut_rate_lim,
                                   limits<double_t> elitist_rate_lim,
                                   std::mt19937 &g) {
  std::uniform_int_distribution<uint_fast16_t>
      pop_size_dis(pop_size_lim.min, pop_size_lim.max);
  std::uniform_int_distribution<uint_fast8_t>
      num_generations_dis(num_generations_lim.min, num_generations_lim.max);
  std::uniform_int_distribution<uint_fast8_t>
      tour_size_dis(tour_size_lim.min, tour_size_lim.max);
  std::uniform_real_distribution<double_t>
      cx_rate_dis(cx_rate_lim.min, cx_rate_lim.max);
  std::uniform_real_distribution<double_t>
      mut_rate_dis(mut_rate_lim.min, mut_rate_lim.max);
  std::uniform_real_distribution<double_t>
      elitist_rate_dis(elitist_rate_lim.min, elitist_rate_lim.max);

  Configurations pop;
  pop.reserve(pop_size);
  for (uint_fast8_t i = 0; i < pop_size; ++i) {
    pop.emplace_back(pop_size_dis(g),
                     num_generations_dis(g),
                     tour_size_dis(g),
                     cx_rate_dis(g),
                     mut_rate_dis(g),
                     elitist_rate_dis(g));
  }
}

Problems generateProblems(){
  Problems probs;
  std::vector<double_t> rewards5x5(canonical5_costmat.size(),1);
  rewards5x5.front() = 0;
  rewards5x5.back() = 0;
  std::vector<double_t> rewards7x7(canonical7_costmat.size(),1);
  rewards7x7.front() = 0;
  rewards7x7.back() = 0;
  std::vector<double_t> rewards9x9(canonical9_costmat.size(),1);
  rewards9x9.front() = 0;
  rewards9x9.back() = 0;
  Problem prob5x5canonical3robs = {canonical5_costmat, rewards5x5, std::vector<double>(3, 3*((canonical5_objective+25.0)/3)/4.0)};
  Problem prob5x5random3robs = {random5_cost_mat, rewards5x5, std::vector<double>(3, 3*((random5_objective+25.0)/3)/4.0)};
  Problem prob5x5canonical4robs = {canonical5_costmat, rewards5x5, std::vector<double>(4, 3*((canonical5_objective+25.0)/4)/4.0)};
  Problem prob5x5random4robs = {random5_cost_mat, rewards5x5, std::vector<double>(4, 3*((random5_objective+25.0)/4)/4.0)};
  Problem prob5x5canonical5robs = {canonical5_costmat, rewards5x5, std::vector<double>(5, 3*((canonical5_objective+25.0)/5)/4.0)};
  Problem prob5x5random5robs = {random5_cost_mat, rewards5x5, std::vector<double>(5, 3*((random5_objective+25.0)/5)/4.0)};
  probs.push_back(std::move(prob5x5canonical3robs));
  probs.push_back(std::move(prob5x5random3robs));
  probs.push_back(std::move(prob5x5canonical4robs));
  probs.push_back(std::move(prob5x5random4robs));
  probs.push_back(std::move(prob5x5canonical5robs));
  probs.push_back(std::move(prob5x5random5robs));
  Problem prob7x7canonical3robs = {canonical7_costmat, rewards7x7, std::vector<double>(3, 3*((canonical7_objective+49.0)/3)/4.0)};
  Problem prob7x7random3robs = {random7_cost_mat, rewards7x7, std::vector<double>(3, 3*((random7_objective+49.0)/3)/4.0)};
  Problem prob7x7canonical4robs = {canonical7_costmat, rewards7x7, std::vector<double>(4, 3*((canonical7_objective+49.0)/4)/4.0)};
  Problem prob7x7random4robs = {random7_cost_mat, rewards7x7, std::vector<double>(4, 3*((random7_objective+49.0)/4)/4.0)};
  Problem prob7x7canonical5robs = {canonical7_costmat, rewards7x7, std::vector<double>(5, 3*((canonical7_objective+49.0)/5)/4.0)};
  Problem prob7x7random5robs = {random7_cost_mat, rewards7x7, std::vector<double>(5, 3*((random7_objective+49.0)/5)/4.0)};
  probs.push_back(std::move(prob7x7canonical3robs));
  probs.push_back(std::move(prob7x7random3robs));
  probs.push_back(std::move(prob7x7canonical4robs));
  probs.push_back(std::move(prob7x7random4robs));
  probs.push_back(std::move(prob7x7canonical5robs));
  probs.push_back(std::move(prob7x7random5robs));
  Problem prob9x9canonical3robs = {canonical9_costmat, rewards9x9, std::vector<double>(3, 3*((canonical9_objective+81.0)/3)/4.0)};
  Problem prob9x9random3robs = {random9_cost_mat, rewards9x9, std::vector<double>(3, 3*((random9_objective+81.0)/3)/4.0)};
  Problem prob9x9canonical4robs = {canonical9_costmat, rewards9x9, std::vector<double>(4, 3*((canonical9_objective+81.0)/4)/4.0)};
  Problem prob9x9random4robs = {random9_cost_mat, rewards9x9, std::vector<double>(4, 3*((random9_objective+81.0)/4)/4.0)};
  Problem prob9x9canonical5robs = {canonical9_costmat, rewards9x9, std::vector<double>(5, 3*((canonical9_objective+81.0)/5)/4.0)};
  Problem prob9x9random5robs = {random9_cost_mat, rewards9x9, std::vector<double>(5, 3*((random9_objective+81.0)/5)/4.0)};
  probs.push_back(std::move(prob9x9canonical3robs));
  probs.push_back(std::move(prob9x9random3robs));
  probs.push_back(std::move(prob9x9canonical4robs));
  probs.push_back(std::move(prob9x9random4robs));
  probs.push_back(std::move(prob9x9canonical5robs));
  probs.push_back(std::move(prob9x9random5robs));
  return probs;
}

int main (int argc, char *argv[]){
  /*
   * 1. Generate all the problem instances.
   * 2. Generate population.
   * 3. Optimise.
   */
  uint_fast8_t pop_size = 100;
  limits<uint_fast16_t> pop_size_lim = {10, 1000};
  limits<uint_fast8_t> num_gen_lim = {5,250};
  limits<uint_fast8_t> tour_size_lim = {3,10};
  limits<double_t> cx_lim = {0.0,1.0};
  limits<double_t> mut_lim = {0.0,1.0};
  limits<double_t> elites_lim = {0.01,0.2};
  std::random_device rd;
  std::mt19937 g(rd());
  generatePopulation(
      pop_size, pop_size_lim, num_gen_lim, tour_size_lim, cx_lim, mut_lim,
      elites_lim, g);
  generateProblems();
  uint8_t max_exp = 10;
  for (uint8_t exp = 0; exp < max_exp; ++exp){

  }
}