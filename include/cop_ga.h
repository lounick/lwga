//
// Created by nick on 31/03/16.
//

#ifndef GUROBI_TESTS_COP_GA_H
#define GUROBI_TESTS_COP_GA_H

#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include <cfloat>
#include <set>
#include <unordered_set>
#include <limits>
#include <chrono>
#include <thread>
#include <future>
#include "ga_types.h"
#include "ga_utils.h"
namespace cop_ga {
class Chromosome {
 public:
  Path path;
  double fitness;
  double cost;
  void calculate_cost(Matrix<double> &cost_mat);
  void evaluate_chromosome(Matrix<double> &cost_mat,
                           std::vector<double> &rewards);
  void mutate(Matrix<double> &cost_mat,
              std::vector<double> &rewards,
              double max_cost,
              std::mt19937 &g);
};

Chromosome generate_chromosome(Matrix<double> &cost_mat,
                               double max_cost,
                               uint_fast32_t idx_start,
                               uint_fast32_t idx_finish,
                               std::mt19937 &g);

Chromosome tournament_select(std::vector<Chromosome> &population,
                             uint_fast32_t tour_size,
                             std::mt19937 &g);

std::pair<Chromosome, Chromosome> cx(Chromosome &c1,
                                     Chromosome &c2,
                                     Matrix<double> &cost_mat,
                                     double max_cost,
                                     std::mt19937 &g);

void cxv(Chromosome &c1,
         Chromosome &c2,
         Matrix<double> &cost_mat,
         std::vector<double> &rewards,
         double max_cost,
         std::mt19937 &g);

void par_cx(std::vector<size_t> indices,
            std::vector<Chromosome> &pop,
            Matrix<double> &cost_mat,
            std::vector<double> &rewards,
            double &max_cost);

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                double &max_cost);

Chromosome ga_cop(std::vector<std::vector<double> > &cost_mat,
                  std::vector<double> &rewards,
                  double max_cost,
                  uint_fast32_t idx_start,
                  uint_fast32_t idx_finish,
                  std::mt19937 &g,
                  size_t population_size,
                  size_t n_generations,
                  uint8_t tour_size,
                  double_t cx_rate,
                  double_t mutation_rate,
                  double_t elite_percent);
}
#endif //GUROBI_TESTS_COP_GA_H
