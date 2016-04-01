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

class Chromosome {
 public:
  Path path;
  double fitness;
  double cost;
};

Chromosome generate_chromosome(Matrix<double> &cost_mat, double max_cost, uint_fast32_t idx_start, uint_fast32_t idx_finish);

Chromosome tournament_select(std::vector<Chromosome> &population, uint_fast32_t tour_size = 3);

std::pair<Chromosome, Chromosome> cx(Chromosome &c1, Chromosome &c2, Matrix<double> &cost_mat, double max_cost);

Chromosome mutate(Chromosome &c, Matrix<double> &cost_mat, std::vector<double> &rewards, double max_cost);

std::pair<size_t, Chromosome> par_mutate(size_t idx,
                                         Chromosome c,
                                         Matrix<double> cost_mat,
                                         std::vector<double> rewards,
                                         double max_cost);

double evaluate_chromosome(Chromosome &c, Matrix<double> &cost_mat, std::vector<double> &rewards);

double get_path_cost(Path &path, Matrix<double> &cost_mat);

Path two_opt_swap(Path &path, size_t &i, size_t &k);

std::pair<Path, double> two_opt(Path &path, Matrix<double> &cost_mat);

std::pair<bool, double> check_feasibility(Chromosome &c, Matrix<double> &cost_mat, double max_cost);

Chromosome ga_cop(std::vector<std::vector<double> > &cost_mat,
                  std::vector<double> &rewards,
                  double max_cost,
                  uint_fast32_t idx_start,
                  uint_fast32_t idx_finish);

std::vector<size_t> get_population_sample(size_t pop_size, int samples);

#endif //GUROBI_TESTS_COP_GA_H
