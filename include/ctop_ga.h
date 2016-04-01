//
// Created by nick on 31/03/16.
//

#ifndef LWGA_CTOP_GA_H
#define LWGA_CTOP_GA_H

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
#include <include/ga_types.h>


class Gene {
 public:
  Path path;
  double fitness;
  double cost;
};

class Chromosome {
 public:
  std::vector<Gene> genes;
  double total_fitness;
};



Chromosome generate_chromosome(Matrix<double> &cost_mat,
                               std::vector<double> &max_cost_v,
                               uint idx_start,
                               uint idx_finish);

Chromosome tournament_select(std::vector<Chromosome> &population, uint tour_size = 3);

std::pair<Chromosome, Chromosome> cx(Chromosome &c1,
                                     Chromosome &c2,
                                     Matrix<double> &cost_mat,
                                     std::vector<double> &max_cost_v);

Chromosome mutate(Chromosome &c,
                  Matrix<double> &cost_mat,
                  std::vector<double> &rewards,
                  std::vector<double> &max_cost_v);

std::pair<size_t, Chromosome> par_mutate(size_t idx,
                                         Chromosome c,
                                         Matrix<double> cost_mat,
                                         std::vector<double> rewards,
                                         std::vector<double> &max_cost_v);

double evaluate_chromosome(Chromosome &c, Matrix<double> &cost_mat, std::vector<double> &rewards);

double get_path_cost(Path &path, Matrix<double> &cost_mat);

Path two_opt_swap(Path &path, size_t &i, size_t &k);

std::pair<Path, double> two_opt(Path &path, Matrix<double> &cost_mat);

std::pair<bool, double> check_feasibility(Chromosome &c, Matrix<double> &cost_mat, std::vector<double>& max_cost_v);

Chromosome ga_cop(Matrix<double> &cost_mat,
                  std::vector<double> &rewards,
                  std::vector<double> max_cost_v,
                  uint idx_start,
                  uint idx_finish);

std::vector<size_t> get_population_sample(size_t pop_size, int samples);

#endif //LWGA_CTOP_GA_H
