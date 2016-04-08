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
#include "ga_types.h"
#include "ga_utils.h"


class Gene {
 public:
  Path path;
  double fitness;
  double cost;
  void evaluate_gene(Matrix<double> &cost_mat, std::vector<double> &rewards, std::vector<uint_fast32_t> &free_vertices);
  void calculate_cost(Matrix<double> &cost_mat);
  void mutate(Matrix<double> &cost_mat, std::vector<double> &rewards, double max_cost, std::mt19937 &g);
};

class Chromosome {
 public:
  std::vector<Gene> genes;
  double total_fitness;
  void evaluate_chromosome();
  void mutate(Matrix<double> &cost_mat, std::vector<double> &rewards, std::vector<double> &max_cost_v, std::mt19937 &g);
};



Chromosome generate_chromosome (Matrix<double> &cost_mat,
                                std::vector<double> &max_cost_v,
                                uint idx_start,
                                uint idx_finish,
                                std::mt19937 &g);

Chromosome tournament_select(std::vector<Chromosome> &population, uint tour_size, std::mt19937 &g);

void cx(Chromosome &c1, Chromosome &c2, Matrix<double> &cost_mat, std::vector<double> &max_cost_v, std::vector<double> &rewards);

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                std::vector<double> &max_cost_v);

Chromosome ga_ctop(Matrix<double> &cost_mat,
                  std::vector<double> &rewards,
                  std::vector<double> max_cost_v,
                  uint idx_start,
                  uint idx_finish,
                  std::mt19937 &g);

#endif //LWGA_CTOP_GA_H
