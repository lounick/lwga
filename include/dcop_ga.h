//
// Created by nick on 12/01/18.
//

#ifndef LWGA_DCOP_GA_H
#define LWGA_DCOP_GA_H

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
#include "dubins.h"
#include "ga_types.h"
#include "ga_utils.h"

namespace dcop_ga {

class Chromosome {
 public:
  Path path;
  Vector<double_t> angles;
  std::shared_ptr<Vector<Point2D>> nodes;
  double_t rho;
  double_t fitness;
  double_t cost;
  std::unordered_set<uint_fast32_t> seen_vertices;
  std::vector<uint_fast32_t> all_vertices;
  std::vector<uint_fast32_t> free_vertices;
  void calculate_cost();
  void evaluate_chromosome(Matrix<double_t> &cost_mat,
                           std::vector<double_t> &rewards);
  void mutate(Matrix<double_t> &cost_mat,
              std::vector<double_t> &rewards,
              double_t max_cost,
              std::mt19937 &g);
};

void expand_neighbours(std::vector<uint_fast32_t> &neighbours,
                       std::vector<uint_fast32_t> &checked,
                       const std::vector<uint_fast32_t> &vertices,
                       const Matrix<double> &cost_mat);

Chromosome generate_chromosome(std::shared_ptr<Vector<Point2D>> nodes,
                               Vector<double_t> std_angles,
                               double_t rho,
                               double_t max_cost,
                               uint_fast32_t idx_start,
                               uint_fast32_t idx_finish,
                               const Matrix<double_t> &cost_mat,
                               std::mt19937 &g);

Chromosome tournament_select(std::vector<Chromosome> &population,
                             uint_fast32_t tour_size,
                             std::mt19937 &g);

std::pair<Chromosome, Chromosome> cx(Chromosome &c1,
                                     Chromosome &c2,
                                     Matrix<double> &cost_mat,
                                     double max_cost,
                                     std::mt19937 &g);

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                double &max_cost);

Chromosome ga_dcop(std::shared_ptr<Vector<Point2D>> nodes,
                   Vector<double_t> std_angles,
                   double_t rho,
                   Matrix<double_t> &cost_mat,
                   std::vector<double_t> &rewards,
                   double_t max_cost,
                   uint_fast32_t idx_start,
                   uint_fast32_t idx_finish,
                   std::mt19937 &g);

}
#endif //LWGA_DCOP_GA_H
