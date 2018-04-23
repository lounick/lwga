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
  Vector<uint_fast32_t> angles;
  std::shared_ptr<Vector<Point2D>> nodes;
  double_t rho;
  double_t fitness;
  double_t cost;
  std::unordered_set<uint_fast32_t> seen_vertices;
  std::vector<uint_fast32_t> all_vertices;
  std::vector<uint_fast32_t> free_vertices;
  void calculate_cost(Matrix<Matrix<double_t>>&dubins_cost_mat);
  void evaluate_chromosome(Matrix<Matrix<double_t>>&dubins_cost_mat,
                           Matrix<double_t> &cost_mat,
                           std::vector<double_t> &rewards);
  void mutate(Matrix<Matrix<double_t>>&dubins_cost_mat,
              Matrix<double_t> &cost_mat,
              Vector<double_t> &std_angles,
              std::vector<double_t> &rewards,
              double_t max_cost,
              std::mt19937 &g, size_t start_idx);
  void remove_vertex(Matrix<Matrix<double_t>>&dubins_cost_mat,
                     Matrix<double_t> &cost_mat,
                     Vector<double_t> &std_angles,
                     std::vector<double_t> &rewards,
                     double_t max_cost,
                     std::mt19937 &g, size_t start_idx);
  double_t max_dist = 0.0; //TODO: [OP] This should be set when chromosome is generated
};

void expand_neighbours(std::vector<uint_fast32_t> &neighbours,
                       std::vector<uint_fast32_t> &checked,
                       const std::vector<uint_fast32_t> &vertices,
                       const Matrix<double> &cost_mat);

Chromosome generate_chromosome(std::shared_ptr<Vector<Point2D>> nodes,
                               const Vector<double_t> &std_angles,
                               double_t rho,
                               double_t max_cost,
                               uint_fast32_t idx_start,
                               uint_fast32_t idx_finish,
                               const Matrix<Matrix<double_t>>&dubins_cost_mat,
                               const Matrix<double_t> &cost_mat,
                               std::mt19937 &g);

Chromosome extend_chromosome(std::shared_ptr<Vector<Point2D>> nodes,
                             const Vector<double_t> &std_angles,
                             double_t rho,
                             double_t max_cost,
                             uint_fast32_t idx_start,
                             uint_fast32_t idx_finish,
                             const Matrix<Matrix<double_t>>&dubins_cost_mat,
                             const Matrix<double_t> &cost_mat,
                             std::mt19937 &g,
                             Chromosome init);

Chromosome tournament_select(std::vector<Chromosome> &population,
                             uint_fast32_t tour_size,
                             std::mt19937 &g);

void cx(Chromosome &c1,
        Chromosome &c2,
        Matrix<Matrix<double_t>>&dubins_cost_mat,
        Vector<double_t> &std_angles,
        Matrix<double_t> &cost_mat,
        Vector<double_t> &rewards,
        double max_cost,
        std::mt19937 &g, size_t start_idx);

void par_cx(std::vector<size_t> indices,
            std::vector<Chromosome> &pop,
            Matrix<Matrix<double_t>>&dubins_cost_mat,
            Vector<double_t> &std_angles,
            Matrix<double_t> &cost_mat,
            std::vector<double_t> &rewards,
            double_t max_cost, size_t start_idx);

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<Matrix<double_t>>&dubins_cost_mat,
                Matrix<double> &cost_mat,
                Vector<double_t> &std_angles,
                std::vector<double> &rewards,
                double &max_cost, size_t start_idx);

Chromosome ga_dcop(std::shared_ptr<Vector<Point2D>> nodes,
                   Vector<double_t> std_angles,
                   double_t rho,
                   Matrix<Matrix<double_t>>&dubins_cost_mat,
                   Matrix<double_t> &cost_mat,
                   std::vector<double_t> &rewards,
                   double_t max_cost,
                   uint_fast32_t idx_start,
                   uint_fast32_t idx_finish,
                   uint_fast16_t pop_size,
                   uint_fast8_t num_gen,
                   uint_fast8_t tour_size,
                   double_t cx_rate,
                   double_t mut_rate,
                   double_t elitist_rate,
                   std::mt19937 &g);

Chromosome replan_ga_dcop(std::shared_ptr<Vector<Point2D>> nodes,
                          Vector<double_t> std_angles,
                          double_t rho,
                          Matrix<Matrix<double_t>>&dubins_cost_mat,
                          Matrix<double_t> &cost_mat,
                          std::vector<double_t> &rewards,
                          double_t max_cost,
                          uint_fast32_t idx_start,
                          uint_fast32_t idx_finish,
                          uint_fast16_t pop_size,
                          uint_fast8_t num_gen,
                          uint_fast8_t tour_size,
                          double_t cx_rate,
                          double_t mut_rate,
                          double_t elitist_rate,
                          std::mt19937 &g, Chromosome init);

}
#endif //LWGA_DCOP_GA_H
