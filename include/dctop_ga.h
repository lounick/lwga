//
// Created by nick on 21/07/18.
//

#ifndef LWGA_DCTOP_GA_H
#define LWGA_DCTOP_GA_H

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

namespace dctop_ga{
class Chromosome {
 public:
  Vector<Path> paths;
  Vector<Path> angles;
  Vector<std::shared_ptr<Matrix<Matrix<double_t>>>> cost_matrices;
  Vector<double_t> radii;
  Vector<double_t> fitnesses;
  Vector<double_t> costs;
  std::unordered_set<uint_fast32_t> seen_vertices;
  std::vector<uint_fast32_t> all_vertices;
  std::vector<uint_fast32_t> free_vertices;
  void calculate_cost(size_t path_idx);
  void evaluate_chromosome(const Matrix<double_t> euclidean_cost,
                           const Vector<double_t> &rewards);
  void mutate(const Matrix<double_t> euclidean_cost,
              const Vector<double_t> std_angles,
              const Vector<double_t> &rewards,
              const Vector<double_t> max_costs,
              const Vector<size_t> start_idxs,
              std::mt19937 &g);
  void remove_vertex(const Matrix<double_t> euclidean_cost,
                     const Vector<double_t> std_angles,
                     const Vector<double_t> &rewards,
                     const Vector<double_t> max_costs,
                     const Vector<size_t> start_idxs,
                     std::mt19937 &g);
};

void expand_neighbours(std::vector<uint_fast32_t> &neighbours,
                       std::vector<uint_fast32_t> &checked,
                       const std::vector<uint_fast32_t> &vertices,
                       const Matrix<double> &cost_mat);

Chromosome generate_chromosome(std::shared_ptr<Vector<Point2D>> nodes,
                               const Vector<double_t> &std_angles,
                               const Vector<double_t> &radii,
                               const Vector<double_t> &max_costs,
                               const Vector<uint_fast32_t> &idx_start,
                               const Vector<uint_fast32_t> &idx_finish,
                               const Vector<std::shared_ptr<Matrix<Matrix<double_t>>>> &dubins_cost_mat,
                               const Matrix<double_t> &cost_mat,
                               std::mt19937 &g);

Chromosome extend_chromosome(std::shared_ptr<Vector<Point2D>> nodes,
                             const Vector<double_t> &std_angles,
                             const Vector<double_t> &radii,
                             const Vector<double_t> &max_costs,
                             const Vector<uint_fast32_t> &idx_start,
                             const Vector<uint_fast32_t> &idx_finish,
                             const Vector<std::shared_ptr<Matrix<Matrix<double_t>>>> &dubins_cost_mat,
                             const Matrix<double_t> &cost_mat,
                             std::mt19937 &g,
                             Chromosome init);

Chromosome tournament_select(std::vector<Chromosome> &population,
                             uint_fast32_t tour_size,
                             std::mt19937 &g);

void cx(Chromosome &c1,
        Chromosome &c2,
        const Vector<double_t> &std_angles,
        const Matrix<double_t> &cost_mat,
        const Vector<double_t> &rewards,
        const Vector<double_t> &max_costs,
        const Vector<uint_fast32_t> &start_idx,
        std::mt19937 &g);

void par_cx(const std::vector<size_t> indices,
            std::vector<Chromosome> &pop,
            const Vector<double_t> &std_angles,
            const Matrix<double_t> &cost_mat,
            const Vector<double_t> &rewards,
            const Vector<double_t> &max_costs,
            const Vector<uint_fast32_t> &start_idx,
            double_t max_cost);

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                const Vector<double_t> &std_angles,
                const Matrix<double_t> &cost_mat,
                const Vector<double_t> &rewards,
                const Vector<double_t> &max_costs,
                const Vector<uint_fast32_t> &start_idx,
                double &max_cost);

Chromosome ga_dctop(std::shared_ptr<Vector<Point2D>> nodes,
                    Vector<double_t> std_angles,
                    Vector<double_t> radii,
                    Vector<std::shared_ptr<Matrix<Matrix<double_t>>>> &dubins_cost_mat,
                    Matrix<double_t> &cost_mat,
                    std::vector<double_t> &rewards,
                    Vector<double_t> &max_costs,
                    Vector<uint_fast32_t> &idx_start,
                    Vector<uint_fast32_t> &idx_finish,
                    uint_fast16_t pop_size,
                    uint_fast8_t num_gen,
                    uint_fast8_t tour_size,
                    double_t cx_rate,
                    double_t mut_rate,
                    double_t elitist_rate,
                    std::mt19937 &g);

Chromosome replan_ga_dctop(Chromosome init,
                           std::shared_ptr<Vector<Point2D>> nodes,
                           Vector<double_t> std_angles,
                           Vector<double_t> radii,
                           Vector<std::shared_ptr<Matrix<Matrix<double_t>>>> &dubins_cost_mat,
                           Matrix<double_t> &cost_mat,
                           std::vector<double_t> &rewards,
                           Vector<double_t> &max_costs,
                           Vector<uint_fast32_t> &idx_start,
                           Vector<uint_fast32_t> &idx_finish,
                           uint_fast16_t pop_size,
                           uint_fast8_t num_gen,
                           uint_fast8_t tour_size,
                           double_t cx_rate,
                           double_t mut_rate,
                           double_t elitist_rate,
                           std::mt19937 &g);
}

#endif //LWGA_DCTOP_GA_H
