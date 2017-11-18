//
// Created by nick on 31/03/16.
//

#ifndef LWGA_CTOP_GA_H
#define LWGA_CTOP_GA_H

#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cfloat>
#include <set>
#include <unordered_set>
#include <limits>
#include <chrono>
#include <thread>
#include <future>
#include "ga_types.h"
#include "ga_utils.h"

bool logically_equal(double a, double b, double error_factor = 1.0);

class Gene {
 public:
  Path path;
  double fitness;
  double cost;
  void evaluate_gene(const Matrix<double> &cost_mat, const std::vector<double> &rewards,
                     const std::vector<uint_fast32_t> &free_vertices);
  void calculate_cost(Matrix<double> &cost_mat);
  void mutate(Matrix<double> &cost_mat, std::vector<double> &rewards, double max_cost, std::mt19937 &g);
};

struct InsertMove {
  InsertMove(uint_fast32_t vertex,
             uint_fast32_t prev_vertex,
             uint_fast32_t next_vertex,
             double cost_increase,
             double total_reward,
             double heuristic);
  uint_fast32_t vertex;
  uint_fast32_t prev_vertex;
  uint_fast32_t next_vertex;
  double cost_increase;
  double total_reward;
  double heuristic;
};

struct Insertion {
  Insertion(uint_fast32_t vertex,
            uint_fast32_t prev_vertex,
            uint_fast32_t next_vertex,
            uint_fast32_t gene,
            double cost_increase,
            double total_reward,
            double heuristic);
  uint_fast32_t vertex;
  uint_fast32_t prev_vertex;
  uint_fast32_t next_vertex;
  uint_fast32_t gene;
  double cost_increase;
  double total_reward;
  double heuristic;
};

class Chromosome {
 public:
  Chromosome();
  Chromosome(size_t num_vertices, size_t num_genes, const uint_fast32_t &start_vertex,
             const uint_fast32_t &end_vertex);
  std::vector<Gene> genes;
  std::unordered_set<uint_fast32_t> seen_vertices;
  std::vector<uint_fast32_t> all_vertices;
  std::vector<uint_fast32_t> free_vertices;
  double total_fitness;
  void insertGene(Gene &gene);
  void removeCommonVertices(const Gene &gene);
  void evaluateGenes(const Matrix<double> &cost_mat, const std::vector<double> &rewards, const
                     std::vector<uint_fast32_t> free_vertices);
  void evaluate_chromosome(Matrix<double> &cost_mat,
                           std::vector<double> &rewards,
                           std::vector<double> &max_cost_v);
  void mutate(Matrix<double> &cost_mat, std::vector<double> &rewards, std::vector<double> &max_cost_v, std::mt19937 &g);
  void GenerateInsertMoves(const Matrix<double_t> &cost_mat,
                           const std::vector<double_t> &rewards,
                           const double_t &max_cost,
                           const Path &path, double_t path_cost,
                           std::vector<InsertMove> &moves,
                           double_t &min, double_t &max);
  double_t GenerateGRASPPath(const Matrix<double_t> &cost_mat,
                             const std::vector<double_t> &rewards,
                             const double_t &max_cost,
                             const uint_fast32_t &start_vertex,
                             const uint_fast32_t &end_vertex,
                             std::mt19937 &g, Path &path);
  double_t GenerateNNGRASPPath(const Matrix<double_t> &cost_mat,
                               const std::vector<double_t> &rewards,
                               const double_t &max_cost,
                               const uint_fast32_t &start_vertex,
                               const uint_fast32_t &end_vertex,
                               std::mt19937 &g, Path &path);
  double_t GenerateRandomPath(const Matrix<double_t> &cost_mat,
                              const std::vector<double_t> &rewards,
                              const double_t &max_cost,
                              const uint_fast32_t &start_vertex,
                              const uint_fast32_t &end_vertex,
                              std::mt19937 &g, Path &path);
  void GenerateGenes(Matrix<double_t> &cost_mat,
                     std::vector<double_t> &rewards,
                     std::vector<double_t> &max_cost_v,
                     const uint_fast32_t &start_vertex,
                     const uint_fast32_t &end_vertex);
 private:
  inline bool logically_equal(double a, double b, double error_factor = 1.0) {
    return a == b ||
        std::abs(a - b) < std::abs(std::min(a, b)) * std::numeric_limits<double>::epsilon() *
            error_factor;
  }
};

Chromosome generate_chromosome(Matrix<double> &cost_mat,
                               std::vector<double> &max_cost_v,
                               std::vector<double> &rewards,
                               uint idx_start,
                               uint idx_finish,
                               std::mt19937 &g,
                               std::string gen_method);

void expand_neighbours(std::vector<uint_fast32_t> &neighbours,
                       std::vector<uint_fast32_t> &checked,
                       const std::vector<uint_fast32_t> &vertices,
                       const Matrix<double> &cost_mat);

Chromosome tournament_select(std::vector<Chromosome> &population, uint tour_size, std::mt19937 &g);

void cx(Chromosome &c1,
        Chromosome &c2,
        Matrix<double> &cost_mat,
        std::vector<double> &max_cost_v,
        std::vector<double> &rewards,
        std::mt19937 &g);

void construct_offspring(Chromosome &offspring,
                         const Chromosome &parent1,
                         const Chromosome &parent2,
                         Matrix<double> &cost_mat,
                         std::vector<double> &max_cost_v,
                         std::vector<double> &rewards,
                         std::mt19937 &g);

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                std::vector<double> &max_cost_v);

void par_cx(std::vector<size_t> indices,
            std::vector<Chromosome> &pop,
            Matrix<double> &cost_mat,
            std::vector<double> &rewards,
            std::vector<double> &max_cost_v);

Chromosome ga_ctop(Matrix<double> &cost_mat,
                   std::vector<double> &rewards,
                   std::vector<double> max_cost_v,
                   uint idx_start,
                   uint idx_finish,
                   std::mt19937 &g,
                   size_t population_size,
                   size_t n_generations,
                   uint8_t tour_size,
                   double_t cx_rate,
                   double_t mutation_rate,
                   std::string generation_method,
                   double_t elite_percent);

#endif //LWGA_CTOP_GA_H
