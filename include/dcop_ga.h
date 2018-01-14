//
// Created by nick on 12/01/18.
//

#ifndef LWGA_DCOP_GA_H
#define LWGA_DCOP_GA_H

#include "ga_types.h"
#include "ga_utils.h"
namespace dcop_ga {

class Chromosome {
 public:
  Path path;
  Vector<double_t> angles;
  double fitness;
  double cost;
  void calculate_cost(Matrix<Matrix<double>> &cost_mat);
  void evaluate_chromosome(Matrix<Matrix<double>> &cost_mat,
                           std::vector<double> &rewards);
  void mutate(Matrix<Matrix<double>> &cost_mat,
              std::vector<double> &rewards,
              double max_cost,
              std::mt19937 &g);
};

Chromosome generate_chromosome (Matrix<double> &cost_mat,
                                double max_cost,
                                uint_fast32_t idx_start,
                                uint_fast32_t idx_finish,
                                std::mt19937 &g);

Chromosome tournament_select(std::vector<Chromosome> &population, uint_fast32_t tour_size, std::mt19937 &g);

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

Chromosome ga_cop(std::vector<std::vector<double> > &cost_mat,
                  std::vector<double> &rewards,
                  double max_cost,
                  uint_fast32_t idx_start,
                  uint_fast32_t idx_finish,
                  std::mt19937 &g);

}
#endif //LWGA_DCOP_GA_H
