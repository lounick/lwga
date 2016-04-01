//
// Created by nick on 30/03/16.
//

#include "ctop_ga.h"

Chromosome generate_chromosome(std::vector<std::vector<double> > &cost_mat,
                               std::vector<double> &max_cost_v,
                               uint idx_start,
                               uint idx_finish) {
  static std::random_device rd;
  static std::mt19937 g(rd());

  uint n_agents = max_cost_v.size();

  Chromosome c;
  c.genes.reserve(n_agents);

  std::vector<uint> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);

  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_start));
  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_finish));

  for (size_t i = 0; i < n_agents; i++) {
    Path path;
    path.reserve(cost_mat.size());

    bool done = false;
    double total_cost = 0;
    path.push_back(idx_start);
    while (!done) {
      size_t vsize = vertices.size();
      if (vsize == 0)
        break;
      size_t rand_idx = g() % vsize;
      uint next_vertex = vertices[rand_idx];
      double next_vertex_cost = cost_mat[path.back()][next_vertex] + 1;
      if (total_cost + next_vertex_cost + cost_mat[next_vertex][idx_finish] <= max_cost_v[i]) {
        total_cost += next_vertex_cost;
        path.push_back(next_vertex);
        vertices.erase(std::remove(vertices.begin(), vertices.end(), next_vertex));
      }
      else {
        done = true;
      }
    }
    path.push_back(idx_finish);
    Gene gene;
    gene.path = path;
    c.genes.push_back(gene);
  }
  return c;
}

Chromosome tournament_select(std::vector<Chromosome> &population, uint tour_size){
  static std::random_device rd;
  static std::mt19937 g(rd());

  std::vector<uint> indices(population.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), g);
  double max_fitness = DBL_MIN;
  Chromosome best_chromosome;

  for (size_t i = 0; i < tour_size; i++) {
    if (population[indices[i]].total_fitness > max_fitness) {
      max_fitness = population[indices[i]].total_fitness;
      best_chromosome = population[indices[i]];
    }
  }

  return best_chromosome;
}

std::pair<Chromosome, Chromosome> cx(Chromosome &c1,
                                     Chromosome &c2,
                                     std::vector<std::vector<double> > &cost_mat,
                                     std::vector<double> &max_cost_v) {
  /*
   * From each parent select the most fit path and exchange with the lowest fit path of the other parent.
   * Then check for double visits etc.
   * Then evaluate fitnesses and costs again.
   */
  std::sort(c1.genes.begin(), c1.genes.end(), [](Gene &g1, Gene &g2){ return g1.fitness > g2.fitness; });
  std::sort(c2.genes.begin(), c2.genes.end(), [](Gene &g1, Gene &g2){ return g1.fitness > g2.fitness; });
  Gene g1_best = c1.fitness
}

Chromosome mutate(Chromosome &c,
                  std::vector<std::vector<double> > &cost_mat,
                  std::vector<double> &rewards,
                  std::vector<double> &max_cost_v) {
  /*
   * Same as mutate for one vehicle but apply sequentially for many vehicles.
   */
}

std::pair<size_t, Chromosome> par_mutate(size_t idx,
                                         Chromosome c,
                                         std::vector<std::vector<double> > cost_mat,
                                         std::vector<double> rewards,
                                         std::vector<double> &max_cost_v){

}

double evaluate_chromosome(Chromosome &c, Matrix<double> &cost_mat, std::vector<double> &rewards){

}

double get_path_cost(Path &path, Matrix<double> &cost_mat){

}

std::vector<uint> two_opt_swap(std::vector<uint> &path, size_t &i, size_t &k){

}

std::pair<std::vector<uint>, double> two_opt(std::vector<uint> &path, Matrix<double> &cost_mat) {

}

std::pair<bool, double> check_feasibility(Chromosome &c, Matrix<double> &cost_mat, std::vector<double>& max_cost_v){

}

Chromosome ga_cop(Matrix<double> &cost_mat,
                  std::vector<double> &rewards,
                  std::vector<double> max_cost_v,
                  uint idx_start,
                  uint idx_finish){

}

std::vector<size_t> get_population_sample(size_t pop_size, int samples){

}
