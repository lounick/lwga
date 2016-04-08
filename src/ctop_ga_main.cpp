//
// Created by nick on 31/03/16.
//

#include "ctop_ga.h"

int main(){
  std::vector<std::pair<double, double> > nodes;
  std::vector<std::vector<double> > cost_mat;
  nodes.push_back(std::make_pair(0, 0));
  for (int i = 1; i < 10; i++) {
    for (int j = -4; j < 5; j++) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(10, 0));

  for (int i = 0; i < 83; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 83; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 83; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  std::vector<double> rewards(cost_mat.size(), 1);
  rewards[0] = 0;
  rewards[82] = 0;
  uint_fast32_t num_robots = 4;
  std::vector<double> max_cost_v(num_robots, 27);
  std::vector<double> fitnesses;
  std::vector<double> times;

  int nexp = 20;

  std::random_device rd;
  std::mt19937 g(rd());

  for (int exp = 0; exp < nexp; exp++) {
    auto start = std::chrono::high_resolution_clock::now();
    Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 82, g);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());

    std::unordered_set<uint_fast32_t> seen;
    for(uint_fast32_t robot = 0; robot < num_robots; ++robot){
      for(size_t path_idx = 0; path_idx < c.genes[robot].path.size(); ++path_idx){
        seen.insert(c.genes[robot].path[path_idx]);
      }
    }

    double fitness = 0;
    std::vector<uint_fast32_t> vertices(cost_mat.size());
    std::iota(vertices.begin(), vertices.end(), 0);
    std::vector<uint_fast32_t> free_vertices;
    std::vector<uint_fast32_t> visited_vertices(seen.begin(), seen.end());
    std::sort(visited_vertices.begin(), visited_vertices.end());
    std::set_difference(vertices.begin(),
                        vertices.end(),
                        visited_vertices.begin(),
                        visited_vertices.end(),
                        std::back_inserter(free_vertices));

    seen.clear();
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
    for(uint_fast32_t robot = 0; robot < num_robots; ++robot) {
      for (size_t i = 1; i < c.genes[robot].path.size() - 1; i++) {
        double extras = 0;
        insert_ret = seen.insert(c.genes[robot].path[i]);
        if (insert_ret.second) {
          for (size_t j = 0; j < free_vertices.size(); j++) {
            if (cost_mat[c.genes[robot].path[i]][free_vertices[j]] < 2) {
              extras += std::exp(-2 * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
            }
          }
          fitness += rewards[c.genes[robot].path[i]] + extras;
        }
      }
    }
    c.evaluate_chromosome();
    std::cout << c.total_fitness << " " << fitness << std::endl;
    fitnesses.push_back(fitness);
  }
  std::cout << "Average fitness: " << std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size()
      << std::endl;
  std::cout << "Average time: " << std::accumulate(times.begin(), times.end(), 0.0) / times.size() << std::endl;
  return 0;
}