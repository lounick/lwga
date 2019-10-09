//
// Created by nick on 31/03/16.
//

#include "cop_ga.h"
#include "ga_utils.h"

int main() {
  std::vector<std::pair<double, double> > nodes;
  std::vector<std::vector<double> > cost_mat;
  nodes.push_back(std::make_pair(0, 0));
  for (int i = 1; i < 6; i++) {
    for (int j = -2; j < 3; j++) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(6, 0));

  for (int i = 0; i < 27; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 27; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 27; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  std::vector<double> rewards(cost_mat.size(), 1);
  rewards[0] = 0;
  rewards[26] = 0;
  std::vector<double> fitnesses;
  std::vector<double> times;

  int nexp = 10;

  std::random_device rd;
  std::mt19937 g(rd());

  for (int exp = 0; exp < nexp; exp++) {
    auto start = std::chrono::high_resolution_clock::now();
    cop_ga::Chromosome c = cop_ga::ga_cop(cost_mat, rewards, 51*3/4, 0, 26, g, 100, 50, 3, 0.1, 0.8, 0.05);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());
//    for (uint_fast32_t vertex:c.path) {
//      std::cout << vertex << " ";
//    }
//    std::cout << std::endl;

    double fitness = 0;
    std::vector<uint_fast32_t> vertices(cost_mat.size());
    std::iota(vertices.begin(), vertices.end(), 0);
    std::vector<uint_fast32_t> free_vertices;
    std::vector<uint_fast32_t> visited_vertices = c.path;
    std::sort(visited_vertices.begin(), visited_vertices.end());
    std::set_difference(vertices.begin(),
                        vertices.end(),
                        visited_vertices.begin(),
                        visited_vertices.end(),
                        std::back_inserter(free_vertices));
    std::unordered_set<uint_fast32_t> seen;
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
    for (size_t i = 1; i < c.path.size() - 1; i++) {
      double extras = 0;
      insert_ret = seen.insert(c.path[i]);
      if (insert_ret.second) {
        for (size_t j = 0; j < free_vertices.size(); j++) {
          if (cost_mat[c.path[i]][free_vertices[j]] < 2) {
            extras += std::exp((log(0.01)/2) * cost_mat[c.path[i]][free_vertices[j]]);
          }
        }
        fitness += rewards[c.path[i]] + extras;
      }
    }
    fitnesses.push_back(fitness);
    std::cout << fitness << std::endl;
    std::cout << get_path_cost(c.path, cost_mat) << std::endl;
    print_path(c.path);
//    std::cout << fitness << std::endl;
//    std::cout << get_path_cost(c.path, cost_mat) << std::endl;
  }
  double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
  std::cout << "Average fitness: " << avg_fit << std::endl;
  double fit_var;
  for (int i = 0; i < fitnesses.size(); ++i){
    fit_var += pow((fitnesses[i]-avg_fit),2);
  }
  fit_var /= fitnesses.size();
  double fit_stddev = sqrt(fit_var);
  std::cout << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

  double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
  double time_var;
  for (int i = 0; i < times.size(); ++i){
    time_var += pow((times[i]-avg_time),2);
  }
  time_var /= times.size();
  double time_stddev = sqrt(time_var);
  std::cout << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
//  std::vector<uint_fast32_t> path(11);
//  std::iota(path.begin(), path.end(), 0);
//
//  for (size_t i = 0; i < 11; i++) {
//    std::cout << path[i] << " ";
//  }
//  std::cout << std::endl;
//
//  std::random_device rd;
//  std::mt19937 g(rd());
//  std::shuffle(path.begin() + 1, path.end() - 1, g);
//
//  for (size_t i = 0; i < 11; i++) {
//    std::cout << path[i] << " ";
//  }
//  std::cout << std::endl;
//
//  std::pair<std::vector<uint_fast32_t>, double> ret = two_opt(path, cost_mat);
//
//  path = ret.first;
//  for (size_t i = 0; i < 11; i++) {
//    std::cout << path[i] << " ";
//  }
//  std::cout << std::endl;
//
//  std::cout << ret.second << std::endl;
//
//  for (int i = 0; i < 20; i++) {
//    std::cout << g() << " ";
//  }



  return 0;
}
