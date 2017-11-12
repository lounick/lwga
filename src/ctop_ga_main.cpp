//
// Created by nick on 31/03/16.
//

#include "ctop_ga.h"

int main(){
  std::vector<std::pair<double, double> > nodes;
  std::vector<std::vector<double> > cost_mat;
  /*
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
  uint_fast32_t num_robots = 2;
  std::vector<double> max_cost_v(num_robots, 0.75+51/4);
   */
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
  uint_fast32_t num_robots = 3;
//  std::vector<double> max_cost_v(num_robots, 4*((num_robots-1)+(82+81.0)/num_robots)/4);
  std::vector<double> max_cost_v(num_robots, 2*((82+81.0)/num_robots)/4.0);
  std::vector<double> fitnesses;
  std::vector<double> times;

  int nexp = 100;

  std::random_device rd;
  std::mt19937 g(rd());

  Chromosome best(cost_mat.size(), num_robots, 0, 82);
  double best_fit = 0;

  for (int exp = 0; exp < nexp; exp++) {
    auto start = std::chrono::high_resolution_clock::now();
//    Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 26, g);
    Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 82, g, 200, 100, 0.0, 1.0, "RANDOM", 0.05);
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
              extras += std::exp((log(0.01) / 2) * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
            }
          }
          fitness += rewards[c.genes[robot].path[i]] + extras;
        }
      }
    }
    c.evaluate_chromosome(cost_mat, rewards, max_cost_v);
    std::cout << exp << " " << c.total_fitness << " " << fitness << std::endl;
    fitnesses.push_back(fitness);
    if(fitness > best_fit){
          best = c;
          best_fit = fitness;
    };
  }
  double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
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

  std::cout << "best_fit: " << best_fit << " max_cost: " << max_cost_v[0] << std::endl;
  for(size_t i = 0; i < num_robots; i++){std::cout << best.genes[i].cost << " "<< get_path_cost(best.genes[i].path, cost_mat) << std::endl;}
  std::cout << "[";
  for(size_t i = 0; i < num_robots; i++){
//    std::cout << best.genes[i].cost << " [";
    std::cout << " [";
    for(size_t j = 0; j < best.genes[i].path.size(); j++){
      if(j < best.genes[i].path.size() - 1)
        std::cout << best.genes[i].path[j] << ", ";
      else {
        if (i < num_robots - 1) {
          std::cout << best.genes[i].path[j] << "], ";
        } else {
          std::cout << best.genes[i].path[j] << "]";
        }
      }
    }
  }
  std::cout << "]" << std::endl;
  return 0;
}