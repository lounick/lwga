//
// Created by nick on 14/01/18.
//

#include "ga_utils.h"
#include "dcop_ga.h"
#include <iostream>
int main(int argc, char *argv[]){
  std::shared_ptr<Vector<Point2D>> nodes = std::make_shared<Vector<Point2D>>();
  nodes->push_back(std::make_pair(0, 0));
  for (int i = 1; i < 6; i++) {
    for (int j = -2; j < 3; j++) {
      nodes->push_back(std::make_pair(i, j));
    }
  }
  nodes->push_back(std::make_pair(6, 0));
  Vector<double_t> std_angles;
  std_angles.reserve(36);
  for(uint_fast8_t i = 0; i < 36; ++i){
    std_angles.push_back(i*10*M_PI/180);
  }
//  Path path_before {0,3,7,12,13,8,2,6,11,17,18,23,26};
//  Vector<double_t> angles_before = {
//      0.0, -M_PI_4, -M_PI_4, M_PI_4, M_PI_2+M_PI_4, -M_PI_2-M_PI_4, -M_PI_4, -M_PI_4, M_PI_4, M_PI_2, M_PI_4,0.0,0.0
//  };
//  Path path_after;
//  Vector<double_t> angles_after;
//  double_t cost_before, cost_after;
//  double_t rho = 2;
//  cost_before = get_dubins_path_cost(nodes, rho, path_before, angles_before);
//  std::tie(path_after, angles_after, cost_after) = dubins_two_opt(nodes, rho
//      , path_before, angles_before, cost_before);
//  std::cout << cost_before << std::endl;
//  std::cout << cost_after << std::endl;
//  std::string comma = "";
//  std::cout << "[";
//  for (size_t i = 0; i < path_after.size(); ++i){
//    std::cout << comma << path_after[i];
//    comma = ", ";
//  }
//  std::cout << "]" << std::endl;
//  comma = "";
//  std::cout << "[";
//  for (size_t i = 0; i < angles_after.size(); ++i){
//    std::cout << comma << angles_after[i];
//    comma = ", ";
//  }
//  std::cout << "]" << std::endl;
  Matrix <double_t> cost_mat;
  for (int i = 0; i < nodes->size(); i++) {
    std::vector<double_t > tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < nodes->size(); i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < nodes->size(); j++) {
      double dist = find_distance(nodes->at(i), nodes->at(j));
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  std::vector<double> rewards(nodes->size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;
  std::vector<double> fitnesses;
  std::vector<double> times;

  int nexp = 1;

  std::random_device rd;
  std::mt19937 g(rd());

  double rho = 0.5;

  for (int exp = 0; exp < nexp; exp++) {
    auto start = std::chrono::high_resolution_clock::now();
    dcop_ga::Chromosome c = dcop_ga::ga_dcop(
        nodes, std_angles, rho, cost_mat, rewards, 30, 0, 26, g);
//    dcop_ga::Chromosome c = dcop_ga::generate_chromosome(
//        nodes, std_angles, rho, 30,
//        0, 26, cost_mat, g);
    std::tie(c.path, c.angles, c.cost) =
        dubins_two_opt(nodes, rho, c.path, c.angles, c.cost);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());
//    for (uint_fast32_t vertex:c.path) {
//      std::cout << vertex << " ";
//    }
//    std::cout << std::endl;
    Path path_after;
    Vector<double_t> angles_after;
    double_t cost_before, cost_after;
    cost_before = get_dubins_path_cost(nodes, rho, c.path, c.angles);
    std::cout << cost_before << std::endl;
    print_path(c.path);
    print_vector<double_t>(c.angles);
    for(int i =0; i < 100; ++i)
    std::tie(path_after, angles_after, cost_after) =
        dubins_two_opt(nodes, rho, c.path, c.angles, cost_before);
    std::cout << cost_after << std::endl;
    print_path(path_after);
    print_vector<double_t>(angles_after);

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
  }
  double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
  std::cout << "Average fitness: " << avg_fit << std::endl;
  double fit_var = 0.0;
  for (int i = 0; i < fitnesses.size(); ++i){
    fit_var += pow((fitnesses[i]-avg_fit),2);
  }
  fit_var /= fitnesses.size();
  double fit_stddev = sqrt(fit_var);
  std::cout << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

  double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
  double time_var = 0.0;
  for (int i = 0; i < times.size(); ++i){
    time_var += pow((times[i]-avg_time),2);
  }
  time_var /= times.size();
  double time_stddev = sqrt(time_var);
  std::cout << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
  return 0;
}