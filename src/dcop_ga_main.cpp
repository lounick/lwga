//
// Created by nick on 14/01/18.
//

#include "ga_utils.h"
#include "dcop_ga.h"
#include <iostream>
int main(int argc, char *argv[]){
  std::shared_ptr<Vector<Point2D>> nodes = std::make_shared<Vector<Point2D>>();
  nodes->push_back(std::make_pair(0, 0));
  for (int i = 1; i < 10; i++) {
    for (int j = -4; j < 5; j++) {
      nodes->push_back(std::make_pair(i, j));
    }
  }
  nodes->push_back(std::make_pair(10, 0));
  Vector<double_t> std_angles;
  uint_fast32_t degrees = 15;
  std_angles.reserve(360/degrees);
  for(uint_fast8_t i = 0; i < 360/degrees; ++i){
    std_angles.push_back(i*degrees*M_PI/180.0);
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
  Matrix <double_t> eucledian_cost_mat;
  for (size_t i = 0; i < nodes->size(); i++) {
    std::vector<double_t > tmp_vec;
    tmp_vec.reserve(nodes->size());
    eucledian_cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < nodes->size(); i++) {
    eucledian_cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < nodes->size(); j++) {
      double dist = find_distance(nodes->at(i), nodes->at(j));
      eucledian_cost_mat[i].push_back(dist);
      eucledian_cost_mat[j].push_back(dist);
    }
  }

  double rho = 0.5;
  Matrix<Matrix<double_t>> dubins_cost_mat;
  dubins_cost_mat.reserve(nodes->size());
  for (size_t i = 0; i < nodes->size(); ++i) {
    Vector<Matrix<double_t>> tmp_vec;
    tmp_vec.reserve(nodes->size());
    for(size_t j = 0; j < nodes->size(); ++j){
      Matrix<double_t> angle_mat;
      angle_mat.reserve(std_angles.size());
      for (size_t k = 0; k < std_angles.size(); ++k){
        Vector<double_t> angle_vec;
        angle_vec.reserve(std_angles.size());
        angle_mat.push_back(angle_vec);
      }
      tmp_vec.push_back(angle_mat);
    }
    dubins_cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < nodes->size(); i++) {
    double_t q0[3] = {nodes->at(i).first, nodes->at(i).second, 0};
    std::unique_ptr<DubinsPath> ptp_path = std::make_unique<DubinsPath>();
    for (size_t j = 0; j < nodes->size(); j++) {
      double_t q1[3] = {nodes->at(j).first, nodes->at(j).second, 0};
      for (size_t k = 0; k < std_angles.size(); ++k){
        q0[2] = std_angles[k];
        for (size_t l = 0; l < std_angles.size(); ++l){
          q1[2] = std_angles[l];
          int ret = dubins_init(q0, q1, rho, ptp_path.get());
          if (ret != 0)
            std::cout << "Dubins ret: " << ret << std::endl;
          dubins_cost_mat[i][j][k].push_back(dubins_path_length(ptp_path.get()));
        }
      }
    }
  }

  std::vector<double> rewards(nodes->size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;
  std::vector<double> fitnesses;
  std::vector<double> times;

  int nexp = 10;

  std::random_device rd;
  std::mt19937 g(rd());
  dcop_ga::Chromosome best;
  double_t best_fitness = std::numeric_limits<double_t>::min();
  for (int exp = 0; exp < nexp; exp++) {
    std::cout << "Experiment: " << exp << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    dcop_ga::Chromosome c = dcop_ga::ga_dcop(
        nodes, std_angles, rho, dubins_cost_mat,
        eucledian_cost_mat, rewards, 81.01, 0, 82, g);
//    dcop_ga::Chromosome c = dcop_ga::generate_chromosome(
//        nodes, std_angles, rho, 30,
//        0, 26, eucledian_cost_mat, g);
//    std::tie(c.path, c.angles, c.cost) =
//        dubins_two_opt(nodes, rho, c.path, c.angles, c.cost);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());
//    for (uint_fast32_t vertex:c.path) {
//      std::cout << vertex << " ";
//    }
//    std::cout << std::endl;
//    Path path_after = {0, 4, 9, 15, 20, 19, 13, 12, 17, 22, 21, 16, 11, 6, 1, 2, 3, 8, 14, 18, 24, 26};
//    Vector<double_t> angles_after = /*{0, 0.349066, 6.10865, 0.349066, 5.93412,
//    4.01426, 5.06145, 5.23599, 0.523599, 0.349066,
//    2.79253, 2.61799, 2.61799, 2.61799, 2.61799,
//    0.698132, 1.5708, 4.71239, 0.698132, 5.58505,
//    0.698132, 4.88692};*/
//        {0, 0.349066, 6.10865, 0.349066, 5.93412,
//        4.01426, 5.06145, 5.23599, 0.523599, 0.349066,
//        M_PI, M_PI, M_PI, M_PI, 2.61799,
//        0.698132, 1.5708, 4.71239, 0.698132, 5.58505,
//        0.698132, 4.88692};
//    double_t cost_before, cost_after;
//    cost_before = get_dubins_path_cost(nodes, rho, c.path, c.angles);
    std::cout << c.cost << std::endl;
    print_path(c.path);
    Vector<double_t> angles;
    for (size_t i = 0; i < c.angles.size(); ++i) {
      angles.push_back(std_angles[c.angles[i]]);
    }
    print_vector<double_t>(angles);
//    std::tie(path_after, angles_after, cost_after) =
//        dubins_two_opt(nodes, rho, path_after, angles_after, cost_before);
//        //dubins_two_opt(nodes, rho, c.path, c.angles, cost_before);
//    std::cout << cost_after << std::endl;
//    print_path(path_after);
//    print_vector<double_t>(angles_after);

    double fitness = 0;
    std::vector<uint_fast32_t> vertices(eucledian_cost_mat.size());
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
          if (eucledian_cost_mat[c.path[i]][free_vertices[j]] < 2) {
            extras += std::exp((log(0.01)/2) * eucledian_cost_mat[c.path[i]][free_vertices[j]]);
          }
        }
        fitness += rewards[c.path[i]] + extras;
      }
    }
    fitnesses.push_back(fitness);
    std::cout << fitness << std::endl;
    if (fitness > best_fitness) {
      best = c;
      best_fitness = fitness;
    }
  }
  std::cout << "Best:" << std::endl;
  print_path(best.path);
  Vector<double_t> angles;
  for (size_t i = 0; i < best.angles.size(); ++i) {
    angles.push_back(std_angles[best.angles[i]]);
  }
  print_vector<double_t>(angles);
  std::cout << best_fitness << std::endl;
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