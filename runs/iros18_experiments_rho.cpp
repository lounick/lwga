//
// Created by nick on 30/01/18.
//

#include <fstream>
#include <cstdlib>
#include <map>

#include "include/dcop_ga.h"

using namespace dcop_ga;

#include <sstream>
#include <iomanip>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 3)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

double_t get_fitness(
    const Chromosome &c, const Matrix<double_t> &cost_mat,
    const Vector<double_t> rewards) {
  std::unordered_set<uint_fast32_t> seen;
  for (size_t path_idx = 0; path_idx < c.path.size(); ++path_idx) {
    seen.insert(c.path[path_idx]);
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
  for (size_t i = 1; i < c.path.size() - 1; i++) {
    double extras = 0;
    insert_ret = seen.insert(c.path[i]);
    if (insert_ret.second) {
      for (size_t j = 0; j < free_vertices.size(); j++) {
        if (cost_mat[c.path[i]][free_vertices[j]] < 2) {
          extras += std::exp((log(0.01) / 2) * cost_mat[c.path[i]][free_vertices[j]]);
        }
      }
      fitness += rewards[c.path[i]] + extras;
    }
  }
}

int main(int argc, char *argv[]) {

  struct Configuration{
    uint16_t pop_size;
    uint8_t num_gen;
    uint8_t tour_size;
    double_t cx_rate;
    double_t mut_rate;
    double_t elite_rate;
  };

  Configuration conf = {450, 35, 3, 0.6, 0.8, 0.07};
  Vector<double_t> std_angles;
  uint_fast32_t degrees = 15;
  std_angles.reserve(360/degrees);
  for(uint_fast8_t i = 0; i < 360/degrees; ++i){
    std_angles.push_back(i*degrees*M_PI/180.0);
  }
  uint_fast32_t num_exp = 1000;
  Vector<double_t> rhos = {0.0000000000001, 0.25, 0.5, 0.75, 1.0};
  //Vector<double_t> cost_coef = {1.0, 0.75, 0.5, 0.25};
  Vector<double_t> cost_coef = {1.0};
  double_t max_cost;
  std::shared_ptr <Vector<Point2D>> nodes;
  Matrix<double_t> eucledian_cost_mat;
  Matrix<Matrix<double_t>> dubins_cost_mat;
  Vector<double_t> rewards;
  Vector<double_t> fitnesses;
  Vector<double_t> mean_fitnesses;
  Vector<double_t> stdev_fitnesses;
  Vector<double_t> times;
  Vector<double_t> mean_times;
  Vector<double_t> stdev_times;
  std::ofstream util_results_file;
  std::ofstream time_results_file;
  std::ofstream best_worst_file;
  std::random_device rd;
  std::mt19937 g(rd());

  /*
   * 5x5 Problem
   * 1.0 Generate nodes
   * 2.0 Calculate eucledian cost matrix
   * 3.0 Calculate rewards
   * 4.0 Assign max cost
   * 5.0 For each rho
   * * 5.1 Calculate dubins cost matrix
   * * 5.2 For each available cost
   * * * 5.2.1 For num_exp
   * * * * 5.2.1.1 Run ga_dcop
   * * * * 5.2.1.2 Store time and fitness
   * * * 5.2.2 Write in file time and fitness
   * 6.0 End
   */
  // Generate nodes
  int size = atoi(argv[1]);
  nodes = std::make_shared<Vector<Point2D >>();
  nodes->push_back(std::make_pair(0, 0));
  for (int i = 1; i < size + 1; i++) {
    for (int j = -(size - 1)/2; j < 1 + (size - 1)/2; j++) {
      nodes->push_back(std::make_pair(i, j));
    }
  }
  nodes->push_back(std::make_pair(size + 1, 0));

  // Calculate eucledian cost matrix
  eucledian_cost_mat.clear();
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

  // Calculate rewards
  rewards.clear();
  rewards = Vector<double_t>(nodes->size(),1);
  rewards.front() = 0.0;
  rewards.back() = 0.0;

  std::map<int,double_t> costs = {{5, 26.01}, {7, 50.9}, {9, 82.01}};
  // Assign max cost
  max_cost = costs[size];

  fitnesses.clear();
  mean_fitnesses.clear();
  stdev_fitnesses.clear();
  times.clear();
  mean_times.clear();
  stdev_times.clear();

  for (double_t rho:rhos) {
    // Calculate dubins cost matrix
    dubins_cost_mat.clear();
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
        if (q0[0] == q1[0] && q0[1] == q1[1]) {
          for (size_t k = 0; k < std_angles.size(); ++k){
            for (size_t l = 0; l < std_angles.size(); ++l){
              dubins_cost_mat[i][j][k].push_back(0);
            }
          }
        } else {
          for (size_t k = 0; k < std_angles.size(); ++k){
            q0[2] = std_angles[k];
            for (size_t l = 0; l < std_angles.size(); ++l){
              q1[2] = std_angles[l];
              int ret = dubins_init(q0, q1, rho, ptp_path.get());
              if (ret != 0) {
                std::cout << "Dubins ret: " << ret << std::endl;
                std::cout << "Dubins path lenght: "
                          << dubins_path_length(ptp_path.get()) << std::endl;
              }
              dubins_cost_mat[i][j][k].push_back(dubins_path_length(ptp_path.get()));
            }
          }
        }
      }
    }
    // For each available cost
    for (double_t cc : cost_coef) {
      double_t max_cost_exp = max_cost*cc;
      Chromosome best, worst;
      double_t best_fitness = std::numeric_limits<double_t>::min();
      double_t worst_fitness = std::numeric_limits<double_t>::max();
      fitnesses.clear();
      times.clear();
      for (uint_fast32_t exp = 0; exp < num_exp; ++exp) {
        auto start = std::chrono::high_resolution_clock::now();
        Chromosome c = ga_dcop(
            nodes, std_angles, rho, dubins_cost_mat,
            eucledian_cost_mat, rewards, max_cost_exp, 0, 1 + size*size,
            conf.pop_size, conf.num_gen, conf.tour_size, conf.cx_rate,
            conf.mut_rate, conf.elite_rate, g);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double_t> diff = end - start;
        times.push_back(diff.count());
        double_t fitness = std::cbrt(c.fitness*c.cost);//get_fitness(c, eucledian_cost_mat, rewards);
        fitnesses.push_back(fitness);
        if (fitness > best_fitness) {
          best_fitness = fitness;
          best = c;
        }
        if (fitness < worst_fitness) {
          worst_fitness = fitness;
          worst = c;
        }
        if (exp > 0 && exp%10 == 0)
          std::cout << "." << std::flush;
      }
      std::cout << std::endl;
      double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
      double fit_var;
      for (int i = 0; i < fitnesses.size(); ++i) {
        fit_var += pow((fitnesses[i] - avg_fit), 2);
      }
      fit_var /= fitnesses.size();
      double fit_stddev = sqrt(fit_var);
      mean_fitnesses.push_back(avg_fit);
      stdev_fitnesses.push_back(fit_stddev);
      double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
      double time_var;
      for (int i = 0; i < times.size(); ++i) {
        time_var += pow((times[i] - avg_time), 2);
      }
      time_var /= times.size();
      double time_stddev = sqrt(time_var);
      mean_times.push_back(avg_time);
      stdev_times.push_back(time_stddev);
      util_results_file.open(
          "/tmp/rho_" + std::to_string(size) + "x" + std::to_string(size) + "_"
              + to_string_with_precision(rho)+"_util.csv");
      std::vector<double_t>::iterator it;
      for (it = fitnesses.begin(); it != fitnesses.end() - 1; ++it){
        util_results_file << *it << ",";
      }
      util_results_file << fitnesses.back() << "\n";
      util_results_file.close();

      time_results_file.open(
          "/tmp/rho_" + std::to_string(size) + "x" + std::to_string(size) + "_"
              + to_string_with_precision(rho)+"_time.csv");
      for (it = times.begin(); it != times.end() - 1; ++it){
        time_results_file << *it << ",";
      }
      time_results_file << fitnesses.back() << "\n";
      time_results_file.close();
      best_worst_file.open(
          "/tmp/rho_" + std::to_string(size) + "x" + std::to_string(size) + "_"
              + to_string_with_precision(rho)+"_best_worst.txt");
      best_worst_file << "Best:\n";
      best_worst_file << "Fitness: "
                      << std::to_string(std::cbrt(best.fitness*best.cost))
                      << "\n";
      best_worst_file << "Path: ";
      print_vector(best.path, best_worst_file);
      best_worst_file << "Worst:\n";
      best_worst_file << "Fitness: "
                      << std::to_string(std::cbrt(worst.fitness*worst.cost))
                      << "\n";
      best_worst_file << "Path: ";
      print_vector(worst.path, best_worst_file);
      best_worst_file.close();
    }
  }
}