//
// Created by nick on 15/04/18.
//

#include <include/top_ga.h>
#include "dcop_ga.h"

int main(int argc, char *argv[]) {

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
  std::vector<double_t> replanned_fitnesses;
  std::vector<double> times;

  double_t total_planning_time;
  double_t max_cost = 82.01;

  int nexp = 50;
//35, 475, 3, 0.6, 0.8, 0.07
  std::random_device rd;
  std::mt19937 g(rd());
  dcop_ga::Chromosome best;
  double_t best_fitness = std::numeric_limits<double_t>::infinity();
  dcop_ga::Chromosome best_replanned;
  dcop_ga::Chromosome worst;
  double_t worst_fitness = std::numeric_limits<double_t>::max();
  for (int exp = 0; exp < nexp; exp++) {
    std::cout << "Experiment: " << exp << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    dcop_ga::Chromosome c = dcop_ga::ga_dcop(
        nodes, std_angles, rho, dubins_cost_mat,
        eucledian_cost_mat, rewards, 82.01, 0, 82, 475, 35, 3, 0.6, 0.8, 0.07, g);

    // TODO: Here insert the code for finding the remaining nodes and re planning.

    // Create and sort a vector of the visited indices.
//    Path visited_nodes;
//    visited_nodes.reserve(c.path.size());
//    Path final_path;
//    final_path.reserve(c.path.size());
//    Vector<uint_fast32_t> final_angles;
//    final_angles.reserve(c.path.size());
//    Point2D new_start;

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
    double_t cost = 0.0;
    dcop_ga::Chromosome new_c = c;
    new_c.path.clear();
    new_c.angles.clear();
    new_c.seen_vertices.clear();
    new_c.free_vertices = c.all_vertices;
    new_c.fitness = 0.0;
    new_c.cost = 0.0;

    for (size_t i = 0; i < c.path.size()-1; i++) {
      double extras = 0;
      insert_ret = seen.insert(c.path[i]);
      if (insert_ret.second) {
        for (size_t j = 0; j < free_vertices.size(); j++) {
          if (eucledian_cost_mat[c.path[i]][free_vertices[j]] < 2 && rewards[c.path[i]] > 0) {
            extras += std::exp((log(0.01)/2) * eucledian_cost_mat[c.path[i]][free_vertices[j]]);
          }
        }
        fitness += rewards[c.path[i]] + extras;
      }
      cost += dubins_cost_mat[c.path[i]][c.path[i+1]][c.angles[i]][c.angles[i+1]];
      if (cost < max_cost/10.0){
//        visited_nodes.push_back(c.path[i]);
//        final_path.push_back(c.path[i]);
//        final_angles.push_back(c.angles[i]);
//        new_start = nodes->at(c.path[i]);
        new_c.path.push_back(c.path[i]);
        new_c.angles.push_back(c.angles[i]);
        new_c.seen_vertices.insert(c.path[i]);
        new_c.free_vertices.erase(std::remove(new_c.free_vertices.begin(),
                                              new_c.free_vertices.end(),
                                              c.path[i]));
      }
      if (cost > max_cost/2.0) {
        break;
      }
    }
    new_c.free_vertices.erase(std::remove(new_c.free_vertices.begin(),
                                          new_c.free_vertices.end(),
                                          c.path.back()));
    new_c.evaluate_chromosome(dubins_cost_mat, eucledian_cost_mat, rewards);
    fitnesses.push_back(fitness);



//    std::sort(visited_nodes.begin(), visited_nodes.end(), std::greater<uint_fast32_t >());

    // Create new nodes vector
//    std::shared_ptr<Vector<Point2D>> new_nodes = std::make_shared<Vector<Point2D>>(nodes->begin(), nodes->end());

    // Remove visited nodes
//    for (uint_fast32_t node : visited_nodes) {
//      new_nodes->erase(new_nodes->begin() + node);
//    }

    // Find the starting idx
//    size_t new_start_idx = 0;
//    for (size_t i = 0; i < new_nodes->size(); ++i){
//      if (new_nodes->at(i) == new_start) {
//        new_start_idx = i;
//        break;
//      }
//    }

    // Calculate new costs etc
//    Matrix <double_t> new_eucledian_cost_mat;
//    for (size_t i = 0; i < new_nodes->size(); i++) {
//      std::vector<double_t > tmp_vec;
//      tmp_vec.reserve(new_nodes->size());
//      new_eucledian_cost_mat.push_back(tmp_vec);
//    }
//
//    for (size_t i = 0; i < new_nodes->size(); i++) {
//      new_eucledian_cost_mat[i].push_back(0);
//      for (size_t j = i + 1; j < new_nodes->size(); j++) {
//        double dist = find_distance(new_nodes->at(i), new_nodes->at(j));
//        new_eucledian_cost_mat[i].push_back(dist);
//        new_eucledian_cost_mat[j].push_back(dist);
//      }
//    }
//
//    Matrix<Matrix<double_t>> new_dubins_cost_mat;
//    new_dubins_cost_mat.reserve(new_nodes->size());
//    for (size_t i = 0; i < new_nodes->size(); ++i) {
//      Vector<Matrix<double_t>> tmp_vec;
//      tmp_vec.reserve(new_nodes->size());
//      for(size_t j = 0; j < new_nodes->size(); ++j){
//        Matrix<double_t> angle_mat;
//        angle_mat.reserve(std_angles.size());
//        for (size_t k = 0; k < std_angles.size(); ++k){
//          Vector<double_t> angle_vec;
//          angle_vec.reserve(std_angles.size());
//          angle_mat.push_back(angle_vec);
//        }
//        tmp_vec.push_back(angle_mat);
//      }
//      new_dubins_cost_mat.push_back(tmp_vec);
//    }

//    for (size_t i = 0; i < new_nodes->size(); i++) {
//      double_t q0[3] = {new_nodes->at(i).first, new_nodes->at(i).second, 0};
//      std::unique_ptr<DubinsPath> ptp_path = std::make_unique<DubinsPath>();
//      for (size_t j = 0; j < new_nodes->size(); j++) {
//        double_t q1[3] = {new_nodes->at(j).first, new_nodes->at(j).second, 0};
//        for (size_t k = 0; k < std_angles.size(); ++k){
//          q0[2] = std_angles[k];
//          for (size_t l = 0; l < std_angles.size(); ++l){
//            q1[2] = std_angles[l];
//            int ret = dubins_init(q0, q1, rho, ptp_path.get());
//            if (ret != 0)
//              std::cout << "Dubins ret: " << ret << std::endl;
//            new_dubins_cost_mat[i][j][k].push_back(dubins_path_length(ptp_path.get()));
//          }
//        }
//      }
//    }

//    std::vector<double> new_rewards(new_nodes->size(), 1);
//    new_rewards.front() = 0;
//    new_rewards.back() = 0;
//
//    double_t new_max_cost = 82.01/4.0;
//    std::cout << new_max_cost << std::endl;
    start = std::chrono::high_resolution_clock::now();
    double_t replanned_fit = -std::numeric_limits<double_t>::infinity();
    dcop_ga::Chromosome br;
    for (size_t cnt = 0; cnt < 3; ++cnt) {
      dcop_ga::Chromosome init = new_c;
      dcop_ga::Chromosome replanned_c = dcop_ga::replan_ga_dcop(
          nodes, std_angles, rho, dubins_cost_mat,
          eucledian_cost_mat, rewards,
          max_cost / 2.0, 0, nodes->size() - 1,
          475, 35, 3, 0.6, 0.8, 0.07, g, init);
      if (replanned_c.fitness > replanned_fit) {
        br = replanned_c;
        replanned_fit = replanned_c.fitness;
      }
    }
    new_c = br;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());

//    std::cout << new_c.cost << std::endl;
//    std::cout << new_c.fitness << std::endl;
//    new_c.calculate_cost(dubins_cost_mat);
//    new_c.evaluate_chromosome(dubins_cost_mat, eucledian_cost_mat, rewards);
//    std::cout << new_c.cost << std::endl;
//    std::cout << new_c.fitness << std::endl;

    // TODO: Calculate the fitness of the replanned path

    // For each node of the path find the according index of the previous node and add it.

//    for (size_t i = 0; i < new_c.path.size(); ++i) {
//      Point2D node = new_nodes->at(new_c.path[i]);
//      for (size_t j = 0; j < nodes->size(); ++j) {
//        if (node == nodes->at(j)) {
//          final_path.push_back(j);
//          break;
//        }
//      }
//      final_angles.push_back(new_c.angles[i]);
//    }
//
    double new_fitness = 0;
    std::iota(vertices.begin(), vertices.end(), 0);
    free_vertices.clear();
    visited_vertices = c.path;
    std::sort(visited_vertices.begin(), visited_vertices.end());
    std::set_difference(vertices.begin(),
                        vertices.end(),
                        visited_vertices.begin(),
                        visited_vertices.end(),
                        std::back_inserter(free_vertices));
    seen.clear();
//    for (size_t i = 0; i < final_path.size(); i++) {
//      double extras = 0;
//      insert_ret = seen.insert(final_path[i]);
//      if (insert_ret.second) {
//        for (size_t j = 0; j < free_vertices.size(); j++) {
//          if (eucledian_cost_mat[final_path[i]][free_vertices[j]] < 2 && rewards[final_path[i]] > 0) {
//            extras += std::exp((log(0.01)/2) * eucledian_cost_mat[final_path[i]][free_vertices[j]]);
//          }
//        }
//        new_fitness += rewards[final_path[i]] + extras;
//      }
//    }
    for (size_t i = 0; i < new_c.path.size()-1; i++) {
      double extras = 0;
      insert_ret = seen.insert(new_c.path[i]);
      if (insert_ret.second) {
        for (size_t j = 0; j < free_vertices.size(); j++) {
          if (eucledian_cost_mat[new_c.path[i]][free_vertices[j]] < 2
              && rewards[new_c.path[i]] > 0) {
            extras += std::exp((log(0.01) / 2)
                                   * eucledian_cost_mat[new_c.path[i]][free_vertices[j]]);
          }
        }
        new_fitness += rewards[new_c.path[i]] + extras;
      }
    }
    if (new_fitness < fitness) {
      std::cout << fitness << " " << new_fitness << std::endl;
    }
    std::flush(std::cout);
    replanned_fitnesses.push_back(new_fitness);
    if (fitness < best_fitness) {
      best_fitness = fitness;
      best = c;
      best_replanned = new_c;
    }
  }

  double_t cost = 0.0;
  Path best_path;
  Vector<double_t> angles;
  std::cout << "best_path = ";
  print_vector(best.path);
  for (size_t i = 0; i < best.path.size(); ++i) {
    angles.push_back(std_angles[best.angles[i]]);
  }
  std::cout << "best_angles = ";
  print_vector(angles);

  angles.clear();
  for (size_t i = 0; i < best.path.size()-1; ++i) {
    cost += dubins_cost_mat[best.path[i]][best.path[i+1]][best.angles[i]][best.angles[i+1]];
    if (cost < 82.01/4.0 + 82.01/4.0) {
      best_path.push_back(best.path[i]);
      angles.push_back(std_angles[best.angles[i]]);
    }
  }
  best_path.push_back(best.path.back());
  angles.push_back(0.0);
  std::cout << "no_replan_path = ";
  print_vector(best_path);
  std::cout << "no_replan_angles = ";
  print_vector(angles);

  angles.clear();
  for (size_t i = 0; i < best_replanned.path.size(); ++i) {
    angles.push_back(std_angles[best_replanned.angles[i]]);
  }
  std::cout << "replan_path = ";
  print_vector(best_replanned.path);
  std::cout << "replan_angles = ";
  print_vector(angles);
//  std::cout << "Best:" << std::endl;
//  print_path(best.path);
//  Vector<double_t> angles;
//  for (size_t i = 0; i < best.angles.size(); ++i) {
//    angles.push_back(std_angles[best.angles[i]]);
//  }
//  print_vector<double_t>(angles);
//  std::cout << best_fitness << std::endl;
//  std::cout << "Worst:" << std::endl;
//  print_path(worst.path);
//  angles.clear();
//  for (size_t i = 0; i < worst.angles.size(); ++i) {
//    angles.push_back(std_angles[worst.angles[i]]);
//  }
//  print_vector<double_t>(angles);
//  std::cout << worst_fitness << std::endl;
  double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
//  std::cout << "Average fitness: " << avg_fit << std::endl;
  double fit_var = 0.0;
  for (int i = 0; i < fitnesses.size(); ++i){
    fit_var += pow((fitnesses[i]-avg_fit),2);
  }
  fit_var /= fitnesses.size();
  double fit_stddev = sqrt(fit_var);
  std::cout << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

  double_t avg_replanned = std::accumulate(replanned_fitnesses.begin(), replanned_fitnesses.end(), 0.0)/replanned_fitnesses.size();
  double_t replaned_var = 0.0;
  for (int i = 0; i < replanned_fitnesses.size(); ++i) {
    replaned_var += pow(replanned_fitnesses[i]-avg_replanned,2);
  }
  replaned_var /= replanned_fitnesses.size();
  double_t replanned_stdev = sqrt(replaned_var);
  std::cout << "Average replanned fitness: " << avg_replanned << " Variance: " << replaned_var << " StdDev: " << replanned_stdev << std::endl;

  double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
  double time_var = 0.0;
  for (int i = 0; i < times.size(); ++i){
    time_var += pow((times[i]-avg_time),2);
  }
  time_var /= times.size();
  double time_stddev = sqrt(time_var);
  std::cout << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;

  std::cout << "fitnesses = { 0.0: ";
  print_vector(fitnesses);
  std::cout << "}" << std::endl;
  std::cout << "replanned_fitnesses = { 0.0: ";
  print_vector(replanned_fitnesses);
  std::cout << "}" << std::endl;

  for (size_t i = 0; i < fitnesses.size(); ++i) {
    replanned_fitnesses[i] = replanned_fitnesses[i]/fitnesses[i] - 1;
  }
  avg_replanned = std::accumulate(replanned_fitnesses.begin(), replanned_fitnesses.end(), 0.0)/replanned_fitnesses.size();
  replaned_var = 0.0;
  for (int i = 0; i < replanned_fitnesses.size(); ++i) {
    replaned_var += pow(replanned_fitnesses[i]-avg_replanned,2);
  }
  replaned_var /= replanned_fitnesses.size();
  replanned_stdev = sqrt(replaned_var);
  std::cout << "Average replanned fitness: " << avg_replanned << " Variance: " << replaned_var << " StdDev: " << replanned_stdev << std::endl;
  return 0;
}