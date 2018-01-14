//
// Created by nick on 31/03/16.
//

#include "ga_utils.h"
#include <cmath>
#include <iostream>
#include "dubins.h"

double find_distance(Point2D a, Point2D b) {
  return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
}

std::vector<size_t> get_population_sample(size_t pop_size, int samples, std::mt19937 &g) {

  std::vector<size_t> indices(pop_size);
  std::iota(indices.begin(), indices.end(), 0);
  size_t max = indices.size() - 1;

  std::vector<size_t> result;

  for (int i = 0; i < samples; ++i) {
    std::uniform_int_distribution<> d(0, max);
    size_t index = d(g);
    std::swap(indices[index], indices[max]);
    result.push_back(indices[max]);
    max--;
  }
  return result;
}

Path two_opt_swap(Path &path, size_t &i, size_t &k) {
  Path new_path = std::vector<uint_fast32_t>(path);
  std::reverse(new_path.begin() + i, new_path.begin() + k + 1);
  return new_path;
}

std::pair<Path, Vector<double_t>> two_opt_swap(
    const Path &path, const Vector<double_t> &angles, size_t &i, size_t &k) {
  Path new_path = std::vector<uint_fast32_t>(path);
  Vector<double_t> new_angles = Vector<double_t> (angles);
  std::reverse(new_path.begin() + i, new_path.begin() + k + 1);
  std::reverse(new_angles.begin() + i, new_angles.begin() + k + 1);
  return std::make_pair(new_path, new_angles);
};

std::pair<Path, double> two_opt(Path &path, const Matrix<double> &cost_mat) {

  bool start_again = true;
  Path tmp_path = std::vector<uint_fast32_t>(path);
  double tmp_path_cost = get_path_cost(tmp_path, cost_mat);
  double best_cost = tmp_path_cost;

  while (start_again) {
    start_again = false;

    for (size_t i = 1; i < tmp_path.size() - 2; ++i) {
      for (size_t k = i + 1; k < tmp_path.size() - 1; ++k) {
        Path new_path = two_opt_swap(tmp_path, i, k);

        //This works only for symmetric costs. for non-symmetric you must change that and calculate the whole reverse path cost.
        double new_cost = best_cost
            - cost_mat[tmp_path[i-1]][tmp_path[i]]
            - cost_mat[tmp_path[k]][tmp_path[k+1]]
            + cost_mat[tmp_path[i-1]][tmp_path[k]]
            + cost_mat[tmp_path[i]][tmp_path[k+1]];

        // Round it to avoid geting stuck in infinite looping due to machine rounding errors.
        new_cost = round( new_cost * 100000.0 ) / 100000.0;

        if (new_cost < best_cost) {
          tmp_path = new_path;
          start_again = true;
          best_cost = new_cost;
        }
      }
    }
  }

  return std::make_pair(tmp_path, best_cost);
}

double get_path_cost(Path &path, const Matrix<double> &cost_mat) {
  double cost = 0;
  for (Path::iterator it = path.begin() + 1; it != path.end(); ++it) {
    cost += cost_mat[*(it - 1)][*it];
    if (it != (path.end() - 1)) {
      cost += 1; //TODO: This assumes a fixed cost per task. Maybe that's not the case in the real world.
    }
  }
  return cost;
}
double calculate_fitness(Path p, std::vector<std::vector<double> > cost_mat, std::vector<double> rewards) {
  double fitness = 0;
  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices = p;
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));
  std::unordered_set<uint_fast32_t> seen;
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
  for (size_t i = 1; i < p.size() - 1; i++) {
    double extras = 0;
    insert_ret = seen.insert(p[i]);
    if (insert_ret.second) {
      for (size_t j = 0; j < free_vertices.size(); j++) {
        if (cost_mat[p[i]][free_vertices[j]] < 2) {
          extras += std::exp(log(0.01)/2 * cost_mat[p[i]][free_vertices[j]])*rewards[free_vertices[j]];
        }
      }
      fitness += rewards[p[i]] + extras;
    }
  }
  return fitness;
}

void print_path(Path p) {
  std::cout << "[";
  std::string comma = "";
  for (Path::iterator it = p.begin(); it != p.end(); ++it){
    std::cout << comma << *it;
    comma = ", ";
  }
  std::cout << "]" << std::endl;
}

std::pair<std::vector<Vertex>, std::vector<double>>
generate_grid(double x_size, double y_size, std::pair<double, double> idx_start) {

  return std::pair<std::vector<Vertex>, std::vector<double>>();
}
void
mutual_two_opt(Path &path1, Path &path2, const Matrix<double_t> &cost_mat, double_t max_cost1, double_t max_cost2) {
  size_t num_vertices = cost_mat.size();
  uint_fast32_t start_vertex = path1.front();
  uint_fast32_t end_vertex = path1.back();

  Path mutual = path1;
  Path::iterator it = path2.end() - 1;
  while (it != path2.begin()){
    mutual.push_back(*it);
    --it;
  }
  mutual.push_back(*it);

  std::pair<Path, double> two_opt_ret = two_opt(mutual,cost_mat);
  mutual = two_opt_ret.first;

  path1.clear();
  path2.clear();

  while (mutual.back() != end_vertex){
    path2.push_back(mutual.back());
    mutual.pop_back();
  }
  path2.push_back(mutual.back());
  mutual.pop_back();
  path1.assign(mutual.begin(), mutual.end());
}

double_t get_dubins_path_cost(
    const std::shared_ptr< const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<double_t> &angles) {
  double_t cost = 0;
  double_t q0[3];
  double_t q1[3];
  std::unique_ptr<DubinsPath> ptp_path = std::make_unique<DubinsPath>();
  for (size_t i = 1; i < path.size(); ++i){
    q0[0] = nodes->at(path[i-1]).first;
    q0[1] = nodes->at(path[i-1]).second;
    q0[2] = angles[i-1];
    q1[0] = nodes->at(path[i]).first;
    q1[1] = nodes->at(path[i]).second;
    q1[2] = angles[i];
    dubins_init(q0, q1, rho, ptp_path.get());
    cost += dubins_path_length(ptp_path.get());
  }
  return cost;
}

std::tuple<Path, Vector<double_t>, double> dubins_two_opt(
    const std::shared_ptr< const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<double_t> &angles, double_t cost) {
  bool start_again = true;
  Path tmp_path = Path(path);
  Vector<double_t> tmp_angles = Vector<double_t>(angles);
  double tmp_path_cost = round( cost * 100000.0 ) / 100000.0;
  double best_cost = tmp_path_cost;
  int count = 0;
  while (start_again) {
    start_again = false;

    for (size_t i = 1; i < tmp_path.size() - 2; ++i) {
      for (size_t k = i + 1; k < tmp_path.size() - 1; ++k) {
        Path new_path;
        Vector<double_t> new_angles;
        std::tie(new_path, new_angles) = two_opt_swap(tmp_path, angles, i, k);
        for (size_t idx = i; idx < k + 1; ++ idx) {
          new_angles[idx] -= M_PI;
        }

        //TODO: Maybe bin the angles to categories
        new_angles[i-1] = atan2(
            nodes->at(new_path[i]).second - nodes->at(new_path[i-1]).second,
            nodes->at(new_path[i]).first - nodes->at(new_path[i-1]).first);
        new_angles[i] = atan2(
            nodes->at(new_path[i+1]).second - nodes->at(new_path[i]).second,
            nodes->at(new_path[i+1]).first - nodes->at(new_path[i]).first);
        new_angles[k-1] = atan2(
            nodes->at(new_path[k]).second - nodes->at(new_path[k-1]).second,
            nodes->at(new_path[k]).first - nodes->at(new_path[k-1]).first);
        new_angles[k] = atan2(
            nodes->at(new_path[k+1]).second - nodes->at(new_path[k]).second,
            nodes->at(new_path[k+1]).first - nodes->at(new_path[k]).first);

        double new_cost = get_dubins_path_cost(nodes, rho, new_path, new_angles);
        // Round it to avoid geting stuck in infinite looping due to machine rounding errors.
        new_cost = round( new_cost * 100000.0 ) / 100000.0;

        if (new_cost < best_cost || logically_equal(new_cost, best_cost)) {
          tmp_path = new_path;
          tmp_angles = new_angles;
          start_again = true;
          best_cost = new_cost;
        }
      }
    }
    count++;
  }
  return std::make_tuple(tmp_path, tmp_angles, best_cost);
}

bool logically_equal(double a, double b, double error_factor) {
  return a == b ||
      std::abs(a - b) < std::abs(std::min(a, b)) * std::numeric_limits<double>::epsilon() *
          error_factor;
}