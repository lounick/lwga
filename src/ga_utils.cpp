//
// Created by nick on 31/03/16.
//

#include "ga_utils.h"
#include <cmath>

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

std::pair<Path, double> two_opt(Path &path, Matrix<double> &cost_mat) {

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

double get_path_cost(Path &path, Matrix<double> &cost_mat) {
  double cost = 0;
  for (Path::iterator it = path.begin() + 1; it != path.end(); ++it) {
    cost += cost_mat[*(it - 1)][*it];
    if (it != (path.end() - 1)) {
      cost += 1; //TODO: This assumes a fixed cost per task. Maybe that's not the case in the real world.
    }
  }
  return cost;
}