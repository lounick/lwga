//
// Created by nick on 16/03/18.
//

#include <mdmtsp_solver.h>

int main(int argc, char *argv[]) {

  uint_fast32_t num_vehicles = 2;
  uint_fast32_t num_depots = 3;
  Vector<Point2D> nodes =
      {{-1.5, 0}, {-1.5, 1}, {-0.5, 2}, {-0.5, 3}, {0.5, 3}, {0.5, 2}, {1.5, 0},
       {1.5, 1}};
  Vector<Point2D> depots = {{0, 4}, {-1.5, 0}, {1.5, 0}};
  depots = {{14, 0}, {0, 0}, {0, 0}};
  nodes = {{3, 3}, {5, -2}, {5, 4}, {6, 6}, {8, 3}, {9, -3}, {9, 5}, {10, 8},
           {11, 3}, {13, 2}, {8,6}, {7,5}, {6,9}};
  Vector<Point2D> problem = depots;
  problem.insert(std::end(problem), std::begin(nodes), std::end(nodes));

  Matrix<double_t> cost_mat;
  for (int i = 0; i < problem.size(); i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < problem.size(); ++i) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < problem.size(); ++j) {
      double dist = find_distance(problem[i], problem[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

//  for (size_t i = num_depots; i < problem.size(); i += 2) {
//    cost_mat[i][i + 1] = 0.0;
//    cost_mat[i + 1][i] = 0.0;
//  }

  Matrix<uint_fast32_t> tours;
  bool success = runMDMTSP(
      num_vehicles, nodes.size(), num_depots, cost_mat, true, tours);

  if (success) {
    for (auto &tour:tours) {
      print_vector(tour, std::cout);
    }
  } else {
    std::cout << "Solution not found" << std::endl;
  }

  return 0;
}