//
// Created by nick on 18/03/18.
//

#include <two_level_solver.h>

int main(int argc, char *argv[]) {
//  Point2D start = {0, 0};
//  Point2D finish = {16, 0};
//  Vector<Point2D> grid1 = generate_sampling_grid(5, 5, {5, 5});
//  Vector<Point2D> grid2 = generate_sampling_grid(5, 5, {5, -5});
//
////  Vector<Point2D> problem;
////  problem.push_back(start);
////  problem.insert(std::end(problem), std::begin(grid1), std::end(grid1));
////  problem.insert(std::end(problem), std::begin(grid2), std::end(grid2));
////  problem.push_back(finish);
////  print_vector(problem);
//
//  SamplingArea s1({8, 5}, 5, 5);
//  s1.nodes(generate_sampling_grid(s1.width(),
//                                  s1.length(),
//                                  {s1.minX() - 1, s1.origin().second}));
////  s1.print(true);
//  SamplingArea s2({8, -5}, 5, 5);
//  s2.nodes(generate_sampling_grid(s2.width(),
//                                  s2.length(),
//                                  {s2.minX() - 1, s2.origin().second}));
////  s2.print(true);
//  Vector<SamplingArea> problem = {s1, s2};
//  two_level_solution s;
//  double_t budget = 14 + 24;

  Vector<SamplingArea>
      problem = generate_random_problem(100, 100, 9, 0.1, 0);
  Point2D start = {0, 0};
  Point2D finish = {100, 100};




  for (SamplingArea &a:problem) {
    Point2D area_start = {a.minX() - 1, a.origin().second};
    a.nodes(generate_sampling_grid(a.width(),a.length(), area_start));
  }

  two_level_solution s;
  double_t budget = 500;

  for (size_t i = 0; i < 10; ++i) {
    s = solveTwoLevel(start, finish, problem, 2, {budget, budget});
    std::cout << s.total_utility << " " << s.total_time << std::endl;
    std::cout << "[";
    print_vector(s.paths[0]);
    std::cout << ",";
    print_vector(s.paths[1]);
    std::cout << "]\n";
  }

}