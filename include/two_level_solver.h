//
// Created by nick on 17/03/18.
//

#ifndef LWGA_TWO_LEVEL_SOLVER_H
#define LWGA_TWO_LEVEL_SOLVER_H

#include <ga_types.h>
#include <ga_utils.h>
#include <cop_ga.h>
#include <mdmtsp_solver.h>

struct two_level_solution {
  double_t total_time;
  double_t total_utility;
  Matrix<uint_fast32_t> paths;
};

double_t getAreaStart(Point2D point, SamplingArea &area);

double_t getAreaStartEnd(SamplingArea &from, SamplingArea &to);

double_t getAreaEnd(Point2D point, SamplingArea &area);

two_level_solution solveTwoLevel(Point2D start, Point2D finish,
                                 Vector<SamplingArea> &problem,
                                 uint_fast32_t num_vehicles,
                                 const Vector<double_t> &budget,
                                 bool add_area_cost = false,
                                 bool fair = false);

#endif //LWGA_TWO_LEVEL_SOLVER_H
