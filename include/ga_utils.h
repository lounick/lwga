//
// Created by nick on 31/03/16.
//

#ifndef LWGA_GA_UTILS_H
#define LWGA_GA_UTILS_H

#include "ga_types.h"
#include <random>
#include <algorithm>

double find_distance(Point2D a, Point2D b);

std::vector<size_t> get_population_sample(size_t pop_size, int samples, std::mt19937 &g);

Path two_opt_swap(Path &path, size_t &i, size_t &k);

std::pair<Path, double> two_opt(Path &path, Matrix<double> &cost_mat);

double get_path_cost(Path &path, Matrix<double> &cost_mat);

#endif //LWGA_GA_UTILS_H
