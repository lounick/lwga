//
// Created by nick on 31/03/16.
//

#ifndef LWGA_GA_UTILS_H
#define LWGA_GA_UTILS_H

#include "ga_types.h"
#include <random>

double find_distance(Point2D a, Point2D b);

std::vector<size_t> get_population_sample(size_t pop_size, int samples, std::mt19937 &g);

#endif //LWGA_GA_UTILS_H
