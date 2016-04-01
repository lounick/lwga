//
// Created by nick on 31/03/16.
//

#ifndef LWGA_GA_TYPES_H
#define LWGA_GA_TYPES_H

#include <vector>
#include <cstdint>

template <typename T> using Matrix = std::vector<std::vector<T> >;
using Path = std::vector<uint_fast32_t>;
using Point2D = std::pair<double, double>;

#endif //LWGA_GA_TYPES_H
