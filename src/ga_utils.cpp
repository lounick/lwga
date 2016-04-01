//
// Created by nick on 31/03/16.
//

#include "ga_utils.h"
#include <cmath>

double find_distance(Point2D a, Point2D b) {
  return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
}