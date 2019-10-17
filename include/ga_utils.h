//
// Created by nick on 31/03/16.
//

#ifndef LWGA_GA_UTILS_H
#define LWGA_GA_UTILS_H

#include "ga_types.h"
#include <random>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <memory>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

double find_distance(Point2D a, Point2D b);

std::vector<size_t> get_population_sample(size_t pop_size,
                                          int samples,
                                          std::mt19937 &g);

Path two_opt_swap(Path &path, size_t &i, size_t &k);

std::pair<Path, Vector<uint_fast32_t >> two_opt_swap(
    const Path &path,
    const Vector<uint_fast32_t> &angles,
    size_t &i,
    size_t &k);

std::pair<Path, double> two_opt(Path &path, const Matrix<double> &cost_mat);

std::tuple<Path, Vector<uint_fast32_t>, double> dubins_two_opt(
    Matrix<Matrix<double_t>> &dubins_cost_mat,
    Vector<double_t> &std_angles,
    const std::shared_ptr<const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<uint_fast32_t> &angles, double_t cost);

void mutual_two_opt(Path &path1,
                    Path &path2,
                    const Matrix<double_t> &cost_mat,
                    double_t max_cost1,
                    double_t max_cost2);

double get_path_cost(Path &path, const Matrix<double> &cost_mat);

double_t get_dubins_path_cost(
    const std::shared_ptr<const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<double_t> &angles);

double calculate_fitness(Path p,
                         std::vector<std::vector<double> > cost_mat,
                         std::vector<double> rewards);

void print_path(Path p);

std::pair<std::vector<Point2D>, std::vector<double>> generate_grid(double x_size,
                                                                  double y_size,
                                                                  std::pair<
                                                                      double,
                                                                      double> idx_start);

bool logically_equal(double a, double b, double error_factor = 1.0);

bool less_equal(double_t a, double_t b, double_t error_factor = 1.0);

template<typename T>
void print_vector(Vector<T> v, std::ostream &out = std::cout) {
  out << "[";
  std::string comma = "";
  for (typename Vector<T>::iterator it = v.begin(); it != v.end(); ++it) {
    out << comma << std::setprecision(10) << *it;
    comma = ", ";
  }
  out << "]" << std::endl;
}

void print_vector(Vector<Point2D> v, std::ostream &out = std::cout);

template<typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

double_t normalise_angle(double_t angle);

std::pair<Vector<uint_fast32_t>, double_t> straighten_path(
    Matrix<Matrix<double_t>> &dubins_cost_mat,
    const std::shared_ptr<const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<uint_fast32_t> &angles, double_t cost);

size_t bin_angle(double_t angle, double_t bin_size);

Vector<Point2D> generate_sampling_grid(uint_fast32_t grid_width,
                                       uint_fast32_t grid_length,
                                       Point2D start = {0, 0},
                                       bool start_end = false);

Vector<SamplingArea> generate_random_problem(uint_fast32_t area_width,
                                             uint_fast32_t area_length,
                                             uint_fast32_t average_region_size,
                                             double_t coverage_percentage,
                                             uint8_t min_regions=0);

Vector<SamplingArea> generate_random_problem(uint_fast32_t area_width,
                                             uint_fast32_t area_length,
                                             uint_fast32_t num_regions,
                                             uint_fast32_t average_region_size,
                                             uint8_t min_regions = 0);

std::vector<std::string> split(const std::string &s, char delim);

bool LineIntersect(
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double x4, double y4,
    double &x, double &y);

#endif //LWGA_GA_UTILS_H
