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

typedef std::pair<double, double> Vertex;

double find_distance(Point2D a, Point2D b);

std::vector<size_t> get_population_sample(size_t pop_size, int samples, std::mt19937 &g);

Path two_opt_swap(Path &path, size_t &i, size_t &k);

std::pair<Path, Vector<uint_fast32_t >> two_opt_swap(
    const Path &path, const Vector<uint_fast32_t > &angles, size_t &i, size_t &k);

std::pair<Path, double> two_opt(Path &path, const Matrix<double> &cost_mat);

std::tuple<Path, Vector<uint_fast32_t >, double> dubins_two_opt(
    Matrix<Matrix<double_t>>&dubins_cost_mat,
    Vector<double_t> &std_angles,
    const std::shared_ptr<const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<uint_fast32_t > &angles, double_t cost);

void mutual_two_opt(Path &path1, Path &path2, const Matrix<double_t> &cost_mat, double_t max_cost1, double_t max_cost2);

double get_path_cost(Path &path, const Matrix<double> &cost_mat);

double_t get_dubins_path_cost(
    const std::shared_ptr< const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<double_t> &angles);

double calculate_fitness(Path p, std::vector<std::vector<double> > cost_mat, std::vector<double> rewards);

void print_path(Path p);

std::pair<std::vector<Vertex>, std::vector<double>> generate_grid (double x_size, double y_size, std::pair<double, double> idx_start);

bool logically_equal(double a, double b, double error_factor = 1.0);

template <typename T>
void print_vector(Vector<T> v){
  std::cout << "[";
  std::string comma = "";
  for (typename Vector<T>::iterator it = v.begin(); it != v.end(); ++it){
    std::cout << comma << std::setprecision(10) << *it;
    comma = ", ";
  }
  std::cout << "]" << std::endl;
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


double_t normalise_angle(double_t angle);

std::pair<Vector<uint_fast32_t>, double_t> straighten_path(
    Matrix<Matrix<double_t>>&dubins_cost_mat,
    const std::shared_ptr< const std::vector<Point2D>> nodes, double_t rho,
    Path &path, Vector<uint_fast32_t> &angles, double_t cost);

size_t bin_angle(double_t angle, double_t bin_size);
#endif //LWGA_GA_UTILS_H
