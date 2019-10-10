#ifndef LWGA_POP_GA_H
#define LWGA_POP_GA_H

#include <cfloat>
#include <cmath>
#include <utility>
#include "ga_types.h"
#include "ga_utils.h"
#include "rng.h"

namespace pop_ga {

enum class GenerationMethod {
  RANDOM,
  GRASP,
};

// TODO: Fill me in
// TODO: Docstring
struct Chromosome {
  Path p;
  Vector<VertexId> free_vertices;
  double_t cost;
  double_t reward;
  double_t fitness;
};

// TODO: Fill me in
// TODO: Docstring
struct Properties {
  GenerationMethod generation_method;
  VertexId start_id;
  VertexId end_id;
  double_t maximum_cost;
  double_t cost_per_time_unit;
  double_t grasp_greediness;
};

struct InsertMove {
  Path::const_iterator vertex_before;
  Path::const_iterator vertex_after;
  double_t score;
  double_t cost_increase;
  double_t heuristic_value;
};

using InsertMoveRet = std::pair<bool, InsertMove>;

// TODO: Docstring
double_t CalculateTotalReward(const Path &p, const Vector<double_t> &rewards);

// TODO: Docstring
double_t CalculateActualReward(const Path &p, const Vector<double_t> &rewards,
                               const Vector<double_t> &service);

// TODO: Docstring
double_t CalculateExpectedReward(const Path &p, const Vector<double_t> &rewards,
                                 const Vector<double_t> &probs);

// TODO: Docstring
double_t CalculateTotalCost(const Path &p, const Matrix<double_t> &costs);

// TODO: Docstring
double_t CalculateActualCost(const Path &p, const Matrix<double_t> &costs,
                             const Vector<double_t> &service);

// TODO: Docstring
double_t CalculateExpectedCost(const Path &p, const Matrix<double_t> &costs,
                               const Vector<double_t> &probs);

// TODO: Docstring
double_t CalculateTMaxObjective(const Path &p, const Vector<double_t> &rewards,
                                const Vector<double_t> &probs,
                                const Matrix<double_t> &costs,
                                double_t cost_per_time_unit);

// TODO: Docstring
double_t CalculateExpectedTimeObjective(const Path &p,
                                        const Vector<double_t> &rewards,
                                        const Vector<double_t> &probs,
                                        const Matrix<double_t> &costs,
                                        double_t cost_per_time_unit);

// TODO: Fill me in
// TODO: Docstring
Vector<Chromosome> InitialisePopulation();

// TODO: Fill me in
// TODO: Docstring
Chromosome GenerateChromosome(const Properties &properties,
                              const Vector<double_t> &rewards,
                              const Vector<double_t> &probs,
                              const Matrix<double_t> &costs,
                              rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
Chromosome GenerateRandomChromosome(const Properties &properties,
                                    const Matrix<double_t> &costs,
                                    rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
Chromosome GenerateGRASPChromosome(const Properties &properties,
                                   const Vector<double_t> &rewards,
                                   const Vector<double_t> &probs,
                                   const Matrix<double_t> &costs,
                                   rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
Vector<InsertMove> GenerateInsertMoves(const Path &p,
                                       const Vector<VertexId> &free_vertices,
                                       const Vector<double_t> &rewards,
                                       const Vector<double_t> &probs,
                                       const Matrix<double_t> &costs,
                                       const double_t max_cost);

// TODO: Fill me in
// TODO: Docstring
InsertMoveRet GenerateInsertMove(
    VertexId free_vertex, Path::const_iterator vertex_before,
    Path::const_iterator vertex_after, const Vector<double_t> &rewards,
    const Vector<double_t> &probs, const Matrix<double_t> &costs,
    const double_t current_cost, const double_t max_cost);

// TODO: Fill me in
// TODO: Docstring
size_t TournamentSelect();

// TODO: Fill me in
// TODO: Docstring
void SelectNewPopulation();

// TODO: Fill me in
// TODO: Docstring
void Crossover();

// TODO: Fill me in
// TODO: Docstring
void Mutate();

// TODO: Fill me in
// TODO: Docstring
void EvaluateChromosome();

// TODO: Fill me in
// TODO: Docstring
Vector<double_t> GenerateActualServiceRequests();
}  // namespace pop_ga

#endif
