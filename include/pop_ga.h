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
  double_t total_cost;
  double_t expected_cost;
  double_t reward;
  double_t expected_reward;
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
  bool grasp_estimated_reward;
  size_t population_size;
  double_t mutate_add_prob;
};

struct InsertMove {
  Path::const_iterator vertex_before;
  Path::const_iterator vertex_after;
  VertexId free_vertex;
  double_t score;
  double_t cost_increase;
  double_t heuristic_value;
};

using InsertMoveRet = std::pair<bool, InsertMove>;

struct CandidateList {
  Vector<InsertMove> insert_moves;
  double_t heuristic_min;
  double_t heuristic_max;
};

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
void RemoveVertex(Vector<VertexId> &v, VertexId vertex);

// TODO: Fill me in
// TODO: Docstring
Vector<VertexId> GetCommonVertices(const Path &p1, const Path &p2);

// TODO: Docstring
Chromosome GenerateCXChromosome(const Chromosome &p1, const Chromosome &p2,
                                VertexId split_vertex);

// TODO: Fill me in
// TODO: Docstring
Vector<Chromosome> InitialisePopulation(const Properties &properties,
                                        const Vector<double_t> &rewards,
                                        const Vector<double_t> &probs,
                                        const Matrix<double_t> &costs,
                                        rng::RandomNumberGenerator &rng);

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
CandidateList GenerateInsertMoves(const Path &p,
                                  const Vector<VertexId> &free_vertices,
                                  const Properties &properties,
                                  const Vector<double_t> &rewards,
                                  const Vector<double_t> &probs,
                                  const Matrix<double_t> &costs,
                                  const double_t current_cost);

// TODO: Fill me in
// TODO: Docstring
InsertMoveRet GenerateInsertMove(
    VertexId free_vertex, Path::const_iterator vertex_before,
    Path::const_iterator vertex_after, const Properties &properties,
    const Vector<double_t> &rewards, const Vector<double_t> &probs,
    const Matrix<double_t> &costs, const double_t current_cost);

// TODO: Fill me in
// TODO: Docstring
Chromosome TournamentSelect(const Vector<Chromosome> &pop,
                            size_t tournament_size,
                            rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
void SelectNewPopulation();

// TODO: Fill me in
// TODO: Docstring
void Crossover(Chromosome &p1, Chromosome &p2, const Properties &properties,
               const Vector<double_t> &rewards, const Vector<double_t> &probs,
               const Matrix<double_t> &costs, rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
void Mutate(Chromosome &c, const Properties &properties,
            const Vector<double_t> &rewards, const Vector<double_t> &probs,
            const Matrix<double_t> &costs, rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
void MutateRemove(Chromosome &c, const Properties &properties,
                  const Vector<double_t> &rewards,
                  const Vector<double_t> &probs, const Matrix<double_t> &costs,
                  rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
void MutateAdd(Chromosome &c, const Properties &properties,
               const Vector<double_t> &rewards, const Vector<double_t> &probs,
               const Matrix<double_t> &costs, rng::RandomNumberGenerator &rng);

// TODO: Fill me in
// TODO: Docstring
void EvaluateChromosome(Chromosome &c, const Properties &properties,
                        const Vector<double_t> &rewards,
                        const Vector<double_t> &probs,
                        const Matrix<double_t> &costs);

// TODO: Fill me in
// TODO: Docstring
Vector<double_t> GenerateActualServiceRequests();
}  // namespace pop_ga

#endif
