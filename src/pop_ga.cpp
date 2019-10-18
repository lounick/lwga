#include "pop_ga.h"
#include <algorithm>
#include <limits>
#include <set>

namespace pop_ga {
// TODO: Fill me in
// TODO: Docstring
// struct Chromosome {}

// TODO: Docstring
double_t CalculateTotalReward(const Path &p, const Vector<double_t> &rewards) {
  double_t reward = 0.0;
  for (const auto &vertex : p) {
    reward += rewards[vertex];
  }
  return reward;
}

// TODO: Docstring
double_t CalculateActualReward(const Path &p, const Vector<double_t> &rewards,
                               const Vector<double_t> &service) {
  double_t reward = 0.0;
  for (const auto &vertex : p) {
    reward += service[vertex] * rewards[vertex];
  }
  return reward;
}

// TODO: Docstring
double_t CalculateExpectedReward(const Path &p, const Vector<double_t> &rewards,
                                 const Vector<double_t> &probs) {
  double_t reward = 0.0;
  for (const auto &vertex : p) {
    reward += probs[vertex] * rewards[vertex];
  }
  return reward;
}

// TODO: Docstring
double_t CalculateTotalCost(const Path &p, const Matrix<double_t> &costs) {
  double_t cost = 0.0;
  Path::const_iterator it = p.begin();
  for (; it != p.end() - 1; ++it) {
    cost += costs[*it][*(it + 1)];
  }
  return cost;
}

// TODO: Docstring
double_t CalculateActualCost(const Path &p, const Matrix<double_t> &costs,
                             const Vector<double_t> &service) {
  return CalculateExpectedCost(p, costs, service);
}

// TODO: Docstring
double_t CalculateExpectedCost(const Path &p, const Matrix<double_t> &costs,
                               const Vector<double_t> &probs) {
  double_t cost = 0.0;
  Path::const_iterator it = p.begin();
  for (; it != p.end() - 1; ++it) {
    double_t sum = 0.0;
    Path::const_iterator sum_it = it + 1;
    for (; sum_it != p.end(); ++sum_it) {
      double_t product = 1.0;
      Path::const_iterator prod_it = it + 1;
      for (; prod_it < sum_it; ++prod_it) {
        product *= (1.0 - probs[*prod_it]);
      }
      sum += costs[*it][*sum_it] * probs[*sum_it] * product;
    }
    cost += probs[*it] * sum;
  }
  return cost;
}

// TODO: Docstring
double_t CalculateTMaxObjective(const Path &p, const Vector<double_t> &rewards,
                                const Vector<double_t> &probs,
                                const Matrix<double_t> &costs,
                                double_t cost_per_time_unit) {
  double_t expected_reward = CalculateExpectedReward(p, rewards, probs);
  double_t total_cost = CalculateTotalCost(p, costs);
  return expected_reward - cost_per_time_unit * total_cost;
}

// TODO: Docstring
double_t CalculateExpectedTimeObjective(const Path &p,
                                        const Vector<double_t> &rewards,
                                        const Vector<double_t> &probs,
                                        const Matrix<double_t> &costs,
                                        double_t cost_per_time_unit) {
  double_t expected_reward = CalculateExpectedReward(p, rewards, probs);
  double_t expected_cost = CalculateExpectedCost(p, costs, probs);
  return expected_reward - cost_per_time_unit * expected_cost;
}

// TODO: Docstring
void RemoveVertex(Vector<VertexId> &v, VertexId vertex) {
  v.erase(std::remove(v.begin(), v.end(), vertex), v.end());
}

// TODO: Docstring
Chromosome InitialiseEmptyChromosome(const Properties &properties,
                                     const Matrix<double_t> costs) {
  Chromosome c;
  c.total_cost = 0.0;
  c.reward = 0.0;
  c.fitness = 0.0;
  c.free_vertices = Vector<VertexId>(costs.size());
  std::iota(c.free_vertices.begin(), c.free_vertices.end(), 0);
  RemoveVertex(c.free_vertices, properties.start_id);
  RemoveVertex(c.free_vertices, properties.end_id);
  c.p.reserve(costs.size());
  return c;
}

// TODO: Docstring
Chromosome GenerateChromosome(const Properties &properties,
                              const Vector<double_t> &rewards,
                              const Vector<double_t> &probs,
                              const Matrix<double_t> &costs,
                              rng::RandomNumberGenerator &rng) {
  switch (properties.generation_method) {
    case GenerationMethod::RANDOM:
      return GenerateRandomChromosome(properties, costs, rng);
    case GenerationMethod::GRASP:
      return GenerateGRASPChromosome(properties, rewards, probs, costs, rng);
  }
}

// TODO: Docstring
Chromosome GenerateRandomChromosome(const Properties &properties,
                                    const Matrix<double_t> &costs,
                                    rng::RandomNumberGenerator &rng) {
  Chromosome c = InitialiseEmptyChromosome(properties, costs);
  // Insert start vertex
  // While not done
  //   Get random vertex
  //   If less than max_cost
  //     Insert
  //   Else
  //     done
  // Insert end vertex
  c.p.push_back(properties.start_id);
  bool done = false;
  while (!done) {
    size_t random_idx = rng.GenerateUniformInt(0, c.free_vertices.size() - 1);
    VertexId vertex = c.free_vertices[random_idx];
    double_t total_cost = c.total_cost + costs[c.p.back()][vertex] +
                          costs[vertex][properties.end_id];
    if (total_cost < properties.maximum_cost ||
        logically_equal(total_cost, properties.maximum_cost)) {
      c.total_cost += costs[c.p.back()][vertex];
      c.p.push_back(vertex);
      RemoveVertex(c.free_vertices, vertex);
    } else {
      done = true;
    }
    if (c.free_vertices.empty()) {
      done = true;
    }
  }
  c.total_cost += costs[c.p.back()][properties.end_id];
  c.p.push_back(properties.end_id);
  return c;
}

// TODO: Docstring
Chromosome GenerateGRASPChromosome(const Properties &properties,
                                   const Vector<double_t> &rewards,
                                   const Vector<double_t> &probs,
                                   const Matrix<double_t> &costs,
                                   rng::RandomNumberGenerator &rng) {
  Chromosome c = InitialiseEmptyChromosome(properties, costs);
  c.p.push_back(properties.start_id);
  c.p.push_back(properties.end_id);
  c.total_cost = costs[properties.start_id][properties.end_id];
  CandidateList cl = GenerateInsertMoves(c.p, c.free_vertices, properties,
                                         rewards, probs, costs, c.total_cost);
  while (!cl.insert_moves.empty()) {
    double_t threshold =
        cl.heuristic_min +
        properties.grasp_greediness * (cl.heuristic_max - cl.heuristic_min);
    Vector<InsertMove> rcl;
    rcl.reserve(cl.insert_moves.size());
    for (const InsertMove &im : cl.insert_moves) {
      if ((im.heuristic_value > threshold) ||
          logically_equal(im.heuristic_value, threshold)) {
        rcl.push_back(im);
      }
    }
    size_t random_im_idx = rng.GenerateUniformInt(0, rcl.size() - 1);
    InsertMove im = rcl[random_im_idx];
    c.p.insert(im.vertex_after, im.free_vertex);
    RemoveVertex(c.free_vertices, im.free_vertex);
    c.total_cost += im.cost_increase;
    cl.insert_moves.clear();
    cl = GenerateInsertMoves(c.p, c.free_vertices, properties, rewards, probs,
                             costs, c.total_cost);
  }
  // TODO: Perform 2-opt
  return c;
}

// TODO: Docstring
CandidateList GenerateInsertMoves(const Path &p,
                                  const Vector<VertexId> &free_vertices,
                                  const Properties &properties,
                                  const Vector<double_t> &rewards,
                                  const Vector<double_t> &probs,
                                  const Matrix<double_t> &costs,
                                  const double_t current_cost) {
  CandidateList cl;
  cl.insert_moves.reserve(free_vertices.size() * (p.size() - 1));
  cl.heuristic_min = std::numeric_limits<double_t>::infinity();
  cl.heuristic_max = -std::numeric_limits<double_t>::infinity();
  for (const VertexId free_vertex : free_vertices) {
    Path::const_iterator path_it = p.begin();
    for (; path_it < p.end() - 1; ++path_it) {
      InsertMoveRet im =
          GenerateInsertMove(free_vertex, path_it, path_it + 1, properties,
                             rewards, probs, costs, current_cost);
      if (im.first) {
        cl.insert_moves.push_back(im.second);
        if (im.second.heuristic_value > cl.heuristic_max) {
          cl.heuristic_max = im.second.heuristic_value;
        }
        if (im.second.heuristic_value < cl.heuristic_min) {
          cl.heuristic_min = im.second.heuristic_value;
        }
      }
    }
  }
  return cl;
}

// TODO: Docstring
InsertMoveRet GenerateInsertMove(
    VertexId free_vertex, Path::const_iterator vertex_before,
    Path::const_iterator vertex_after, const Properties &properties,
    const Vector<double_t> &rewards, const Vector<double_t> &probs,
    const Matrix<double_t> &costs, const double_t current_cost) {
  InsertMoveRet im;
  im.first = false;
  // Calcuclate cost increase
  im.second.cost_increase = costs[*vertex_before][free_vertex] +
                            costs[free_vertex][*vertex_after] -
                            costs[*vertex_before][*vertex_after];
  if ((current_cost + im.second.cost_increase < properties.maximum_cost) ||
      logically_equal(current_cost + im.second.cost_increase,
                      properties.maximum_cost)) {
    // Insert can be performed
    im.first = true;
    // TODO: Check different policies for score. I.e. if we take the
    // probabilities into account
    im.second.vertex_before = vertex_before;
    im.second.vertex_after = vertex_after;
    im.second.free_vertex = free_vertex;
    im.second.score = rewards[free_vertex];
    im.second.heuristic_value = im.second.score / im.second.cost_increase;
  } else {
    // Insert can't be performed
    im.first = false;
  }
  return im;
}

// TODO: Fill me in
// TODO: Docstring
Vector<Chromosome> InitialisePopulation(const Properties &properties,
                                        const Vector<double_t> &rewards,
                                        const Vector<double_t> &probs,
                                        const Matrix<double_t> &costs,
                                        rng::RandomNumberGenerator &rng) {
  Vector<Chromosome> pop;
  pop.reserve(properties.population_size);
  for (size_t i = 0; i < properties.population_size; ++i) {
    pop.push_back(GenerateChromosome(properties, rewards, probs, costs, rng));
  }
  return pop;
}

// TODO: Fill me in
// TODO: Docstring
Chromosome TournamentSelect(const Vector<Chromosome> &pop,
                            size_t tournament_size,
                            rng::RandomNumberGenerator &rng) {
  Chromosome best;
  double_t best_fitness = -std::numeric_limits<double_t>::infinity();
  Vector<size_t> indices =
      rng.SampleRandomIndices(0, pop.size() - 1, tournament_size);
  for (const size_t &idx : indices) {
    if (pop[idx].fitness > best_fitness) {
      best_fitness = pop[idx].fitness;
      best = pop[idx];
    }
  }
  return best;
}

// TODO: Fill me in
// TODO: Docstring
Vector<VertexId> GetCommonVertices(const Path &p1, const Path &p2) {
  Vector<VertexId> intersection;

  if (p1.size() < 3 || p2.size() < 3) {
    return intersection;
  }

  std::set<VertexId> p1_seen{p1.begin() + 1, p1.end() - 1};
  std::set<VertexId> p2_seen{p2.begin() + 1, p2.end() - 1};

  std::set_intersection(p1_seen.begin(), p1_seen.end(), p2_seen.begin(),
                        p2_seen.end(), std::back_inserter(intersection));

  return intersection;
}

// TODO: Docstring
Chromosome GenerateCXChromosome(const Chromosome &p1, const Chromosome &p2,
                                VertexId split_vertex) {
  Chromosome ret;
  Path::const_iterator p1_vertex_pos =
      std::find(p1.p.begin(), p1.p.end(), split_vertex);
  Path::const_iterator p2_vertex_pos =
      std::find(p2.p.begin(), p2.p.end(), split_vertex);

  ret.free_vertices = p1.free_vertices;

  ret.p.insert(ret.p.end(), p1.p.begin(), p1_vertex_pos);

  Path::const_iterator p1_iter = p1_vertex_pos + 1;
  for (; p1_iter < p1.p.end() - 1; ++p1_iter) {
    ret.free_vertices.push_back(*p1_iter);
  }

  std::set<VertexId> ret_seen(ret.p.begin(), ret.p.end());
  for (; p2_vertex_pos < p2.p.end(); ++p2_vertex_pos) {
    auto insert_res = ret_seen.insert(*p2_vertex_pos);
    if (insert_res.second) {
      ret.p.push_back(*p2_vertex_pos);
      RemoveVertex(ret.free_vertices, *p2_vertex_pos);
    }
  }
  return ret;
}

// TODO: Fill me in
// TODO: Docstring
void Crossover(Chromosome &p1, Chromosome &p2, const Properties &properties,
               const Vector<double_t> &rewards, const Vector<double_t> &probs,
               const Matrix<double_t> &costs, rng::RandomNumberGenerator &rng) {
  Vector<VertexId> common_vertices = GetCommonVertices(p1.p, p2.p);

  if (common_vertices.empty()) {
    return;
  }

  size_t random_vertex_idx =
      rng.GenerateUniformInt(0, common_vertices.size() - 1);
  VertexId vertex = common_vertices[random_vertex_idx];

  Chromosome off1 = GenerateCXChromosome(p1, p2, vertex);
  Chromosome off2 = GenerateCXChromosome(p2, p1, vertex);

  EvaluateChromosome(off1, properties, rewards, probs, costs);
  EvaluateChromosome(off2, properties, rewards, probs, costs);

  if (less_equal(off1.total_cost, properties.maximum_cost)) {
    p1 = off1;
  }

  if (less_equal(off2.total_cost, properties.maximum_cost)) {
    p2 = off2;
  }
}

// TODO: Fill me in
// TODO: Docstring
void Mutate(Chromosome &c, const Properties &properties,
            const Vector<double_t> &rewards, const Vector<double_t> &probs,
            const Matrix<double_t> &costs, rng::RandomNumberGenerator &rng) {
  double_t add_prob = rng.GenerateUniformDouble(0.0, 1.0);
  if (add_prob < properties.mutate_add_prob) {
    MutateAdd(c, properties, rewards, probs, costs, rng);
  } else {
    MutateRemove(c, properties, rewards, probs, costs, rng);
  }
}

// TODO: Fill me in
// TODO: Docstring
void MutateRemove(Chromosome &c, const Properties &properties,
                  const Vector<double_t> &rewards,
                  const Vector<double_t> &probs, const Matrix<double_t> &costs,
                  rng::RandomNumberGenerator &rng) {
  // For each internal vertex calculate the loss over travel decrease
  // Remove the minimum one

  // No vertex to remove
  if (c.p.size() < 3) {
    return;
  }

  double_t min_loss = std::numeric_limits<double_t>::infinity();
  Path::iterator min_vertex = c.p.end();
  Path::iterator it = c.p.begin() + 1;
  for (; it < c.p.end() - 1; ++it) {
    double_t cost_reduction = costs[*(it - 1)][*it] + costs[*it][*(it + 1)] -
                              costs[*(it - 1)][*(it + 1)];
    double_t reward_reduction = rewards[*it];
    double_t loss = reward_reduction / cost_reduction;
    if (loss < min_loss) {
      min_loss = loss;
      min_vertex = it;
    }
  }
  c.free_vertices.push_back(*min_vertex);
  c.p.erase(min_vertex);
  EvaluateChromosome(c, properties, rewards, probs, costs);
}

// TODO: Fill me in
// TODO: Docstring
void MutateAdd(Chromosome &c, const Properties &properties,
               const Vector<double_t> &rewards, const Vector<double_t> &probs,
               const Matrix<double_t> &costs, rng::RandomNumberGenerator &rng) {
  // Choose a random vertex from the free ones.
  // Form all insert moves.
  // Choose the best one (or maybe a random one?)
  // If feasible to insert do it
  if (c.free_vertices.empty()) {
    return;
  }

  size_t random_vertex_idx =
      rng.GenerateUniformInt(0, c.free_vertices.size() - 1);
  VertexId random_vertex = c.free_vertices[random_vertex_idx];
  Path::iterator it = c.p.begin();
  InsertMove best_insert;
  best_insert.heuristic_value = -std::numeric_limits<double_t>::infinity();
  best_insert.vertex_before = c.p.end();
  for (; it < c.p.end() - 1; ++it) {
    InsertMoveRet ret =
        GenerateInsertMove(random_vertex, it, it + 1, properties, rewards,
                           probs, costs, c.total_cost);
    if (ret.first) {
      if (ret.second.heuristic_value > best_insert.heuristic_value) {
        best_insert = ret.second;
      }
    }
  }
  if (best_insert.vertex_before != c.p.end()) {
    c.p.insert(best_insert.vertex_after, best_insert.free_vertex);
    RemoveVertex(c.free_vertices, best_insert.free_vertex);
    EvaluateChromosome(c, properties, rewards, probs, costs);
  }
}

// TODO: Fill me in
// TODO: Docstring
void EvaluateChromosome(Chromosome &c, const Properties &properties,
                        const Vector<double_t> &rewards,
                        const Vector<double_t> &probs,
                        const Matrix<double_t> &costs) {
  c.total_cost = CalculateTotalCost(c.p, costs);
  c.expected_cost = CalculateExpectedCost(c.p, costs, probs);
  c.expected_reward = CalculateExpectedReward(c.p, rewards, probs);
  c.fitness =
      c.expected_reward - properties.cost_per_time_unit * c.expected_cost;
}
}  // namespace pop_ga

