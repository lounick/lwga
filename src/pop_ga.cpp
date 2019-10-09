#include "pop_ga.h"
#include <algorithm>

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
Chromosome GenerateChromosome(const Properties &properties,
                              const Matrix<double_t> &costs,
                              rng::RandomNumberGenerator &rng) {
  switch (properties.generation_method) {
    case GenerationMethod::RANDOM:
      return GenerateRandomChromosome(properties, costs, rng);
    case GenerationMethod::GRASP:
      return GenerateGRASPChromosome(properties, costs, rng);
  }
}

// TODO: Docstring
Chromosome GenerateRandomChromosome(const Properties &properties,
                                    const Matrix<double_t> &costs,
                                    rng::RandomNumberGenerator &rng) {
  Chromosome c;
  // Insert start vertex
  // While not done
  //   Get random vertex
  //   If less than max_cost
  //     Insert
  //   Else
  //     done
  // Insert end vertex
  c.cost = 0.0;
  c.reward = 0.0;
  c.fitness = 0.0;
  c.free_vertices = Vector<VertexId>(costs.size());
  std::iota(c.free_vertices.begin(), c.free_vertices.end(), 0);
  c.free_vertices.erase(std::remove(
      c.free_vertices.begin(), c.free_vertices.end(), properties.start_id));
  c.free_vertices.erase(std::remove(c.free_vertices.begin(),
                                    c.free_vertices.end(), properties.end_id));
  c.p.reserve(costs.size());
  c.p.push_back(properties.start_id);
  bool done = false;
  while (!done) {
    size_t random_idx = rng.GenerateUniformInt(0, c.free_vertices.size() - 1);
    VertexId vertex = c.free_vertices[random_idx];
    double_t total_cost =
        c.cost + costs[c.p.back()][vertex] + costs[vertex][properties.end_id];
    if (total_cost < properties.maximum_cost ||
        logically_equal(total_cost, properties.maximum_cost)) {
      c.cost += costs[c.p.back()][vertex];
      c.p.push_back(vertex);
      c.free_vertices.erase(
          std::remove(c.free_vertices.begin(), c.free_vertices.end(), vertex));
    } else {
      done = true;
    }
    if (c.free_vertices.empty()) {
      done = true;
    }
  }
  c.cost += costs[c.p.back()][properties.end_id];
  c.p.push_back(properties.end_id);
  return c;
}

// TODO: Docstring
Chromosome GenerateGRASPChromosome(const Properties &properties,
                                   const Matrix<double_t> &costs,
                                   rng::RandomNumberGenerator &rng) {
  Chromosome c;
  return c;
}
}  // namespace pop_ga

