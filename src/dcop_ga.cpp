//
// Created by nick on 12/01/18.
//

#include "dcop_ga.h"

namespace dcop_ga {

void Chromosome::calculate_cost(Matrix<Matrix<double_t>>&dubins_cost_mat) {
  cost = 0.0;
  for (size_t i = 1; i < path.size(); ++i){
    cost += dubins_cost_mat[path[i-1]][path[i]][angles[i-1]][angles[i]];
  }
}

void Chromosome::evaluate_chromosome(Matrix<Matrix<double_t>>&dubins_cost_mat,
                                     Matrix<double_t> &cost_mat,
                                     std::vector<double_t> &rewards) {
  fitness = 0;
  static std::vector<uint_fast32_t> vertices;
  if(vertices.size() == 0) {
    vertices.reserve(cost_mat.size());
    for (uint_fast32_t i = 0; i < cost_mat.size(); ++i)
      vertices.push_back(i);
  }

  size_t pathsize = path.size() - 1;
  size_t fsize = free_vertices.size();

  for (size_t i = 1; i < pathsize; ++i) {
    double extras = 0;
    uint_fast32_t vertex = path[i];
    if (seen_vertices.count(vertex)) {
      for (size_t j = 0; j < fsize; ++j) {
        double dist = cost_mat[vertex][free_vertices[j]];
        if (dist < 2) {
          extras += std::exp((log(0.01)/2) * dist)*rewards[free_vertices[j]]; //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
        }
      }
      fitness += rewards[vertex] + extras;
    }
  }
  fitness = pow(fitness,3)/cost;
}

void Chromosome::mutate(
    Matrix<Matrix<double_t>>&dubins_cost_mat, Matrix<double_t> &cost_mat,
    Vector<double_t> &std_angles,
    std::vector<double_t> &rewards, double_t max_cost, std::mt19937 &g){
  Chromosome mutated(*this);
  std::tie(mutated.path, mutated.angles, mutated.cost) =
      dubins_two_opt(
          dubins_cost_mat, std_angles, mutated.nodes, mutated.rho, mutated.path,
          mutated.angles, mutated.cost);
  mutated.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
  for (uint_fast32_t iter = 0; iter < 10; ++iter){
    if (std::generate_canonical<double, 10>(g) < 0.99) {
      //Try to add
      if (mutated.cost > 0.9*max_cost) {
        std::uniform_int_distribution<> dis(1,mutated.path.size()-2);
        size_t rand_vertex_idx = dis(g);
        uint_fast32_t rand_vertex = mutated.path[rand_vertex_idx];
        uint_fast32_t rand_angle = mutated.angles[rand_vertex_idx];
        Vector<uint_fast32_t> available_vertices;
        available_vertices.reserve(free_vertices.size());
        for (uint_fast32_t i:mutated.free_vertices) {
          if (cost_mat[rand_vertex][i] < 2) {
            available_vertices.push_back(i);
          }
        }

        double_t best_cost = mutated.cost;
        Vector<uint_fast32_t> best_path = mutated.path;
        Vector<uint_fast32_t> best_angles = mutated.angles;
        double_t best_fitness = mutated.fitness;
        uint_fast32_t best_vertex = rand_vertex;
        uint_fast32_t prev_vertex = mutated.path[rand_vertex_idx - 1];
        uint_fast32_t prev_angle = mutated.angles[rand_vertex_idx - 1];
        uint_fast32_t next_vertex = mutated.path[rand_vertex_idx + 1];
        uint_fast32_t next_angle = mutated.angles[rand_vertex_idx + 1];
        double_t cost_removed = mutated.cost
            - dubins_cost_mat[prev_vertex][rand_vertex][prev_angle][rand_angle]
            - dubins_cost_mat[rand_vertex][next_vertex][rand_angle][next_angle];
        for (uint_fast32_t i:available_vertices) {
          Chromosome tmp_c(mutated);
//          tmp_c.path.reserve(cost_mat.size());
//          tmp_c.path = mutated.path;
//          tmp_c.angles.reserve(cost_mat.size());
//          tmp_c.angles = mutated.angles;
//          tmp_c.path.erase(tmp_c.path.begin() + rand_vertex_idx);
//          tmp_c.path.insert(tmp_c.path.begin() + rand_vertex_idx, i);
          tmp_c.path[rand_vertex_idx] = i;
          uint8_t angles_num = 5;
          Vector<uint_fast32_t> p_angles, n_angles;
          p_angles.reserve(angles_num);
          n_angles.reserve(angles_num);
          for (int32_t angle_cnt = -2; angle_cnt < 3; ++angle_cnt) {
            if (int32_t(prev_angle) + angle_cnt < 0) {
              p_angles.push_back(std_angles.size() + prev_angle + angle_cnt);
            } else if (prev_angle + angle_cnt >= std_angles.size()) {
              p_angles.push_back(prev_angle + angle_cnt - std_angles.size());
            } else {
              p_angles.push_back(prev_angle + angle_cnt);
            }
            if (int32_t(next_angle) + angle_cnt < 0) {
              n_angles.push_back(std_angles.size() + next_angle + angle_cnt);
            } else if (next_angle + angle_cnt >= std_angles.size()) {
              n_angles.push_back(next_angle + angle_cnt - std_angles.size());
            } else {
              n_angles.push_back(next_angle + angle_cnt);
            }
          }
          double_t best_cost_added = std::numeric_limits<double_t>::max();
          uint_fast32_t best_pa, best_a, best_na;
          for (uint_fast32_t pa:p_angles) {
            for (uint_fast32_t na:n_angles) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a) {
                double_t cost_added =
                    dubins_cost_mat[prev_vertex][i][pa][a]
                        + dubins_cost_mat[i][next_vertex][a][na];
                if (cost_added < best_cost_added) {
                  best_cost_added = cost_added;
                  best_pa = pa;
                  best_na = na;
                  best_a = a;
                }
              }
            }
          }
          tmp_c.angles[rand_vertex_idx-1] = best_pa;
          tmp_c.angles[rand_vertex_idx] = best_a;
          tmp_c.angles[rand_vertex_idx+1] = best_na;
//          std::tie(tmp_c.path, tmp_c.angles, tmp_c.cost) =
//              dubins_two_opt(dubins_cost_mat, std_angles, tmp_c.nodes,
//                             tmp_c.rho, tmp_c.path, tmp_c.angles, tmp_c.cost);
          tmp_c.calculate_cost(dubins_cost_mat);
          if (tmp_c.cost < max_cost || logically_equal(tmp_c.cost, max_cost)) {
            tmp_c.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
            if (tmp_c.fitness > best_fitness) {
              best_fitness = tmp_c.fitness;
              best_cost = tmp_c.cost;
              best_path = tmp_c.path;
              best_angles = tmp_c.angles;
              best_vertex = i;
            }
          }
        }
        mutated.path = best_path;
        mutated.angles = best_angles;
        mutated.cost = best_cost;
        mutated.fitness = best_fitness;
        if (best_vertex != rand_vertex) {
          mutated.free_vertices.push_back(rand_vertex);
          mutated.free_vertices.erase(std::find(
              mutated.free_vertices.begin(), mutated.free_vertices.end(), best_vertex));
          mutated.seen_vertices.insert(best_vertex);
        }
      } else {
        if (free_vertices.size() == 0)
          continue;
        // Choose the vertex that has the most free neighbours
        // Choose the neighbour that has the most free neighbours
        // Try to add the neighbour before and after the vertex
        // If max cost is not violated
        // Choose addition with minimal cost increase
      }
    } else {
      //Remove
      if (mutated.cost >= 0.9 * max_cost) {
        size_t to_remove = 0;
        double min_loss = std::numeric_limits<double>::infinity();
        for (size_t i = 1; i < mutated.path.size() - 1; ++i) {
          uint_fast32_t v = mutated.path[i];
          uint_fast32_t a = mutated.angles[i];
          uint_fast32_t pv= mutated.path[i-1];
          uint_fast32_t pa = mutated.angles[i-1];
          uint_fast32_t nv = mutated.path[i+1];
          uint_fast32_t na = mutated.angles[i+1];

          double travel_decrease =dubins_cost_mat[pv][v][pa][a]
              + dubins_cost_mat[v][nv][a][na] - dubins_cost_mat[pv][nv][pa][na];
          double loss = rewards[v];
          double extras = 0;
          for (size_t j = 0; j < free_vertices.size(); ++j) {
            if (cost_mat[mutated.path[i]][free_vertices[j]] < 2) {
              extras += std::exp((log(0.01)/2) * cost_mat[mutated.path[i]][free_vertices[j]])*rewards[free_vertices[j]];
            }
          }
          loss += extras;

          if (travel_decrease != 0) {
            loss /= travel_decrease;
          }
          else {
            loss = std::numeric_limits<double>::infinity();
          }

          if (loss <= min_loss) {
            min_loss = loss;
            to_remove = i;
          }
        }
        if (to_remove > 0) {
          mutated.free_vertices.push_back(mutated.path[to_remove]);
          mutated.seen_vertices.erase(to_remove);
          mutated.path.erase(mutated.path.begin() + to_remove);
          mutated.calculate_cost(dubins_cost_mat);
          mutated.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
        }
      }
    }
  }
  path = mutated.path;
  angles = mutated.angles;
  cost = mutated.cost;
  fitness = mutated.fitness;
  free_vertices = mutated.free_vertices;
  seen_vertices = mutated.seen_vertices;
}

void expand_neighbours(std::vector<uint_fast32_t> &neighbours,
                       std::vector<uint_fast32_t> &checked,
                       const std::vector<uint_fast32_t> &vertices,
                       const Matrix<double> &cost_mat) {

  std::unordered_set<uint_fast32_t> nset;
  nset.reserve(vertices.size());

  checked.insert(checked.end(), neighbours.begin(), neighbours.end());
  std::vector<uint_fast32_t> free_vertices;
  free_vertices.reserve(vertices.size());
  std::sort(checked.begin(), checked.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      checked.begin(),
                      checked.end(),
                      std::back_inserter(free_vertices));

  std::vector<uint_fast32_t>::iterator it;
  for (it = neighbours.begin(); it != neighbours.end(); ++it) {
    std::vector<uint_fast32_t>::iterator free_it;
    for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
      double dist = cost_mat[*it][*free_it];
      if (dist < 2) {
        nset.insert(*free_it);
      }
    }
  }
  neighbours.assign(nset.begin(), nset.end());
}

Chromosome generate_chromosome(std::shared_ptr<Vector<Point2D>> nodes,
                               Vector<double_t> std_angles,
                               double_t rho,
                               double_t max_cost,
                               uint_fast32_t idx_start,
                               uint_fast32_t idx_finish,
                               Matrix<Matrix<double_t>>&dubins_cost_mat,
                               const Matrix<double_t> &cost_mat,
                               std::mt19937 &g) {
  Chromosome c;
  c.nodes = nodes;
  c.rho = rho;
  c.fitness = 0;
  c.cost = 0;
  c.all_vertices = std::vector<uint_fast32_t>(cost_mat.size());
  std::iota(c.all_vertices.begin(), c.all_vertices.end(), 0);
  c.seen_vertices.reserve(cost_mat.size());
  c.free_vertices = std::vector<uint_fast32_t>(
      c.all_vertices.begin(), c.all_vertices.end());
  c.free_vertices.erase(
      std::remove(c.free_vertices.begin(), c.free_vertices.end(), idx_start));
  c.free_vertices.erase(
      std::remove(c.free_vertices.begin(), c.free_vertices.end(), idx_finish));
  c.path.reserve(cost_mat.size());
  c.path.push_back(idx_start);
  c.angles.reserve(cost_mat.size());
  c.angles.push_back(0.0);

  std::vector<uint_fast32_t> used_vertices;
  used_vertices.reserve(cost_mat.size());
  double_t path_cost = 0;
  bool done = false;
  std::vector<uint_fast32_t> checked;
  checked.reserve(cost_mat.size());
  checked.clear();
  while(!done){
    std::vector<uint_fast32_t> neighbours;
    neighbours.reserve(cost_mat.size());
    neighbours.clear();
    neighbours.push_back(c.path.back());
    std::vector<uint_fast32_t> available;
    available.reserve(cost_mat.size());
    available.clear();
    bool neighbours_found = false;
    while (!neighbours_found) {
      expand_neighbours(neighbours, checked, c.all_vertices, cost_mat);
      if (neighbours.size() == 0) { neighbours_found = true; }
      else {
        std::sort(neighbours.begin(), neighbours.end());
        std::set_intersection(neighbours.begin(),
                              neighbours.end(),
                              c.free_vertices.begin(),
                              c.free_vertices.end(),
                              std::back_inserter(available));
        if (available.size() > 0) {
          neighbours_found = true;
        }
      }
    }
    if (available.size() > 0) {
      std::uniform_real_distribution<> dis(0, 1);
      std::vector<double_t> distance_weights;
      std::vector<double_t> neighbour_weights;
      distance_weights.reserve(available.size());
      neighbour_weights.reserve(available.size());
      std::vector<double_t> combined;
      combined.reserve(available.size());
      std::vector<double_t> cdf;
      cdf.reserve(available.size());
      double_t tot_dist = 0;
      double_t tot_neighbours = 0;
      if (c.seen_vertices.size() > 0) {
        std::vector<uint_fast32_t>::iterator avail_it;
        for (avail_it = available.begin(); avail_it != available.end(); ++avail_it) {
          double mean_dist = 0;
          double free_neighbours = 0;
          std::unordered_set<uint_fast32_t>::iterator used_it;
          for (used_it = c.seen_vertices.begin(); used_it != c.seen_vertices.end(); ++used_it) {
            mean_dist += cost_mat[*avail_it][*used_it];
          }
          std::vector<uint_fast32_t>::iterator free_it;
          for (free_it = c.free_vertices.begin(); free_it != c.free_vertices.end(); ++free_it) {
            if (cost_mat[*avail_it][*free_it] < 2) {
              ++free_neighbours;
            }
          }
          mean_dist = mean_dist / double(c.seen_vertices.size());
          tot_dist += mean_dist;
          tot_neighbours += free_neighbours;
          distance_weights.push_back(mean_dist);
          neighbour_weights.push_back(free_neighbours);
        }

        std::vector<double>::iterator dist_weight_it = distance_weights.begin();
        std::vector<double>::iterator neigh_weight_it = neighbour_weights.begin();
        double tot_combined = 0;
        for (; dist_weight_it != distance_weights.end(); ++dist_weight_it, ++neigh_weight_it) {
          combined.push_back(*dist_weight_it / double_t(tot_dist) + *neigh_weight_it / double_t(tot_neighbours));
          tot_combined += *dist_weight_it / double_t(tot_dist) + *neigh_weight_it / double_t(tot_neighbours);
        }
        std::vector<double>::iterator combined_it;
        cdf.push_back(combined.front() / tot_combined);
        for (combined_it = combined.begin() + 1; combined_it != combined.end(); ++combined_it) {
          cdf.push_back(cdf.back() + (*combined_it / tot_combined));
        }
      } else {
        //No prior all is equally probable;
        double weight = 1.0 / available.size();
        cdf.push_back(weight);
        std::vector<uint_fast32_t>::iterator avail_it;
        for (avail_it = available.begin() + 1; avail_it != available.end(); ++avail_it) {
          cdf.push_back(cdf.back() + weight);
        }
      }
      double prob = dis(g);
      auto vertex_pos = std::lower_bound(cdf.cbegin(), cdf.cend(), prob);
      uint_fast32_t next_vertex = available[std::distance(cdf.cbegin(), vertex_pos)];
      uint_fast32_t min_angle_next, min_angle_end;
      double_t min_distance = std::numeric_limits<double_t>::max();
      for (size_t next_angle_idx = 0; next_angle_idx < std_angles.size(); ++next_angle_idx) {
        for (size_t end_angle_idx = 0; end_angle_idx < std_angles.size(); ++end_angle_idx){
          double_t dist = 0.0;
          dist += dubins_cost_mat[c.path.back()][next_vertex][c.angles.back()][next_angle_idx];
          dist += dubins_cost_mat[next_vertex][idx_finish][next_angle_idx][end_angle_idx];
          if(dist < min_distance){
            min_distance = dist;
            min_angle_next = next_angle_idx;
            min_angle_end = end_angle_idx;
          }
        }
      }

      double next_vertex_end_cost = path_cost + min_distance;
      if (next_vertex_end_cost < max_cost) {
        path_cost += dubins_cost_mat[c.path.back()][next_vertex][c.angles.back()][min_angle_next];
        c.path.push_back(next_vertex);
        c.angles.push_back(min_angle_next);
        c.free_vertices.erase(std::remove(c.free_vertices.begin(), c.free_vertices.end(), next_vertex));
        c.seen_vertices.insert(next_vertex);
        c.cost = path_cost;
      } else {
        done = true;
      }
    } else {
      done = true;
    }
  }

  uint_fast32_t min_angle_end;
  double_t min_distance = std::numeric_limits<double_t>::max();
  for (size_t end_angle_idx = 0;
       end_angle_idx < std_angles.size(); ++end_angle_idx) {
    double_t dist = dubins_cost_mat[c.path.back()][idx_finish][c.angles.back()][end_angle_idx];
    if (dist < min_distance){
      min_distance = dist;
      min_angle_end = end_angle_idx;
    }
  }
  c.cost += min_distance;
  c.path.push_back(idx_finish);
  c.angles.push_back(min_angle_end);
  return c;
}

Chromosome tournament_select(std::vector<Chromosome> &population,
                             uint_fast32_t tour_size,
                             std::mt19937 &g){
  Chromosome best_chromosome;
  if (population.size() < tour_size) {
    best_chromosome = population[0];
  } else {
    std::vector<uint> indices(population.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), g);
    double_t max_fitness = std::numeric_limits<double_t>::min();
    for (size_t i = 0; i < tour_size; i++) {
      if (population[indices[i]].fitness > max_fitness) {
        max_fitness = population[indices[i]].fitness;
        best_chromosome = population[indices[i]];
      }
    }
  }
  return best_chromosome;
}

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<Matrix<double_t>>&dubins_cost_mat,
                Matrix<double> &cost_mat,
                Vector<double_t> &std_angles,
                std::vector<double> &rewards,
                double &max_cost){
  // Get hash of thread id for the seed of the generator.
  std::hash<std::thread::id> hasher;
  static thread_local std::mt19937 g(hasher(std::this_thread::get_id()));

  for(size_t idx:indices)
    pop[idx].mutate(dubins_cost_mat, cost_mat, std_angles, rewards, max_cost, g);
}

Chromosome ga_dcop(std::shared_ptr<Vector<Point2D>> nodes,
                   Vector<double_t> std_angles,
                   double_t rho,
                   Matrix<Matrix<double_t>>&dubins_cost_mat,
                   Matrix<double_t> &cost_mat,
                   std::vector<double_t> &rewards,
                   double_t max_cost,
                   uint_fast32_t idx_start,
                   uint_fast32_t idx_finish,
                   std::mt19937 &g){
  /*
   * Initialise population
   * While gen < max_gen
   *  Select new pop
   *  Cx
   *  Mutate
   * Select fittest
  */

  uint_fast32_t pop_size = 200;
  int tour_size = 3;
  int max_gen = 75;

  // Initialise population
  std::vector<Chromosome> pop;
  pop.reserve(pop_size);

  for (uint_fast32_t i = 0; i < pop_size; ++i) {
//    std::cout << "Generating chromosome: " << i << std::endl;
    Chromosome c = generate_chromosome(
        nodes, std_angles, rho, max_cost,
        idx_start, idx_finish, dubins_cost_mat, cost_mat, g);
//    std::tie(c.angles, c.cost) = straighten_path(dubins_cost_mat, nodes, rho, c.path, c.angles, c.cost);
    std::tie(c.path, c.angles, c.cost) =
        dubins_two_opt(dubins_cost_mat, std_angles, nodes, rho, c.path, c.angles, c.cost);
//    std::tie(c.angles, c.cost) = straighten_path(dubins_cost_mat, nodes, rho, c.path, c.angles, c.cost);
    c.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
    pop.push_back(c);
  }

  for (int gen = 0; gen < max_gen; ++gen) {
    // Select new populations
    std::vector<Chromosome> new_pop;
    new_pop.reserve(pop_size);

    for (uint_fast32_t i = 0; i < pop_size; ++i) {
      new_pop.push_back(tournament_select(pop, tour_size, g));
    }

    // Cx
//    for (int i = 0; i < 20; ++i) {
//      std::vector<size_t> indices = get_population_sample(new_pop.size(), 2, g);
//      std::pair<Chromosome, Chromosome> cx_ret = cx(new_pop[indices[0]], new_pop[indices[1]], cost_mat, max_cost, g);
//      new_pop[indices[0]] = cx_ret.first;
//      new_pop[indices[1]] = cx_ret.second;
//      new_pop[indices[0]].evaluate_chromosome(cost_mat, rewards);
//      new_pop[indices[1]].evaluate_chromosome(cost_mat, rewards);
//    }
    // Mutate
    std::vector<size_t> indices = get_population_sample(new_pop.size(), 50, g);
    uint_fast32_t M = 6; //number of cores
    uint_fast32_t chunk_size = indices.size()/M;

//    for (size_t idx : indices) {
//      new_pop[idx].mutate(dubins_cost_mat, cost_mat, std_angles, rewards, max_cost, g);
//    }
    std::vector< std::future<void> > future_v;
    for (uint_fast32_t thread_count = 0; thread_count < M; ++thread_count) {
      //std::launch::deferred|std::launch::async;
      std::vector<size_t > tmpv(indices.begin()+thread_count*chunk_size,indices.begin()+(thread_count+1)*chunk_size);
      future_v.push_back(std::async(std::launch::async, par_mutate, tmpv, std::ref(new_pop), std::ref(dubins_cost_mat), std::ref(cost_mat), std::ref(std_angles), std::ref(rewards), std::ref(max_cost)));
    }
    if(indices.size()%M != 0){
      std::vector<size_t> tmpv(indices.begin()+(M)*chunk_size,indices.end());
      future_v.push_back(std::async(std::launch::async, par_mutate, tmpv, std::ref(new_pop), std::ref(dubins_cost_mat), std::ref(cost_mat), std::ref(std_angles), std::ref(rewards), std::ref(max_cost)));
    }
    for(auto &f: future_v){
      f.get();
    }

    pop = new_pop;
  }
  std::sort(pop.begin(), pop.end(), [](Chromosome c1, Chromosome c2) { return c1.fitness > c2.fitness; });
  Chromosome best = pop[0];
  return best;
}

}