#include "cop_ga.h"

void Chromosome::calculate_cost(Matrix<double> &cost_mat) {
  cost = get_path_cost(path, cost_mat);
}

void Chromosome::evaluate_chromosome(Matrix<double> &cost_mat, std::vector<double> &rewards) {
  fitness = 0;
  static std::vector<uint_fast32_t> vertices;
  if(vertices.size() == 0) {
    vertices.reserve(cost_mat.size());
    for (uint_fast32_t i = 0; i < cost_mat.size(); ++i)
      vertices.push_back(i);
  }
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices = path;
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));
  std::unordered_set<uint_fast32_t> seen;
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;

  size_t pathsize = path.size() - 1;
  size_t fsize = free_vertices.size();

  for (size_t i = 1; i < pathsize; ++i) {
    double extras = 0;
    uint_fast32_t vertex = path[i];
    insert_ret = seen.insert(vertex);
    if (insert_ret.second) {
      for (size_t j = 0; j < fsize; ++j) {
        double dist = cost_mat[vertex][free_vertices[j]];
        if (dist < 2) {
          extras += std::exp(-2 * dist); //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
        }
      }
      fitness += rewards[vertex] + extras;
    }
  }
  fitness = pow(fitness,3)/cost;
}

void Chromosome::mutate(Matrix<double> &cost_mat, std::vector<double> &rewards, double max_cost, std::mt19937 &g) {
  Chromosome mutated;
  mutated.path.reserve(cost_mat.size());
  std::pair<Path, double> two_opt_return;
  two_opt_return = two_opt(path, cost_mat);
  mutated.path = two_opt_return.first;
  mutated.cost = two_opt_return.second;
  mutated.evaluate_chromosome(cost_mat, rewards);

  std::unordered_set<uint_fast32_t> seen;
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_return;
  Path new_path;
  for (size_t i = 0; i < mutated.path.size(); ++i) {
    insert_return = seen.insert(mutated.path[i]);
    if (insert_return.second) {
      new_path.push_back(mutated.path[i]);
    }
  }
  double cost_difference = mutated.path.size() - new_path.size();
  mutated.cost -= cost_difference; //TODO: This assumes a fixed cost per task. Maybe that's not the case in the real world.
  mutated.path = new_path;
  mutated.evaluate_chromosome(cost_mat, rewards);

  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices = mutated.path;
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));

  for (uint_fast32_t iter = 0; iter < 10; ++iter) {
    if (std::generate_canonical<double, 10>(g) < 0.9) {
      if (mutated.cost >= 0.95 * max_cost) {
        //TODO: This is bound to have different effect in different grid sizes and budgets.
        //TODO: Should come up with something including the average travel cost or something smarter.
        /*
         * Pick random vertex from solution, apart from start and end.
         * Check it's free neighbours
         * Choose neighbour that maximises the fitness after removing vertex
         * If new fitness is better than old then accept solution
         */
        std::uniform_int_distribution<> dis(1, mutated.path.size() - 2);
        size_t rand_vertex_idx = dis(g);
        uint_fast32_t rand_vertex = mutated.path[rand_vertex_idx];
        std::vector<uint_fast32_t> available_vertices;
        for (uint_fast32_t i:free_vertices) {
          if (cost_mat[rand_vertex][i] < 2) { //TODO: fix the constant sensing distace
            available_vertices.push_back(i);
          }
        }

        double best_cost = mutated.cost;
        std::vector<uint_fast32_t> best_path = mutated.path;
        double best_fitness = mutated.fitness;
        uint_fast32_t best_vertex = rand_vertex;
        uint_fast32_t prev_vertex = mutated.path[rand_vertex_idx - 1];
        uint_fast32_t next_vertex = mutated.path[rand_vertex_idx + 1];
        double cost_removed = mutated.cost
            - cost_mat[prev_vertex][rand_vertex]
            - cost_mat[rand_vertex][next_vertex];
        for (uint_fast32_t i:available_vertices) {
          Chromosome tmp_c;
          tmp_c.path.reserve(cost_mat.size());
          tmp_c.path = mutated.path;
          tmp_c.cost = cost_removed
              + cost_mat[prev_vertex][i]
              + cost_mat[i][next_vertex];
          tmp_c.path.erase(tmp_c.path.begin() + rand_vertex_idx);
          tmp_c.path.insert(tmp_c.path.begin() + rand_vertex_idx, i);
          if (tmp_c.cost <= max_cost) {
            tmp_c.evaluate_chromosome(cost_mat, rewards);
            if (tmp_c.fitness > best_fitness) {
              best_fitness = tmp_c.fitness;
              best_cost = tmp_c.cost;
              best_path = tmp_c.path;
              best_vertex = i;
            }
          }
        }
        mutated.path = best_path;
        mutated.fitness = best_fitness;
        mutated.cost = best_cost;
        if (best_vertex != rand_vertex) {
          free_vertices.push_back(rand_vertex);
          free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), best_vertex));
        }
      }
      else {
        if (free_vertices.size() == 0)
          continue;

        std::uniform_int_distribution<> dis(0, free_vertices.size() - 1);
        size_t vertex_idx = dis(g);
        uint_fast32_t vertex = free_vertices[vertex_idx];

        double best_travel_increase = 0;
        double best_fitness = 0;
        size_t ins_pos = 0;

        for (size_t i = 1; i < mutated.path.size(); ++i) {
          double travel_increase = cost_mat[mutated.path[i - 1]][vertex] + cost_mat[vertex][mutated.path[i]]
              - cost_mat[mutated.path[i - 1]][mutated.path[i]];
          double fit = rewards[vertex]; //TODO: Should we add the extras to the fit?
          if (travel_increase != 0) {
            fit /= travel_increase;
          }
          else {
            fit = std::numeric_limits<double>::infinity();
          }
          if(fit > best_fitness) {
            if (mutated.cost + travel_increase + 1 <= max_cost) {
              best_travel_increase = travel_increase;
              best_fitness = fit;
              ins_pos = i;
            }
          }
        }
        if (ins_pos > 0) {
          mutated.path.insert(mutated.path.begin() + ins_pos, vertex);
          mutated.cost += best_travel_increase + 1;
          mutated.evaluate_chromosome(cost_mat, rewards);
          free_vertices.erase(free_vertices.begin() + vertex_idx);
        }
      }
    }
    else {
      std::unordered_set<uint_fast32_t> duplicates;
      std::unordered_set<uint_fast32_t> seen_vertices;
      std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
      for (uint_fast32_t i:mutated.path) {
        insert_ret = seen_vertices.insert(i);
        if (!insert_ret.second) {
          duplicates.insert(i);
        }
      }
      if (duplicates.size() > 0) {
        // Choose random duplicate
        std::vector<uint_fast32_t> duplicate_v(duplicates.begin(), duplicates.end());
        std::uniform_int_distribution<> dis(0, duplicate_v.size());
        size_t duplcate_idx = dis(g);
        uint_fast32_t duplicate = duplicate_v[duplcate_idx];

        // Find the indices of the duplicate in the path
        std::vector<size_t> indices;
        std::vector<uint_fast32_t>::iterator it = mutated.path.begin();
        while ((it = std::find(it, mutated.path.end(), duplicate)) != mutated.path.end()) {
          indices.push_back(std::distance(mutated.path.begin(), it));
          ++it;
        }

        // Remove the one with the minimum loss
        double min_loss = std::numeric_limits<double>::infinity();
        size_t to_remove = 0;

        for (size_t idx:indices) {
          double travel_decrease = cost_mat[mutated.path[idx - 1]][mutated.path[idx]]
              + cost_mat[mutated.path[idx]][mutated.path[idx + 1]]
              - cost_mat[mutated.path[idx - 1]][mutated.path[idx + 1]];

          double loss = rewards[mutated.path[idx]];

          double extras = 0;
          for (size_t j = 0; j < free_vertices.size(); ++j) {
            if (cost_mat[mutated.path[idx]][free_vertices[j]] < 2) {
              extras += std::exp(-2 * cost_mat[mutated.path[idx]][free_vertices[j]]);
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
            to_remove = idx;
          }
        }

        if (to_remove > 0) {
          mutated.path.erase(mutated.path.begin() + to_remove);
          mutated.calculate_cost(cost_mat);
          mutated.evaluate_chromosome(cost_mat, rewards);
        }
      }
      else {
        if (mutated.cost >= 0.95 * max_cost) {
          size_t to_remove = 0;
          double min_loss = std::numeric_limits<double>::infinity();
          for (size_t i = 1; i < mutated.path.size() - 1; ++i) {
            double travel_decrease =
                cost_mat[mutated.path[i - 1]][mutated.path[i]] + cost_mat[mutated.path[i]][mutated.path[i + 1]]
                    - cost_mat[mutated.path[i - 1]][mutated.path[i + 1]];
            double loss = rewards[mutated.path[i]];
            double extras = 0;
            for (size_t j = 0; j < free_vertices.size(); ++j) {
              if (cost_mat[mutated.path[i]][free_vertices[j]] < 2) {
                extras += std::exp(-2 * cost_mat[mutated.path[i]][free_vertices[j]]);
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
            free_vertices.push_back(mutated.path[to_remove]);
            mutated.path.erase(mutated.path.begin() + to_remove);
            mutated.calculate_cost(cost_mat);
            mutated.evaluate_chromosome(cost_mat, rewards);
          }
        }
      }
    }
  }
  path = mutated.path;
  cost = mutated.cost;
  fitness = mutated.fitness;
}

std::pair<Chromosome, Chromosome> cx(Chromosome &c1, Chromosome &c2, Matrix<double> &cost_mat, double max_cost, std::mt19937 &g) {

  Chromosome off1, off2;
  off1.path.reserve(cost_mat.size());
  off2.path.reserve(cost_mat.size());
  std::unordered_set<int> s1(c1.path.begin(), c1.path.end());
  std::unordered_set<int> s2(c2.path.begin(), c2.path.end());
  std::vector<int> intersection;
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(intersection));

  if (intersection.size() == 0) {
    off1 = c1;
    off2 = c2;
    return std::make_pair(off1, off2);
  }
  else {
    size_t rand_idx = g() % intersection.size();
    std::vector<uint_fast32_t>::iterator idx1 = std::find(c1.path.begin(), c1.path.end(), intersection[rand_idx]);
    std::vector<uint_fast32_t>::iterator idx2 = std::find(c2.path.begin(), c2.path.end(), intersection[rand_idx]);
    std::copy(c1.path.begin(), idx1, std::back_inserter(off1.path));
    std::copy(idx2, c2.path.end(), std::back_inserter(off1.path));
    std::copy(c2.path.begin(), idx2, std::back_inserter(off2.path));
    std::copy(idx1, c1.path.end(), std::back_inserter(off2.path));
    /*
     * 1. Do a two-opt in each offspring.
     * 2. Remove any duplicates
     * 3. Check feasibility and decide what to return.
     */
    std::pair<std::vector<uint_fast32_t>, double> two_opt_ret = two_opt(off1.path, cost_mat);
    off1.path = two_opt_ret.first;
    off1.cost = two_opt_ret.second;
    two_opt_ret = two_opt(off2.path, cost_mat);
    off2.path = two_opt_ret.first;
    off2.cost = two_opt_ret.second;

    std::unordered_set<uint_fast32_t> seen;
    std::vector<uint_fast32_t> new_path;
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_result;
    for (size_t i = 0; i < off1.path.size(); ++i) {
      insert_result = seen.insert(off1.path[i]);
      if (insert_result.second) {
        new_path.push_back(off1.path[i]);
      }
    }
    double cost_difference = off1.path.size() - new_path.size();
    off1.cost -= cost_difference;
    off1.path = new_path;

    seen.clear();
    new_path.clear();
    for (size_t i = 0; i < off2.path.size(); ++i) {
      insert_result = seen.insert(off2.path[i]);
      if (insert_result.second) {
        new_path.push_back(off2.path[i]);
      }
    }
    cost_difference = off2.path.size() - new_path.size();
    off2.cost -= cost_difference;
    off2.path = new_path;

    bool feas1 = (off1.cost <= max_cost);
    bool feas2 = (off2.cost <= max_cost);

    if (feas1 && feas2) {
      return std::make_pair(off1, off2);
    }
    else {
      if (feas1) {
        if (c1.fitness > c2.fitness) {
          off2.path = c1.path;
          off2.cost = c1.cost;
          off2.fitness = c1.fitness;
        }
        else {
          off2.path = c2.path;
          off2.cost = c2.cost;
          off2.fitness = c2.fitness;
        }
        return std::make_pair(off1, off2);
      }
      else {
        if (feas2) {
          if (c1.fitness > c2.fitness) {
            off1.path = c1.path;
            off1.cost = c1.cost;
            off1.fitness = c1.fitness;
          }
          else {
            off1.path = c2.path;
            off1.cost = c2.cost;
            off1.fitness = c2.fitness;
          }
          return std::make_pair(off1, off2);
        }
        else {
          off1.path = c1.path;
          off1.cost = c1.cost;
          off1.fitness = c1.fitness;
          off2.path = c2.path;
          off2.cost = c2.cost;
          off2.fitness = c2.fitness;
          return std::make_pair(off1, off2);
        }
      }
    }
  }
}

Chromosome tournament_select(std::vector<Chromosome> &population, uint_fast32_t tour_size, std::mt19937 &g) {

  std::vector<uint_fast32_t> indices(population.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), g);
  double max_fitness = DBL_MIN;
  Chromosome best_chromosome;

  for (size_t i = 0; i < tour_size; ++i) {
    if (population[indices[i]].fitness > max_fitness) {
      max_fitness = population[indices[i]].fitness;
      best_chromosome = population[indices[i]];
    }
  }

  return best_chromosome;
}

Chromosome generate_chromosome(Matrix<double> &cost_mat, double max_cost, uint_fast32_t idx_start, uint_fast32_t idx_finish, std::mt19937 &g) {

  Chromosome c;
  c.path.reserve(cost_mat.size());
  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);

  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_start));
  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_finish));

  bool done = false;
  double total_cost = 0;
  c.path.push_back(idx_start);

  while (!done) {
    if (vertices.size() == 0)
      break;

    size_t rand_idx = g() % vertices.size();
    uint_fast32_t next_vertex = vertices[rand_idx];

    if (total_cost + cost_mat[c.path.back()][next_vertex] + cost_mat[next_vertex][idx_finish] + 1 <= max_cost) {
      total_cost += cost_mat[c.path.back()][next_vertex] + 1;
      c.path.push_back(next_vertex);
      vertices.erase(std::remove(vertices.begin(), vertices.end(), next_vertex));
    }
    else {
      done = true;
    }
  }

  c.path.push_back(idx_finish);
  return c;
}

std::pair<bool, double> check_feasibility(Chromosome &c, Matrix<double> &cost_mat, double max_cost) {
  double path_cost = get_path_cost(c.path, cost_mat);
  bool is_feasible = path_cost <= max_cost;
  return std::make_pair(is_feasible, path_cost);
}

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                double &max_cost){

  // Get hash of thread id for the seed of the generator.
  std::hash<std::thread::id> hasher;
  static thread_local std::mt19937 g(hasher(std::this_thread::get_id()));

  for(size_t idx:indices)
    pop[idx].mutate(cost_mat, rewards, max_cost, g);
}

Chromosome ga_cop(std::vector<std::vector<double> > &cost_mat,
                  std::vector<double> &rewards,
                  double max_cost,
                  uint_fast32_t idx_start,
                  uint_fast32_t idx_finish,
                  std::mt19937 &g) {

  /*
   * Initialise population
   * While gen < max_gen
   *  Select new pop
   *  Cx
   *  Mutate
   * Select fittest
  */

  uint_fast32_t pop_size = 100;
  int tour_size = 3;
  int max_gen = 50;

  // Initialise population
  std::vector<Chromosome> pop;
  pop.reserve(pop_size);

  for (uint_fast32_t i = 0; i < pop_size; ++i) {
    Chromosome c = generate_chromosome(cost_mat, max_cost, idx_start, idx_finish, g);
    std::pair<std::vector<uint_fast32_t>, double> two_opt_ret = two_opt(c.path, cost_mat);
    c.path = two_opt_ret.first;
    c.cost = two_opt_ret.second;
    c.evaluate_chromosome(cost_mat, rewards);
    pop.push_back(c);
  }

  for (int gen = 0; gen < max_gen; ++gen) {
    // Select new population
//    std::cout << "Calculating generation " << gen << std::endl;
    std::vector<Chromosome> new_pop;
    new_pop.reserve(pop_size);

    for (uint_fast32_t i = 0; i < pop_size; ++i) {
      new_pop.push_back(tournament_select(pop, tour_size, g));
    }

    // Cx
    for (int i = 0; i < 20; ++i) {
      std::vector<size_t> indices = get_population_sample(new_pop.size(), 2, g);
      std::pair<Chromosome, Chromosome> cx_ret = cx(new_pop[indices[0]], new_pop[indices[1]], cost_mat, max_cost, g);
      new_pop[indices[0]] = cx_ret.first;
      new_pop[indices[1]] = cx_ret.second;
      new_pop[indices[0]].evaluate_chromosome(cost_mat, rewards);
      new_pop[indices[1]].evaluate_chromosome(cost_mat, rewards);
    }

    // Mutate
    std::vector<size_t> indices = get_population_sample(new_pop.size(), 25, g);
    uint_fast32_t M = 8; //number of cores
    uint_fast32_t chunk_size = indices.size()/M;

    /*for (size_t idx : indices) {
      new_pop[idx] = mutate(new_pop[idx], cost_mat, rewards, max_cost);
    }*/
    std::vector< std::future<void> > future_v;
    for (uint_fast32_t thread_count = 0; thread_count < M; ++thread_count) {
      //std::launch::deferred|std::launch::async;
      std::vector<size_t > tmpv(indices.begin()+thread_count*chunk_size,indices.begin()+(thread_count+1)*chunk_size);
      future_v.push_back(std::async(std::launch::async, par_mutate, tmpv, std::ref(new_pop), std::ref(cost_mat), std::ref(rewards), std::ref(max_cost)));
    }
    if(indices.size()%M != 0){
      std::vector<size_t> tmpv(indices.begin()+(M)*chunk_size,indices.end());
      future_v.push_back(std::async(std::launch::async, par_mutate, tmpv, std::ref(new_pop), std::ref(cost_mat), std::ref(rewards), std::ref(max_cost)));
    }
    for(auto &f: future_v){
      f.get();
    }


    pop = new_pop;
  }
  std::sort(pop.begin(), pop.end(), [](Chromosome c1, Chromosome c2) { return c1.fitness > c2.fitness; });
  Chromosome best = pop[0];
  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices = best.path;
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));
  for (size_t vertex_idx = 1; vertex_idx < best.path.size() - 2; ++vertex_idx) {
    std::vector<uint_fast32_t> available_vertices;
    for (uint_fast32_t fv:free_vertices) {
      if (cost_mat[best.path[vertex_idx]][fv] < 2) { //TODO: Fix the constant sensor distance
        available_vertices.push_back(fv);
      }
    }

    double best_cost = best.cost;
    std::vector<uint_fast32_t> best_path = best.path;
    double best_fitness = best.fitness;
    uint_fast32_t best_vertex = best.path[vertex_idx];

    for (uint_fast32_t i:available_vertices) {
      Chromosome tmp_c;
      tmp_c.path = best.path;
      tmp_c.cost = best.cost
          - cost_mat[best.path[vertex_idx - 1]][best.path[vertex_idx]]
          - cost_mat[best.path[vertex_idx]][best.path[vertex_idx + 1]]
          + cost_mat[best.path[vertex_idx - 1]][i]
          + cost_mat[i][best.path[vertex_idx + 1]];
      tmp_c.path.erase(tmp_c.path.begin() + vertex_idx);
      tmp_c.path.insert(tmp_c.path.begin() + vertex_idx, i);
      if (tmp_c.cost <= max_cost) {
        tmp_c.evaluate_chromosome(cost_mat, rewards);
        if (tmp_c.fitness > best_fitness) {
          best_fitness = tmp_c.fitness;
          best_cost = tmp_c.cost;
          best_path = tmp_c.path;
          best_vertex = i;
        }
      }
    }
    best.path = best_path;
    best.fitness = best_fitness;
    best.cost = best_cost;
    if (best_vertex != best.path[vertex_idx]) {
      free_vertices.push_back(best.path[vertex_idx]);
      free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), best_vertex));
    }
  }
  return best;
}
