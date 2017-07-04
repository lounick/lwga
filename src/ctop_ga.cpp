//
// Created by nick on 30/03/16.
//

#include "ctop_ga.h"

void
Gene::evaluate_gene(Matrix<double> &cost_mat, std::vector<double> &rewards, std::vector<uint_fast32_t> &free_vertices) {
  fitness = 0;

//  static std::vector<uint_fast32_t> vertices;
//  if(vertices.size() == 0) {
//    vertices.reserve(cost_mat.size());
//    for (uint_fast32_t i = 0; i < cost_mat.size(); ++i)
//      vertices.push_back(i);
//  }

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
          extras += std::exp((log(0.01) / 2) * dist)
              * rewards[free_vertices[j]];; //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
        }
      }
      fitness += rewards[vertex] + extras;
    }
  }
  fitness = double(pow(fitness, 3)) / double(cost);
}

void Gene::calculate_cost(Matrix<double> &cost_mat) {
//  cost = 0;
//  for (Path::iterator it = path.begin() + 1; it != path.end(); ++it) {
//    cost += cost_mat[*(it - 1)][*it];
//    if (it != (path.end() - 1)) {
//      cost += 1;
//    }
//  }
  cost = get_path_cost(path, cost_mat);
}

void Gene::mutate(Matrix<double> &cost_mat, std::vector<double> &rewards, double max_cost, std::mt19937 &g) {

}

void Chromosome::evaluate_chromosome(Matrix<double> &cost_mat,
                                     std::vector<double> &rewards,
                                     std::vector<double> &max_cost_v) {

  uint_fast32_t num_robots = max_cost_v.size();
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    std::pair<Path, double> two_opt_return;
    two_opt_return = two_opt(genes[robot].path, cost_mat);
    genes[robot].path = two_opt_return.first;
    genes[robot].cost = two_opt_return.second + genes[robot].path.size() - 2;
  }

  std::unordered_set<uint_fast32_t> seen;
  seen.reserve(cost_mat.size());
  seen.insert(genes[0].path.front());
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
  for (std::vector<Gene>::iterator it = genes.begin(); it != genes.end(); ++it) {
    std::vector<uint_fast32_t>::iterator path_it = it->path.begin() + 1;
    while (path_it != it->path.end() - 1) {
      insert_ret = seen.insert(*path_it);
      if (!insert_ret.second) {
        // This is inefficient TODO:Optimise
        // See https://stackoverflow.com/questions/3747691/stdvector-iterator-invalidation
        path_it = it->path.erase(path_it);
      } else {
        ++path_it;
      }
    }
  }
  seen.insert(genes[0].path.back());

  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices(seen.begin(), seen.end());
//  visited_vertices.push_back(mutated.genes[0].path.front());
//  visited_vertices.push_back(mutated.genes[0].path.back());
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
//    mutated.genes[robot].calculate_cost(cost_mat);
    genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
  }

  total_fitness =
      std::accumulate(genes.begin(), genes.end(), 0.0, [](double total, const Gene &g) { return total + g.fitness; });
}

void Chromosome::mutate(Matrix<double> &cost_mat,
                        std::vector<double> &rewards,
                        std::vector<double> &max_cost_v,
                        std::mt19937 &g) {
  /*
   * Same as mutate for one vehicle but apply sequentially for many vehicles.
   */
  uint_fast32_t num_robots = max_cost_v.size();
  Chromosome mutated;
  mutated.genes.reserve(num_robots);

  //First perform a 2-opt in each gene.
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    Gene gene;
    gene.path.reserve(cost_mat.size());
    std::pair<Path, double> two_opt_return;
    two_opt_return = two_opt(genes[robot].path, cost_mat);
    gene.path = two_opt_return.first;
    gene.cost = two_opt_return.second + gene.path.size() - 2;
//    gene.evaluate_gene(cost_mat,rewards);
    mutated.genes.push_back(gene);

  }

  std::unordered_set<uint_fast32_t> seen1;
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret1;
  for (std::vector<Gene>::iterator it = mutated.genes.begin(); it != mutated.genes.end(); ++it) {
    for (std::vector<uint_fast32_t>::iterator path_it = it->path.begin() + 1; path_it != it->path.end() - 1;
         ++path_it) {
      insert_ret1 = seen1.insert(*path_it);
    }
  }

  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices(seen1.begin(), seen1.end());
  visited_vertices.push_back(mutated.genes[0].path.front());
  visited_vertices.push_back(mutated.genes[0].path.back());
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
//    mutated.genes[robot].calculate_cost(cost_mat);
    mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
  }

  // Remove any duplicates. Though there shouldn't be any.
  std::sort(mutated.genes.begin(),
            mutated.genes.end(),
            [](const Gene &g1, const Gene &g2) -> bool { return g1.fitness > g2.fitness; });
  std::unordered_set<uint_fast32_t> seen(mutated.genes.front().path.begin(), mutated.genes.front().path.end());
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
  for (std::vector<Gene>::iterator it = mutated.genes.begin() + 1; it != mutated.genes.end(); ++it) {
    std::vector<uint_fast32_t> new_path;
    new_path.reserve(it->path.size());
    new_path.push_back(it->path.front());
    for (std::vector<uint_fast32_t>::iterator path_it = it->path.begin() + 1; path_it != it->path.end() - 1;
         ++path_it) {
      insert_ret = seen.insert(*path_it);
      if (insert_ret.second) {
        new_path.push_back(*path_it);
      }
    }
    new_path.push_back(it->path.back());
    it->path = new_path;
  }

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    mutated.genes[robot].calculate_cost(cost_mat);
    mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
  }

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    for (uint_fast32_t iter = 0; iter < 10; ++iter) {
      if (std::generate_canonical<double, 10>(g) < 0.9) {
        if (mutated.genes[robot].cost >= 0.95 * max_cost_v[robot]) {
          //TODO: This is bound to have different effect in different grid sizes and budgets.
          //TODO: Should come up with something including the average travel cost or something smarter.
          /*
           * Pick random vertex from solution, apart from start and end.
           * Check it's free neighbours
           * Choose neighbour that maximises the fitness after removing vertex
           * If new fitness is better than old then accept solution
           */
          std::uniform_int_distribution<> dis(1, mutated.genes[robot].path.size() - 2);
          size_t rand_vertex_idx = dis(g);
          uint_fast32_t rand_vertex = mutated.genes[robot].path[rand_vertex_idx];
          std::vector<uint_fast32_t> available_vertices;
          for (const uint_fast32_t &i:free_vertices) {
            if (cost_mat[rand_vertex][i] < 2) { //TODO: fix the constant sensing distace
              available_vertices.push_back(i);
            }
          }

          double best_cost = mutated.genes[robot].cost;
          std::vector<uint_fast32_t> best_path = mutated.genes[robot].path;
          double best_fitness = mutated.genes[robot].fitness;
          uint_fast32_t best_vertex = rand_vertex;
          uint_fast32_t prev_vertex = mutated.genes[robot].path[rand_vertex_idx - 1];
          uint_fast32_t next_vertex = mutated.genes[robot].path[rand_vertex_idx + 1];
          double cost_removed = mutated.genes[robot].cost
              - cost_mat[prev_vertex][rand_vertex]
              - cost_mat[rand_vertex][next_vertex];
          for (uint_fast32_t i:available_vertices) {
            Gene tmp_c;
            tmp_c.path.reserve(cost_mat.size());
            tmp_c.path = mutated.genes[robot].path;
            tmp_c.cost = cost_removed
                + cost_mat[prev_vertex][i]
                + cost_mat[i][next_vertex];
            tmp_c.path.erase(tmp_c.path.begin() + rand_vertex_idx);
            tmp_c.path.insert(tmp_c.path.begin() + rand_vertex_idx, i);
            if (tmp_c.cost <= max_cost_v[robot]) {
              tmp_c.evaluate_gene(cost_mat, rewards, free_vertices);
              if (tmp_c.fitness > best_fitness) {
                best_fitness = tmp_c.fitness;
                best_cost = tmp_c.cost;
                best_path = tmp_c.path;
                best_vertex = i;
              }
            }
          }
          mutated.genes[robot].path = best_path;
          mutated.genes[robot].fitness = best_fitness;
          mutated.genes[robot].cost = best_cost;
          if (best_vertex != rand_vertex) {
            free_vertices.push_back(rand_vertex);
            free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), best_vertex));
          }
        } else {
          if (free_vertices.size() == 0)
            continue;
          //============================================================================================================
          /*if (mutated.genes[robot].path.size() == 2)
            continue;
          std::uniform_int_distribution<> dis(1, mutated.genes[robot].path.size() - 2);
          size_t vertex_idx = dis(g);
          uint_fast32_t vertex = mutated.genes[robot].path[vertex_idx];
          std::vector<uint_fast32_t>::iterator free_it;
          std::vector<uint_fast32_t> neighbors;
          neighbors.reserve(free_vertices.size());
          for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
            if (cost_mat[vertex][*free_it] <= 2) {
              neighbors.push_back(*free_it);
            }
          }


//           Iterate over neighbours and find each ones weight
//           Choose neighbour to add based on weight
//           Add him.

          if (neighbors.size() > 0) {
            std::uniform_real_distribution<> dis(0, 1);
            std::vector<double> distance_weights;
            std::vector<double> neighbour_weights;
            distance_weights.reserve(neighbors.size());
            neighbour_weights.reserve(neighbors.size());
            std::vector<double> combined;
            combined.reserve(neighbors.size());
            std::vector<double> cdf;
            cdf.reserve(neighbors.size());
            double tot_dist = 0;
            double tot_neighbours = 0;
            if (visited_vertices.size() > 0) {
              std::vector<uint_fast32_t>::iterator avail_it;
              for (avail_it = neighbors.begin(); avail_it != neighbors.end(); ++avail_it) {
                double mean_dist = 0;
                double free_neighbours = 0;
                std::vector<uint_fast32_t>::iterator used_it;
                for (used_it = visited_vertices.begin(); used_it != visited_vertices.end(); ++used_it) {
                  mean_dist += cost_mat[*avail_it][*used_it];
                }
                std::vector<uint_fast32_t>::iterator free_it;
                for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
                  if (cost_mat[*avail_it][*free_it] < 2) {
                    ++free_neighbours;
                  }
                }
                mean_dist = mean_dist / double(visited_vertices.size());
                tot_dist += mean_dist;
                tot_neighbours += free_neighbours;
                distance_weights.push_back(mean_dist);
                neighbour_weights.push_back(free_neighbours);
              }

              std::vector<double>::iterator dist_weight_it = distance_weights.begin();
              std::vector<double>::iterator neigh_weight_it = neighbour_weights.begin();
              double tot_combined = 0;
              for (; dist_weight_it != distance_weights.end(); ++dist_weight_it, ++neigh_weight_it) {
                combined.push_back(*dist_weight_it + *neigh_weight_it);
                tot_combined += *dist_weight_it + *neigh_weight_it;
              }
              std::vector<double>::iterator combined_it;
              cdf.push_back(combined.front() / tot_combined);
              for (combined_it = combined.begin() + 1; combined_it != combined.end(); ++combined_it) {
                cdf.push_back(cdf.back() + (*combined_it / tot_combined));
              }
            }
            double prob = dis(g);
            auto vertex_pos = std::lower_bound(cdf.cbegin(), cdf.cend(), prob);
            uint_fast32_t next_vertex = neighbors[std::distance(cdf.cbegin(), vertex_pos)];
//            double next_vertex_cost = cost_mat[path.back()][next_vertex] + 1;
            double travel_increase = cost_mat[vertex][next_vertex]
                + cost_mat[next_vertex][mutated.genes[robot].path[vertex_idx + 1]]
                - cost_mat[vertex][mutated.genes[robot].path[vertex_idx + 1]];
            if (mutated.genes[robot].cost + travel_increase + 1 <= max_cost_v[robot]) {
              mutated.genes[robot].path.insert(mutated.genes[robot].path.begin() + vertex_idx + 1, next_vertex);
              mutated.genes[robot].cost += travel_increase + 1;
              mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
              free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), next_vertex));
              visited_vertices.push_back(next_vertex);
            }
          }*/
          //============================================================================================================
          std::uniform_int_distribution<> dis(0, free_vertices.size() - 1);
          size_t vertex_idx = dis(g);
          uint_fast32_t vertex = free_vertices[vertex_idx];

          double best_travel_increase = 0;
          double best_fitness = 0;
          size_t ins_pos = 0;

          for (size_t i = 1; i < mutated.genes[robot].path.size(); ++i) {
            double travel_increase = cost_mat[mutated.genes[robot].path[i - 1]][vertex]
                + cost_mat[vertex][mutated.genes[robot].path[i]]
                - cost_mat[mutated.genes[robot].path[i - 1]][mutated.genes[robot].path[i]];
            double fit = rewards[vertex]; //TODO: Should we add the extras to the fit?
            if (travel_increase != 0) {
              fit /= travel_increase;
            } else {
              fit = std::numeric_limits<double>::infinity();
            }
            if (fit > best_fitness) {
              if (mutated.genes[robot].cost + travel_increase + 1 <= max_cost_v[robot]) {
                best_travel_increase = travel_increase;
                best_fitness = fit;
                ins_pos = i;
              }
            }
          }
          if (ins_pos > 0) {
            mutated.genes[robot].path.insert(mutated.genes[robot].path.begin() + ins_pos, vertex);
            mutated.genes[robot].cost += best_travel_increase + 1;
            mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
            free_vertices.erase(free_vertices.begin() + vertex_idx);
          }
          //============================================================================================================
        }
      } else {
        std::unordered_set<uint_fast32_t> duplicates;
        std::unordered_set<uint_fast32_t> seen_vertices;
        std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
        for (uint_fast32_t i:mutated.genes[robot].path) {
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
          std::vector<uint_fast32_t>::iterator it = mutated.genes[robot].path.begin();
          while ((it = std::find(it, mutated.genes[robot].path.end(), duplicate)) != mutated.genes[robot].path.end()) {
            indices.push_back(std::distance(mutated.genes[robot].path.begin(), it));
            ++it;
          }

          // Remove the one with the minimum loss
          double min_loss = std::numeric_limits<double>::infinity();
          size_t to_remove = 0;

          for (size_t idx:indices) {
            double travel_decrease;
            if (idx < mutated.genes[robot].path.size() - 1)
              travel_decrease = cost_mat[mutated.genes[robot].path[idx - 1]][mutated.genes[robot].path[idx]]
                  + cost_mat[mutated.genes[robot].path[idx]][mutated.genes[robot].path[idx + 1]]
                  - cost_mat[mutated.genes[robot].path[idx - 1]][mutated.genes[robot].path[idx + 1]];
            else
              travel_decrease = 0;

            double loss = rewards[mutated.genes[robot].path[idx]];

            double extras = 0;
            for (size_t j = 0; j < free_vertices.size(); ++j) {
              if (cost_mat[mutated.genes[robot].path[idx]][free_vertices[j]] < 2) {
                extras += std::exp(-2 * cost_mat[mutated.genes[robot].path[idx]][free_vertices[j]]);
              }
            }
            loss += extras;

            if (travel_decrease != 0) {
              loss /= travel_decrease;
            } else {
              loss = std::numeric_limits<double>::infinity();
            }

            if (loss <= min_loss) {
              min_loss = loss;
              to_remove = idx;
            }
          }

          if (to_remove > 0) {
            mutated.genes[robot].path.erase(mutated.genes[robot].path.begin() + to_remove);
            mutated.genes[robot].calculate_cost(cost_mat);
            mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
          }
        } else {
          if (mutated.genes[robot].cost >= 0.95 * max_cost_v[robot]) {
            size_t to_remove = 0;
            double min_loss = std::numeric_limits<double>::infinity();
            for (size_t i = 1; i < mutated.genes[robot].path.size() - 1; ++i) {
              double travel_decrease =
                  cost_mat[mutated.genes[robot].path[i - 1]][mutated.genes[robot].path[i]]
                      + cost_mat[mutated.genes[robot].path[i]][mutated.genes[robot].path[i + 1]]
                      - cost_mat[mutated.genes[robot].path[i - 1]][mutated.genes[robot].path[i + 1]];
              double loss = rewards[mutated.genes[robot].path[i]];
              double extras = 0;
              for (size_t j = 0; j < free_vertices.size(); ++j) {
                if (cost_mat[mutated.genes[robot].path[i]][free_vertices[j]] < 2) {
                  extras += std::exp(-2 * cost_mat[mutated.genes[robot].path[i]][free_vertices[j]]);
                }
              }
              loss += extras;

              if (travel_decrease != 0) {
                loss /= travel_decrease;
              } else {
                loss = std::numeric_limits<double>::infinity();
              }

              if (loss <= min_loss) {
                min_loss = loss;
                to_remove = i;
              }
            }
            if (to_remove > 0) {
              free_vertices.push_back(mutated.genes[robot].path[to_remove]);
              mutated.genes[robot].path.erase(mutated.genes[robot].path.begin() + to_remove);
              mutated.genes[robot].calculate_cost(cost_mat);
              mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
            }
          }
        }
      }
    }
  }
  genes = mutated.genes;
  evaluate_chromosome(cost_mat, rewards, max_cost_v);
}

Chromosome generate_chromosome(Matrix<double> &cost_mat,
                               std::vector<double> &max_cost_v,
                               uint idx_start,
                               uint idx_finish,
                               std::mt19937 &g) {
  uint n_agents = max_cost_v.size();

  Chromosome c;
  c.genes.reserve(n_agents);

  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);

//  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_start));
//  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_finish));

  std::vector<uint_fast32_t> free_vertices = vertices;
  free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), idx_start));
  free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), idx_finish));
  std::vector<uint_fast32_t> used_vertices;
  used_vertices.reserve(cost_mat.size());
  for (size_t i = 0; i < n_agents; i++) {
    Path path;
    path.reserve(cost_mat.size());

    bool done = false;
    double total_cost = 0;
    path.push_back(idx_start);
    std::vector<uint_fast32_t> checked;
    checked.reserve(cost_mat.size());
    checked.clear();
    while (!done) {
      std::vector<uint_fast32_t> neighbours;
      neighbours.reserve(cost_mat.size());
      neighbours.clear();
      neighbours.push_back(path.back());
      std::vector<uint_fast32_t> available;
      available.reserve(cost_mat.size());
      available.clear();
      bool neighbours_found = false;
      while (!neighbours_found) {
        expand_neighbours(neighbours, checked, vertices, cost_mat);
        if (neighbours.size() == 0) { neighbours_found = true; }
        else {
          std::sort(neighbours.begin(), neighbours.end());
          std::set_intersection(neighbours.begin(),
                                neighbours.end(),
                                free_vertices.begin(),
                                free_vertices.end(),
                                std::back_inserter(available));
          if (available.size() > 0) {
            neighbours_found = true;
          }
        }
      }
      if (available.size() > 0) {
        //==============================================================================================================
        //choose randomly
//        size_t rand_idx = g() % available.size();
//        uint_fast32_t next_vertex = available[rand_idx];
        //==============================================================================================================
        //chose randomly based on weight
        std::uniform_real_distribution<> dis(0, 1);
        std::vector<double> distance_weights;
        std::vector<double> neighbour_weights;
        distance_weights.reserve(available.size());
        neighbour_weights.reserve(available.size());
        std::vector<double> combined;
        combined.reserve(available.size());
        std::vector<double> cdf;
        cdf.reserve(available.size());
        double tot_dist = 0;
        double tot_neighbours = 0;
        if (used_vertices.size() > 0) {
          std::vector<uint_fast32_t>::iterator avail_it;
          for (avail_it = available.begin(); avail_it != available.end(); ++avail_it) {
            double mean_dist = 0;
            double free_neighbours = 0;
            std::vector<uint_fast32_t>::iterator used_it;
            for (used_it = used_vertices.begin(); used_it != used_vertices.end(); ++used_it) {
              mean_dist += cost_mat[*avail_it][*used_it];
            }
            std::vector<uint_fast32_t>::iterator free_it;
            for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
              if (cost_mat[*avail_it][*free_it] < 2) {
                ++free_neighbours;
              }
            }
            mean_dist = mean_dist / double(used_vertices.size());
            tot_dist += mean_dist;
            tot_neighbours += free_neighbours;
            distance_weights.push_back(mean_dist);
            neighbour_weights.push_back(free_neighbours);
          }

          std::vector<double>::iterator dist_weight_it = distance_weights.begin();
          std::vector<double>::iterator neigh_weight_it = neighbour_weights.begin();
          double tot_combined = 0;
          for (; dist_weight_it != distance_weights.end(); ++dist_weight_it, ++neigh_weight_it) {
            combined.push_back(*dist_weight_it + *neigh_weight_it);
            tot_combined += *dist_weight_it + *neigh_weight_it;
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
        //==============================================================================================================
        double next_vertex_cost = cost_mat[path.back()][next_vertex] + 1;
        if (total_cost + next_vertex_cost + cost_mat[next_vertex][idx_finish] <= max_cost_v[i]) {
          total_cost += next_vertex_cost;
          path.push_back(next_vertex);
          free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), next_vertex));
          used_vertices.push_back(next_vertex);
        } else {
          done = true;
        }
      } else {
        done = true;
      }
    }
    path.push_back(idx_finish);
    Gene gene;
    gene.path = path;
    c.genes.push_back(gene);
  }
  return c;
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
      if (dist <= 2) {
        nset.insert(*free_it);
      }
    }
  }
  neighbours.assign(nset.begin(), nset.end());
}


//Chromosome generate_chromosome(Matrix<double> &cost_mat, std::vector<double> &max_cost_v, uint idx_start, uint idx_finish, std::mt19937 &g) {
//
//  uint n_agents = max_cost_v.size();
//
//  Chromosome c;
//  c.genes.reserve(n_agents);
//
//  std::vector<uint> vertices(cost_mat.size());
//  std::iota(vertices.begin(), vertices.end(), 0);
//
//  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_start));
//  vertices.erase(std::remove(vertices.begin(), vertices.end(), idx_finish));
//
//  for (size_t i = 0; i < n_agents; i++) {
//    Path path;
//    path.reserve(cost_mat.size());
//
//    bool done = false;
//    double total_cost = 0;
//    path.push_back(idx_start);
//    while (!done) {
//      size_t vsize = vertices.size();
//      if (vsize == 0)
//        break;
//      size_t rand_idx = g() % vsize;
//      uint next_vertex = vertices[rand_idx];
//      double next_vertex_cost = cost_mat[path.back()][next_vertex] + 1;
//      if (total_cost + next_vertex_cost + cost_mat[next_vertex][idx_finish] <= max_cost_v[i]) {
//        total_cost += next_vertex_cost;
//        path.push_back(next_vertex);
//        vertices.erase(std::remove(vertices.begin(), vertices.end(), next_vertex));
//      }
//      else {
//        done = true;
//      }
//    }
//    path.push_back(idx_finish);
//    Gene gene;
//    gene.path = path;
//    c.genes.push_back(gene);
//  }
//  return c;
//}

Chromosome tournament_select(std::vector<Chromosome> &population, uint tour_size, std::mt19937 &g) {

  std::vector<uint> indices(population.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), g);
  double max_fitness = DBL_MIN;
  Chromosome best_chromosome;

  for (size_t i = 0; i < tour_size; i++) {
    if (population[indices[i]].total_fitness > max_fitness) {
      max_fitness = population[indices[i]].total_fitness;
      best_chromosome = population[indices[i]];
    }
  }

  return best_chromosome;
}

void cx(Chromosome &c1,
        Chromosome &c2,
        std::vector<std::vector<double> > &cost_mat,
        std::vector<double> &max_cost_v,
        std::vector<double> &rewards) {
  /*
   * From each parent select the most fit path and exchange with the lowest fit path of the other parent.
   * Then check for double visits etc.
   * Then evaluate fitnesses and costs again.
   */
  uint_fast32_t num_robots = max_cost_v.size();

  std::sort(c1.genes.begin(),
            c1.genes.end(),
            [](const Gene &g1, const Gene &g2) -> bool { return g1.fitness > g2.fitness; });
  std::sort(c2.genes.begin(),
            c2.genes.end(),
            [](const Gene &g1, const Gene &g2) -> bool { return g1.fitness > g2.fitness; });
  //Gene g1_best = c1.genes.front();
  c1.genes.back() = c2.genes.front();
  c2.genes.back() = c1.genes.front();
  std::sort(c1.genes.begin(),
            c1.genes.end(),
            [](const Gene &g1, const Gene &g2) -> bool { return g1.fitness > g2.fitness; });
  std::sort(c2.genes.begin(),
            c2.genes.end(),
            [](const Gene &g1, const Gene &g2) -> bool { return g1.fitness > g2.fitness; });

  std::unordered_set<uint_fast32_t> seen(c1.genes.front().path.begin(), c1.genes.front().path.end());
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
  for (std::vector<Gene>::iterator it = c1.genes.begin() + 1; it != c1.genes.end(); ++it) {
    std::vector<uint_fast32_t> new_path;
    new_path.reserve(it->path.size());
    new_path.push_back(it->path.front());
    for (std::vector<uint_fast32_t>::iterator path_it = it->path.begin() + 1; path_it != it->path.end() - 1;
         ++path_it) {
      insert_ret = seen.insert(*path_it);
      if (insert_ret.second) {
        new_path.push_back(*path_it);
      }
    }
    new_path.push_back(it->path.back());
    it->path = new_path;
  }

  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::vector<uint_fast32_t> visited_vertices(seen.begin(), seen.end());
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    c1.genes[robot].calculate_cost(cost_mat);
    c1.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
  }
  c1.evaluate_chromosome(cost_mat, rewards, max_cost_v);

  seen.clear();
  for (std::vector<uint_fast32_t>::iterator path_it = c2.genes.front().path.begin();
       path_it != c2.genes.front().path.end(); ++path_it) {
    seen.insert(*path_it);
  }
  for (std::vector<Gene>::iterator it = c2.genes.begin() + 1; it != c2.genes.end(); ++it) {
    std::vector<uint_fast32_t> new_path;
    new_path.reserve(it->path.size());
    new_path.push_back(it->path.front());
    for (std::vector<uint_fast32_t>::iterator path_it = it->path.begin() + 1; path_it != it->path.end() - 1;
         ++path_it) {
      insert_ret = seen.insert(*path_it);
      if (insert_ret.second) {
        new_path.push_back(*path_it);
      }
    }
    new_path.push_back(it->path.back());
    it->path = new_path;
  }
  visited_vertices.clear();
  visited_vertices = std::vector<uint_fast32_t>(seen.begin(), seen.end());
  std::sort(visited_vertices.begin(), visited_vertices.end());
  free_vertices.clear();
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    c2.genes[robot].calculate_cost(cost_mat);
    c2.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
  }
  c2.evaluate_chromosome(cost_mat, rewards, max_cost_v);
}

void par_mutate(std::vector<size_t> indices,
                std::vector<Chromosome> &pop,
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                std::vector<double> &max_cost_v) {
  // Get hash of thread id for the seed of the generator.
  std::hash<std::thread::id> hasher;
  static thread_local std::mt19937 g(hasher(std::this_thread::get_id()));
  for (size_t idx:indices)
    pop[idx].mutate(cost_mat, rewards, max_cost_v, g);
}

Chromosome ga_ctop(Matrix<double> &cost_mat,
                   std::vector<double> &rewards,
                   std::vector<double> max_cost_v,
                   uint idx_start,
                   uint idx_finish,
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
  int tour_size = 5;
  int max_gen = 50;
  uint_fast32_t num_robots = max_cost_v.size();

  // Initialise population
  std::vector<Chromosome> pop;
  pop.reserve(pop_size);

  for (uint_fast32_t i = 0; i < pop_size; ++i) {
    Chromosome c = generate_chromosome(cost_mat, max_cost_v, idx_start, idx_finish, g);
    std::vector<uint_fast32_t> vertices(cost_mat.size());
    std::iota(vertices.begin(), vertices.end(), 0);
    std::vector<uint_fast32_t> free_vertices;
    free_vertices.reserve(cost_mat.size());
    std::unordered_set<uint_fast32_t> seen;
    seen.reserve(cost_mat.size());
    std::vector<uint_fast32_t> visited_vertices;
    visited_vertices.reserve(cost_mat.size());
//    for(int i = 0; i < num_robots-1; ++i){
//      for (int j = i+1; j < num_robots; ++j){
//        std::vector<uint_fast32_t> mutual_opt;
//        mutual_opt.reserve(cost_mat.size()*2);
//        mutual_opt.assign(c.genes[i].path.begin(),c.genes[i].path.end());
//        for (int k = c.genes[j].path.size()-1; k >= 0; --k)
//          mutual_opt.push_back(c.genes[j].path[k]);
//        std::pair<std::vector<uint_fast32_t>, double> mutual_opt_ret = two_opt(mutual_opt, cost_mat);
//        std::vector<uint_fast32_t>::iterator it = mutual_opt_ret.first.begin();
//        c.genes[i].path.clear();
//        while (*it != idx_finish){
//          c.genes[i].path.push_back(*it);
//          ++it;
//        }
//        c.genes[i].path.push_back(idx_finish);
//        c.genes[j].path.clear();
//        it = mutual_opt_ret.first.end()-1;
//        while (*it != idx_finish){
//          c.genes[j].path.push_back(*it);
//          --it;
//        }
//        c.genes[j].path.push_back(idx_finish);
//      }
//    }
    for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
      std::pair<std::vector<uint_fast32_t>, double> two_opt_ret = two_opt(c.genes[robot].path, cost_mat);
      c.genes[robot].path = two_opt_ret.first;
      c.genes[robot].cost = two_opt_ret.second;

      for (size_t path_idx = 0; path_idx < c.genes[robot].path.size(); ++path_idx) {
        seen.insert(c.genes[robot].path[path_idx]);
      }
      free_vertices.clear();
      visited_vertices.clear();
      visited_vertices = std::vector<uint_fast32_t>(seen.begin(), seen.end());
      std::sort(visited_vertices.begin(), visited_vertices.end());
      std::set_difference(vertices.begin(),
                          vertices.end(),
                          visited_vertices.begin(),
                          visited_vertices.end(),
                          std::back_inserter(free_vertices));
      c.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
    }
    c.evaluate_chromosome(cost_mat, rewards, max_cost_v);
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
//      cx(new_pop[indices[0]], new_pop[indices[1]], cost_mat, max_cost_v, rewards);
    }

    // Mutate
    std::vector<size_t> indices = get_population_sample(new_pop.size(), 25, g);
    uint_fast32_t M = 1; //number of cores
    uint_fast32_t chunk_size = indices.size() / M;

    /*for (size_t idx : indices) {
      new_pop[idx] = mutate(new_pop[idx], cost_mat, rewards, max_cost);
    }*/
    std::vector<std::future<void> > future_v;
    for (uint_fast32_t thread_count = 0; thread_count < M; ++thread_count) {
      //std::launch::deferred|std::launch::async;
      std::vector<size_t>
          tmpv(indices.begin() + thread_count * chunk_size, indices.begin() + (thread_count + 1) * chunk_size);
      future_v.push_back(std::async(std::launch::async,
                                    par_mutate,
                                    tmpv,
                                    std::ref(new_pop),
                                    std::ref(cost_mat),
                                    std::ref(rewards),
                                    std::ref(max_cost_v)));
    }
    if (indices.size() % M != 0) {
      std::vector<size_t> tmpv(indices.begin() + (M) * chunk_size, indices.end());
      future_v.push_back(std::async(std::launch::async,
                                    par_mutate,
                                    tmpv,
                                    std::ref(new_pop),
                                    std::ref(cost_mat),
                                    std::ref(rewards),
                                    std::ref(max_cost_v)));
    }
    for (auto &f: future_v) {
      f.get();
    }


    pop = new_pop;
  }
  std::sort(pop.begin(), pop.end(), [](Chromosome c1, Chromosome c2) { return c1.total_fitness > c2.total_fitness; });
  Chromosome best = pop[0];
  std::vector<uint_fast32_t> vertices(cost_mat.size());
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<uint_fast32_t> free_vertices;
  std::unordered_set<uint_fast32_t> seen(best.genes[0].path.begin(), best.genes[0].path.end());
  for (uint_fast32_t robot = 1; robot < num_robots; ++robot) {
    for (uint_fast32_t vertex = 0; vertex < best.genes[robot].path.size(); ++vertex) {
      seen.insert(vertex);
    }
  }
  std::vector<uint_fast32_t> visited_vertices(seen.begin(), seen.end());
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(vertices.begin(),
                      vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    for (size_t vertex_idx = 1; vertex_idx < best.genes[robot].path.size() - 2; ++vertex_idx) {
      std::vector<uint_fast32_t> available_vertices;
      for (uint_fast32_t fv:free_vertices) {
        if (cost_mat[best.genes[robot].path[vertex_idx]][fv] < 2) { //TODO: Fix the constant sensor distance
          available_vertices.push_back(fv);
        }
      }

      double best_cost = best.genes[robot].cost;
      std::vector<uint_fast32_t> best_path = best.genes[robot].path;
      double best_fitness = best.genes[robot].fitness;
      uint_fast32_t best_vertex = best.genes[robot].path[vertex_idx];

      for (uint_fast32_t i:available_vertices) {
        Gene tmp_c;
        tmp_c.path = best.genes[robot].path;
        tmp_c.cost = best.genes[robot].cost
            - cost_mat[best.genes[robot].path[vertex_idx - 1]][best.genes[robot].path[vertex_idx]]
            - cost_mat[best.genes[robot].path[vertex_idx]][best.genes[robot].path[vertex_idx + 1]]
            + cost_mat[best.genes[robot].path[vertex_idx - 1]][i]
            + cost_mat[i][best.genes[robot].path[vertex_idx + 1]];
        tmp_c.path.erase(tmp_c.path.begin() + vertex_idx);
        tmp_c.path.insert(tmp_c.path.begin() + vertex_idx, i);
        if (tmp_c.cost <= max_cost_v[robot]) {
          tmp_c.evaluate_gene(cost_mat, rewards, free_vertices);
          if (tmp_c.fitness > best_fitness) {
            best_fitness = tmp_c.fitness;
            best_cost = tmp_c.cost;
            best_path = tmp_c.path;
            best_vertex = i;
          }
        }
      }
      best.genes[robot].path = best_path;
      best.genes[robot].fitness = best_fitness;
      best.genes[robot].cost = best_cost;
      if (best_vertex != best.genes[robot].path[vertex_idx]) {
        free_vertices.push_back(best.genes[robot].path[vertex_idx]);
        free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), best_vertex));
      }
    }
  }
  return best;
}
