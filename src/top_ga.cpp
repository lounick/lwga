//
// Created by nick on 21/07/17.
//

#include <include/top_ga.h>

void Gene::calculate_cost(Matrix<double> &cost_mat) {
  cost = get_path_cost(path, cost_mat);
}

void Gene::evaluate_gene(std::vector<double> &rewards) {
  fitness = 0;
  for(int i = 0; i < path.size(); ++i){
    fitness += rewards[path[i]];
  }
  fitness = double_t(pow(fitness,3))/double_t(cost);
}

void Chromosome::evaluate_chromosome(Matrix<double> &cost_mat,
                                     std::vector<double> &rewards,
                                     std::vector<double> &max_cost_v) {

  uint_fast32_t num_robots = max_cost_v.size();
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    std::pair<Path, double> two_opt_return;
    two_opt_return = two_opt(genes[robot].path, cost_mat);
    genes[robot].path = two_opt_return.first;
    genes[robot].cost = two_opt_return.second;
  }

  seen_vertices.clear();
  free_vertices.clear();

  seen_vertices.insert(genes[0].path.front());
  std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
  for (std::vector<Gene>::iterator it = genes.begin(); it != genes.end(); ++it) {
    std::vector<uint_fast32_t>::iterator path_it = it->path.begin() + 1;
    while (path_it != it->path.end() - 1) {
      insert_ret = seen_vertices.insert(*path_it);
      if (!insert_ret.second) {
        // This is inefficient TODO:Optimise
        // See https://stackoverflow.com/questions/3747691/stdvector-iterator-invalidation
        path_it = it->path.erase(path_it);
      } else {
        ++path_it;
      }
    }
  }
  seen_vertices.insert(genes[0].path.back());


  std::vector<uint_fast32_t> visited_vertices(seen_vertices.begin(), seen_vertices.end());
  std::sort(visited_vertices.begin(), visited_vertices.end());
  std::set_difference(all_vertices.begin(),
                      all_vertices.end(),
                      visited_vertices.begin(),
                      visited_vertices.end(),
                      std::back_inserter(free_vertices));

  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    genes[robot].calculate_cost(cost_mat);
    genes[robot].evaluate_gene(rewards);
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
  if(genes.size() == 0 || genes.front().path.size() == 0)
    return;
  uint_fast32_t start_vertex = genes.front().path.front();
  uint_fast32_t end_vertex = genes.front().path.back();
  Chromosome mutated(cost_mat.size(), num_robots, start_vertex, end_vertex);
  mutated.genes.reserve(num_robots);
  mutated.genes = genes;

  std::uniform_real_distribution<> dis(0, 1);
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    for (uint_fast32_t iter = 0; iter < 10; ++iter) {
      if (dis(g) < 0.5) {
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
              tmp_c.evaluate_gene(rewards);
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
            seen_vertices.insert(best_vertex);
            seen_vertices.erase(rand_vertex);
          }
        } else {
          if (free_vertices.size() == 0)
            continue;
          //============================================================================================================
//          if (mutated.genes[robot].path.size() == 2)
//            continue;
//          std::uniform_int_distribution<> dis(1, mutated.genes[robot].path.size() - 2);
//          size_t vertex_idx = dis(g);
//          uint_fast32_t vertex = mutated.genes[robot].path[vertex_idx];
//          std::vector<uint_fast32_t>::iterator free_it;
//          std::vector<uint_fast32_t> neighbors;
//          neighbors.reserve(free_vertices.size());
//          for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
//            if (cost_mat[vertex][*free_it] < 2) {
//              neighbors.push_back(*free_it);
//            }
//          }
//
//
////           Iterate over neighbours and find each ones weight
////           Choose neighbour to add based on weight
////           Add him.
//
//          if (neighbors.size() > 0) {
//            std::uniform_real_distribution<> dis(0, 1);
//            std::vector<double> distance_weights;
//            std::vector<double> neighbour_weights;
//            distance_weights.reserve(neighbors.size());
//            neighbour_weights.reserve(neighbors.size());
//            std::vector<double> combined;
//            combined.reserve(neighbors.size());
//            std::vector<double> cdf;
//            cdf.reserve(neighbors.size());
//            double tot_dist = 0;
//            double tot_neighbours = 0;
//            if (seen_vertices.size() > 0) {
//              std::vector<uint_fast32_t>::iterator avail_it;
//              for (avail_it = neighbors.begin(); avail_it != neighbors.end(); ++avail_it) {
//                double mean_dist = 0;
//                double free_neighbours = 0;
//                std::unordered_set<uint_fast32_t>::iterator used_it;
//                for (used_it = seen_vertices.begin(); used_it != seen_vertices.end(); ++used_it) {
//                  mean_dist += cost_mat[*avail_it][*used_it];
//                }
//                std::vector<uint_fast32_t>::iterator free_it;
//                for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
//                  if (cost_mat[*avail_it][*free_it] < 2) {
//                    ++free_neighbours;
//                  }
//                }
//                mean_dist = mean_dist / double(seen_vertices.size());
//                tot_dist += mean_dist;
//                tot_neighbours += free_neighbours;
//                distance_weights.push_back(mean_dist);
//                neighbour_weights.push_back(free_neighbours);
//              }
//
//              std::vector<double>::iterator dist_weight_it = distance_weights.begin();
//              std::vector<double>::iterator neigh_weight_it = neighbour_weights.begin();
//              double tot_combined = 0;
//              for (; dist_weight_it != distance_weights.end(); ++dist_weight_it, ++neigh_weight_it) {
//                combined.push_back(*dist_weight_it + *neigh_weight_it);
//                tot_combined += *dist_weight_it + *neigh_weight_it;
//              }
//              std::vector<double>::iterator combined_it;
//              cdf.push_back(combined.front() / tot_combined);
//              for (combined_it = combined.begin() + 1; combined_it != combined.end(); ++combined_it) {
//                cdf.push_back(cdf.back() + (*combined_it / tot_combined));
//              }
//            }
//            double prob = dis(g);
//            auto vertex_pos = std::lower_bound(cdf.cbegin(), cdf.cend(), prob);
//            uint_fast32_t next_vertex = neighbors[std::distance(cdf.cbegin(), vertex_pos)];
////            double next_vertex_cost = cost_mat[path.back()][next_vertex] + 1;
//            double travel_increase = cost_mat[vertex][next_vertex]
//                + cost_mat[next_vertex][mutated.genes[robot].path[vertex_idx + 1]]
//                - cost_mat[vertex][mutated.genes[robot].path[vertex_idx + 1]];
//            if (mutated.genes[robot].cost + travel_increase + 1 <= max_cost_v[robot]) {
//              mutated.genes[robot].path.insert(mutated.genes[robot].path.begin() + vertex_idx + 1, next_vertex);
//              mutated.genes[robot].cost += travel_increase + 1;
//              mutated.genes[robot].evaluate_gene(cost_mat, rewards, free_vertices);
//              free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), next_vertex));
//              seen_vertices.insert(next_vertex);
//            }
//          }
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
            mutated.genes[robot].evaluate_gene(rewards);
            free_vertices.erase(free_vertices.begin() + vertex_idx);
          }
          //============================================================================================================
        }
      } else {
        if(mutated.genes[robot].path.size() < 3){
          continue;
        }
//        if (mutated.genes[robot].cost >= 0.95 * max_cost_v[robot]) {
        size_t to_remove = 0;
        double min_loss = std::numeric_limits<double>::infinity();
        for (size_t i = 1; i < mutated.genes[robot].path.size() - 1; ++i) {
          double travel_decrease =
              cost_mat[mutated.genes[robot].path[i - 1]][mutated.genes[robot].path[i]]
                  + cost_mat[mutated.genes[robot].path[i]][mutated.genes[robot].path[i + 1]]
                  - cost_mat[mutated.genes[robot].path[i - 1]][mutated.genes[robot].path[i + 1]];
          double loss = rewards[mutated.genes[robot].path[i]];

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
          mutated.genes[robot].evaluate_gene(rewards);
        }
//        }
      }
    }
  }
  genes = mutated.genes;
  evaluate_chromosome(cost_mat, rewards, max_cost_v);
}

Chromosome::Chromosome(size_t num_vertices, size_t num_genes, const uint_fast32_t &start_vertex,
                       const uint_fast32_t &end_vertex) {
  std::cout << num_vertices << " " << num_genes << " " << start_vertex << " " << end_vertex << std::endl;
  all_vertices = std::vector<uint_fast32_t>(num_vertices);
  std::iota(all_vertices.begin(), all_vertices.end(), 0);
  seen_vertices.reserve(num_vertices);
  free_vertices = std::vector<uint_fast32_t>(all_vertices.begin(), all_vertices.end());
  free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), start_vertex));
  free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), end_vertex));
  genes.reserve(num_genes);
  for (size_t i = 0; i < num_genes; ++i) {
    genes.push_back(Gene());
    genes[i].cost = 0;
    genes[i].fitness = 0;
    genes[i].path.reserve(num_vertices);
  }
}

Chromosome::Chromosome() {
  ;
}

void Chromosome::GenerateInsertMoves(const Matrix<double_t> &cost_mat,
                                     const std::vector<double_t> &rewards,
                                     const double_t &max_cost,
                                     const Path &path, double_t path_cost,
                                     std::vector<InsertMove> &moves,
                                     double_t &min, double_t &max) {

  moves.clear();
  min = std::numeric_limits<double_t>::max();
  max = std::numeric_limits<double_t>::min();

  for (const uint_fast32_t &v : free_vertices) {
    //Calculate the reward of the free vertex
    double_t reward = 0;

    reward = rewards[v];

    //Generate all the possible insert moves of this vertex to this path
    Path::const_iterator p;
    for (p = path.begin() + 1; p != path.end(); ++p) {
      uint_fast32_t prev_vertex = *(p - 1);
      uint_fast32_t next_vertex = *p;
      double_t cost_increase = cost_mat[prev_vertex][v] + cost_mat[v][next_vertex] - cost_mat[prev_vertex][next_vertex] + 1;
      if ((path_cost + cost_increase) < max_cost) {
        double_t heuristic = reward / cost_increase;
        if (heuristic < min) {
          min = heuristic;
        }
        if (heuristic > max) {
          max = heuristic;
        }
        moves.emplace_back(v, prev_vertex, next_vertex, cost_increase, reward, heuristic);
      }
    }
  }
}

double_t Chromosome::GenerateGRASPPath(const Matrix<double_t> &cost_mat,
                                       const std::vector<double_t> &rewards,
                                       const double_t &max_cost,
                                       const uint_fast32_t &start_vertex,
                                       const uint_fast32_t &end_vertex,
                                       std::mt19937 &g, Path &path) {
  std::uniform_real_distribution<> dis(0, 1);
  double_t greediness = dis(g);
  path.push_back(start_vertex);
  path.push_back(end_vertex);
  double_t path_cost = cost_mat[start_vertex][end_vertex];
  std::vector<InsertMove> insert_candidates;
  insert_candidates.reserve(free_vertices.size());
  double_t min, max;
  GenerateInsertMoves(cost_mat, rewards, max_cost, path, path_cost, insert_candidates, min, max);
  while (insert_candidates.size() > 0) {
    double_t threshold = min + greediness * (max - min);
    std::vector<uint_fast32_t> restricted_candidates_indexes;
    restricted_candidates_indexes.reserve(insert_candidates.size());

    for (size_t i = 0; i < insert_candidates.size(); ++i) {
      if (insert_candidates[i].heuristic >= threshold) {
        restricted_candidates_indexes.push_back(i);
      }
    }
//    std::cout << "insert candidates: " << insert_candidates.size() << " RIC: " << restricted_candidates_indexes.size() << std::endl;
    if (restricted_candidates_indexes.size() > 0) {
      std::uniform_int_distribution<size_t> index_dis(0, restricted_candidates_indexes.size() - 1);
      size_t idx = index_dis(g);
      path.insert(std::find(path.begin(), path.end(), insert_candidates[restricted_candidates_indexes[idx]].next_vertex),
                  insert_candidates[restricted_candidates_indexes[idx]].vertex);
      path_cost += insert_candidates[restricted_candidates_indexes[idx]].cost_increase;
      free_vertices.erase(std::find(free_vertices.begin(),free_vertices.end(),insert_candidates[restricted_candidates_indexes[idx]].vertex));
    }
    GenerateInsertMoves(cost_mat, rewards, max_cost, path, path_cost, insert_candidates, min, max);
  }
  return path_cost;
}

void Chromosome::GenerateGenes(Matrix<double_t> &cost_mat,
                               std::vector<double_t> &rewards,
                               std::vector<double_t> &max_cost_v,
                               const uint_fast32_t &start_vertex,
                               const uint_fast32_t &end_vertex) {
  // Initialise the genes
  std::vector<Gene>::iterator it;
  for (it = genes.begin(); it != genes.end(); ++it) {
    it->path.push_back(start_vertex);
    it->path.push_back(end_vertex);
    it->cost = cost_mat[start_vertex][end_vertex];
  }

  std::list<uint_fast32_t> free_list(free_vertices.begin(), free_vertices.end());
  std::list<uint_fast32_t>::iterator list_it = free_list.begin();
  while ((free_list.size() > 0) && (list_it != free_list.end())) {
    double_t reward = 0;
    reward = rewards[*list_it];
    double_t DoubleNAN = std::numeric_limits<double>::quiet_NaN();
    Insertion best(*list_it, 0, 0, 0, DoubleNAN, DoubleNAN, DoubleNAN);
    std::vector<uint_fast32_t>::iterator free_it;
    std::vector<Gene>::iterator gene_it;
    for (gene_it = genes.begin(); gene_it != genes.end(); ++gene_it) {
      Path::iterator path_it;
      for (path_it = gene_it->path.begin() + 1; path_it != gene_it->path.end(); ++path_it) {
        uint_fast32_t prev_vertex = *(path_it - 1);
        uint_fast32_t next_vertex = *path_it;
        double_t cost_increase =
            cost_mat[prev_vertex][*list_it] + cost_mat[*list_it][next_vertex] - cost_mat[prev_vertex][next_vertex] + 1;
        if (gene_it->cost + cost_increase < max_cost_v[0]) { // Ugly hack. Assumes all max costs are the same. Fixme!
          double_t avail_cost = max_cost_v[0] - (gene_it->cost + cost_increase);
          uint_fast32_t q = 0;
          for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
            if (*free_it != *list_it) {
              double_t dist = cost_mat[*list_it][*free_it];
              if (dist < avail_cost) {
                ++q;
              }
            }
          }
          double_t heuristic = q * (reward / cost_increase);
          if (std::isnan(best.heuristic) || best.heuristic < heuristic) {
            best.heuristic = heuristic;
            best.total_reward = reward;
            best.cost_increase = cost_increase;
            best.next_vertex = next_vertex;
            best.prev_vertex = prev_vertex;
            best.vertex = *list_it;
            best.gene = std::distance(genes.begin(), gene_it);
          }
        }
      }
    }
    if (std::isnan(best.heuristic)) {
      ++list_it;
    } else {
      genes[best.gene].path.insert(std::find(genes[best.gene].path.begin(),
                                             genes[best.gene].path.end(),
                                             best.next_vertex), best.vertex);
      genes[best.gene].cost += best.cost_increase;
      list_it = free_list.erase(list_it);
    }
  }
  free_vertices = std::vector<uint_fast32_t>(free_list.begin(), free_list.end());
}

double_t Chromosome::GenerateNNGRASPPath(const Matrix<double_t> &cost_mat,
                                         const std::vector<double_t> &rewards,
                                         const double_t &max_cost,
                                         const uint_fast32_t &start_vertex,
                                         const uint_fast32_t &end_vertex,
                                         std::mt19937 &g,
                                         Path &path) {
  std::uniform_real_distribution<> dis(0, 1);
  path.push_back(start_vertex);
  std::vector<uint_fast32_t> used_vertices;
  used_vertices.reserve(cost_mat.size());
  double_t path_cost = 0;
  bool done = false;
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
      expand_neighbours(neighbours, checked, all_vertices, cost_mat);
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
//        available = free_vertices;
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
      if (seen_vertices.size() > 0) {
        std::vector<uint_fast32_t>::iterator avail_it;
        for (avail_it = available.begin(); avail_it != available.end(); ++avail_it) {
          double mean_dist = 0;
          double free_neighbours = 0;
          std::unordered_set<uint_fast32_t>::iterator used_it;
          for (used_it = seen_vertices.begin(); used_it != seen_vertices.end(); ++used_it) {
            mean_dist += cost_mat[*avail_it][*used_it];
          }
          std::vector<uint_fast32_t>::iterator free_it;
          for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
            if (cost_mat[*avail_it][*free_it] < 2) {
              ++free_neighbours;
            }
          }
          mean_dist = mean_dist / double(seen_vertices.size());
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
      if (path_cost + next_vertex_cost + cost_mat[next_vertex][end_vertex] <= max_cost) {
        path_cost += next_vertex_cost;
        path.push_back(next_vertex);
        free_vertices.erase(std::remove(free_vertices.begin(), free_vertices.end(), next_vertex));
        seen_vertices.insert(next_vertex);
      } else {
        done = true;
      }
    } else {
      done = true;
    }
  }
  path.push_back(end_vertex);
  return path_cost;
}

Chromosome generate_chromosome(Matrix<double> &cost_mat,
                               std::vector<double> &max_cost_v,
                               std::vector<double> &rewards,
                               uint idx_start,
                               uint idx_finish,
                               std::mt19937 &g,
                               std::string gen_method) {
  uint n_agents = max_cost_v.size();
  size_t num_vertices = cost_mat.size();
  Chromosome c(num_vertices, max_cost_v.size(), idx_start, idx_finish);
  std::string method = "GRASP";
  if (method == "NNGRASP") {
    for (int i = 0; i < n_agents; ++i) {
      c.genes[i].cost =
          c.GenerateNNGRASPPath(cost_mat, rewards, max_cost_v[i], idx_start, idx_finish, g, c.genes[i].path);
    }
  } else if (method == "GRASP") {
    for (int i = 0; i < n_agents; ++i) {
      c.genes[i].cost =
          c.GenerateGRASPPath(cost_mat, rewards, max_cost_v[i], idx_start, idx_finish, g, c.genes[i].path);
    }
  } else {
    //method == GRASP2
    c.GenerateGenes(cost_mat, rewards, max_cost_v, idx_start, idx_finish);
  }

//  for (size_t i = 0; i < n_agents-1; ++i){
//    mutual_two_opt(c.genes[i].path,c.genes[i+1].path, cost_mat, max_cost_v[i], max_cost_v[i+1]);
//  }
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
      if (dist < 2) {
        nset.insert(*free_it);
      }
    }
  }
  neighbours.assign(nset.begin(), nset.end());
}

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

Chromosome ga_top(Matrix<double> &cost_mat,
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
  int tour_size = 3;
  int max_gen = 100;
  uint_fast32_t num_robots = max_cost_v.size();

  // Initialise population
  std::vector<Chromosome> pop;
  pop.reserve(pop_size);

  for (uint_fast32_t i = 0; i < pop_size; ++i) {
    Chromosome c = generate_chromosome(cost_mat, max_cost_v, rewards, idx_start, idx_finish, g, std::__cxx11::string());
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
//    for (int i = 0; i < 20; ++i) {
//      std::vector<size_t> indices = get_population_sample(new_pop.size(), 2, g);
//      cx(new_pop[indices[0]], new_pop[indices[1]], cost_mat, max_cost_v, rewards);
//    }

    // Mutate
    std::vector<size_t> indices = get_population_sample(new_pop.size(), 3*new_pop.size()/4, g);
    uint_fast32_t M = 8; //number of cores
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
//  std::vector<uint_fast32_t> vertices(cost_mat.size());
//  std::iota(vertices.begin(), vertices.end(), 0);
//  std::vector<uint_fast32_t> free_vertices;
//  std::unordered_set<uint_fast32_t> seen(best.genes[0].path.begin(), best.genes[0].path.end());
//  for (uint_fast32_t robot = 1; robot < num_robots; ++robot) {
//    for (uint_fast32_t vertex = 0; vertex < best.genes[robot].path.size(); ++vertex) {
//      seen.insert(vertex);
//    }
//  }
//  std::vector<uint_fast32_t> visited_vertices(seen.begin(), seen.end());
//  std::sort(visited_vertices.begin(), visited_vertices.end());
//  std::set_difference(vertices.begin(),
//                      vertices.end(),
//                      visited_vertices.begin(),
//                      visited_vertices.end(),
//                      std::back_inserter(free_vertices));
//  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
//    for (size_t vertex_idx = 1; vertex_idx < best.genes[robot].path.size() - 2; ++vertex_idx) {
//      std::vector<uint_fast32_t> available_vertices;
//      for (uint_fast32_t fv:free_vertices) {
//        if (cost_mat[best.genes[robot].path[vertex_idx]][fv] < 2) { //TODO: Fix the constant sensor distance
//          available_vertices.push_back(fv);
//        }
//      }
//
//      double best_cost = best.genes[robot].cost;
//      std::vector<uint_fast32_t> best_path = best.genes[robot].path;
//      double best_fitness = best.genes[robot].fitness;
//      uint_fast32_t best_vertex = best.genes[robot].path[vertex_idx];
//
//      for (uint_fast32_t i:available_vertices) {
//        Gene tmp_c;
//        tmp_c.path = best.genes[robot].path;
//        tmp_c.cost = best.genes[robot].cost
//            - cost_mat[best.genes[robot].path[vertex_idx - 1]][best.genes[robot].path[vertex_idx]]
//            - cost_mat[best.genes[robot].path[vertex_idx]][best.genes[robot].path[vertex_idx + 1]]
//            + cost_mat[best.genes[robot].path[vertex_idx - 1]][i]
//            + cost_mat[i][best.genes[robot].path[vertex_idx + 1]];
//        tmp_c.path.erase(tmp_c.path.begin() + vertex_idx);
//        tmp_c.path.insert(tmp_c.path.begin() + vertex_idx, i);
//        if (tmp_c.cost <= max_cost_v[robot]) {
//          tmp_c.evaluate_gene(cost_mat, rewards, free_vertices);
//          if (tmp_c.fitness > best_fitness) {
//            best_fitness = tmp_c.fitness;
//            best_cost = tmp_c.cost;
//            best_path = tmp_c.path;
//            best_vertex = i;
//          }
//        }
//      }
//      best.genes[robot].path = best_path;
//      best.genes[robot].fitness = best_fitness;
//      best.genes[robot].cost = best_cost;
//      if (best_vertex != best.genes[robot].path[vertex_idx]) {
//        free_vertices.push_back(best.genes[robot].path[vertex_idx]);
//        free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), best_vertex));
//      }
//    }
//  }

//  for (int j = 0; j < 100; ++j) {
//    best.mutate(cost_mat, rewards, max_cost_v, g);
//  }

  return best;
}
InsertMove::InsertMove(uint_fast32_t vertex,
                       uint_fast32_t prev_vertex,
                       uint_fast32_t next_vertex,
                       double cost_increase,
                       double total_reward,
                       double heuristic)
    : vertex(vertex), prev_vertex(prev_vertex), next_vertex(next_vertex), cost_increase(cost_increase),
      total_reward(total_reward), heuristic(heuristic) {
  ;
}
Insertion::Insertion(uint_fast32_t vertex,
                     uint_fast32_t prev_vertex,
                     uint_fast32_t next_vertex,
                     uint_fast32_t gene,
                     double cost_increase,
                     double total_reward,
                     double heuristic)
    : vertex(vertex), prev_vertex(prev_vertex), next_vertex(next_vertex), gene(gene), cost_increase(cost_increase),
      total_reward(total_reward), heuristic(heuristic) {

}