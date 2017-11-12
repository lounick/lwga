//
// Created by nick on 30/03/16.
//

#include <include/ctop_ga.h>

void
Gene::evaluate_gene(const Matrix<double> &cost_mat, const std::vector<double> &rewards, const
std::vector<uint_fast32_t> &free_vertices) {
  fitness = 0;

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
              * rewards[free_vertices[j]]; //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
        }
      }
      fitness += rewards[vertex] + extras;
    }
  }
  fitness = double(pow(fitness, 3)) / double(cost);
}

void Gene::calculate_cost(Matrix<double> &cost_mat) {
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
  Chromosome mutated(*this);
//  Chromosome mutated(cost_mat.size(), num_robots, genes.front().path.front(), genes.front().path.back());
//  mutated.genes.reserve(num_robots);
//  mutated.genes = genes;
  //mutated.evaluate_chromosome(cost_mat, rewards, max_cost_v);
//  mutated.total_fitness = total_fitness;

  std::uniform_real_distribution<> dis(0, 1);
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
    for (uint_fast32_t iter = 0; iter < 10; ++iter) {
      if (dis(g) < 0.9) {
        double max_cost = 0.95 * max_cost_v[robot];
        if (mutated.genes[robot].cost > max_cost || logically_equal(mutated.genes[robot].cost, max_cost)) {
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
          double best_fitness = mutated.genes[robot].fitness;//mutated.total_fitness;//mutated.genes[robot].fitness;
//          double best_path_fitness = mutated.genes[robot].fitness;
          uint_fast32_t best_vertex = rand_vertex;
          uint_fast32_t prev_vertex = mutated.genes[robot].path[rand_vertex_idx - 1];
          uint_fast32_t next_vertex = mutated.genes[robot].path[rand_vertex_idx + 1];
          double cost_removed = mutated.genes[robot].cost
              - cost_mat[prev_vertex][rand_vertex]
              - cost_mat[rand_vertex][next_vertex];
          free_vertices.push_back(rand_vertex);
          seen_vertices.erase(rand_vertex);
          for (uint_fast32_t i:available_vertices) {
            Gene tmp_c;
            tmp_c.path.reserve(cost_mat.size());
            tmp_c.path = mutated.genes[robot].path;
            tmp_c.cost = cost_removed
                + cost_mat[prev_vertex][i]
                + cost_mat[i][next_vertex];
            tmp_c.path.erase(tmp_c.path.begin() + rand_vertex_idx);
            tmp_c.path.insert(tmp_c.path.begin() + rand_vertex_idx, i);
//            free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), i));
//            seen_vertices.insert(i);
            if (tmp_c.cost < max_cost_v[robot] || logically_equal(tmp_c.cost, max_cost_v[robot])) {
//              double total_fitness = 0.0;
              tmp_c.evaluate_gene(cost_mat, rewards, free_vertices);
//              total_fitness += tmp_c.fitness;
//              for (uint_fast32_t gene_c = 0; gene_c < mutated.genes.size(); ++gene_c){
//                if (gene_c == robot)
//                  continue;
//                mutated.genes[gene_c].evaluate_gene(cost_mat, rewards, free_vertices);
//                total_fitness += mutated.genes[gene_c].fitness;
//              }
              if (tmp_c.fitness > best_fitness) {//total_fitness > best_fitness) {//tmp_c.fitness > best_fitness) {
                best_fitness = tmp_c.fitness;//total_fitness;//tmp_c.fitness;
//                best_path_fitness = tmp_c.fitness;
                best_cost = tmp_c.cost;
                best_path = tmp_c.path;
                best_vertex = i;
              }
            }
//            free_vertices.push_back(i);
//            seen_vertices.erase(i);
          }
          mutated.genes[robot].path = best_path;
          mutated.genes[robot].fitness = best_fitness;
          mutated.genes[robot].cost = best_cost;
//          mutated.total_fitness = best_fitness;
          if (best_vertex != rand_vertex) {
//            free_vertices.push_back(rand_vertex);
            free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), best_vertex));
            seen_vertices.insert(best_vertex);
//            seen_vertices.erase(rand_vertex);
          } else {
            free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), rand_vertex));
            seen_vertices.insert(rand_vertex);
          }
        } else {
          if (free_vertices.size() == 0)
            continue;
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
            double extras = 0;

            for (size_t j = 0; j < free_vertices.size(); ++j) {
              double dist = cost_mat[vertex][free_vertices[j]];
              if (dist < 2) {
                extras += std::exp((log(0.01) / 2) * dist)
                    * rewards[free_vertices[j]]; //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
              }
            }
            fit += extras;
            if (travel_increase != 0) {
              fit /= travel_increase;
            } else {
              fit = std::numeric_limits<double>::infinity();
            }
            if (fit > best_fitness) {
              double total_increase = mutated.genes[robot].cost + travel_increase + 1;
              if (total_increase < max_cost_v[robot] || logically_equal(total_increase, max_cost_v[robot])) {
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
        }
      } else {
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
              extras += std::exp((log(0.01) / 2) * cost_mat[mutated.genes[robot].path[i]][free_vertices[j]]);
            }
          }
          loss += extras;

          if (travel_decrease != 0) {
            loss /= travel_decrease;
          } else {
            loss = std::numeric_limits<double>::infinity();
          }

          if (loss < min_loss || logically_equal(loss, min_loss)) {
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
  genes = mutated.genes;
  evaluate_chromosome(cost_mat, rewards, max_cost_v);
}

Chromosome::Chromosome(size_t num_vertices, size_t num_genes, const uint_fast32_t &start_vertex,
                       const uint_fast32_t &end_vertex) {
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
    double_t extras = 0;
    std::vector<uint_fast32_t>::iterator it;
    for (it = free_vertices.begin(); it != free_vertices.end(); ++it) {
      if (*it != v) {
        double_t dist = cost_mat[v][*it];
        if (dist < 2) {
          //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
          extras += std::exp((log(0.01) / 2) * dist) * rewards[*it];
        }
      }
    }
    reward += rewards[v] + extras;

    //Generate all the possible insert moves of this vertex to this path
    Path::const_iterator p;
    for (p = path.begin() + 1; p != path.end(); ++p) {
      uint_fast32_t prev_vertex = *(p - 1);
      uint_fast32_t next_vertex = *p;
      double_t
          cost_increase = cost_mat[prev_vertex][v] + cost_mat[v][next_vertex] - cost_mat[prev_vertex][next_vertex] + 1;
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
      path.insert(std::find(path.begin(),
                            path.end(),
                            insert_candidates[restricted_candidates_indexes[idx]].next_vertex),
                  insert_candidates[restricted_candidates_indexes[idx]].vertex);
      path_cost += insert_candidates[restricted_candidates_indexes[idx]].cost_increase;
      free_vertices.erase(std::find(free_vertices.begin(),
                                    free_vertices.end(),
                                    insert_candidates[restricted_candidates_indexes[idx]].vertex));
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
    double_t extras = 0;
    std::vector<uint_fast32_t>::iterator free_it;
    for (free_it = free_vertices.begin(); free_it != free_vertices.end(); ++free_it) {
      if (*free_it != *list_it) {
        double_t dist = cost_mat[*list_it][*free_it];
        if (dist < 2) {
          //TODO: This assumes a fixed sensor range. Change it to what Valerio did in the quadratic COP problem solved by Gurobi.
          extras += std::exp((log(0.01) / 2) * dist) * rewards[*free_it];
        }
      }
    }
    reward += rewards[*list_it] + extras;
    double_t DoubleNAN = std::numeric_limits<double>::quiet_NaN();
    Insertion best(*list_it, 0, 0, 0, DoubleNAN, DoubleNAN, DoubleNAN);

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

void Chromosome::insertGene(Gene &gene){
  Path::iterator path_it;
  genes.push_back(gene);
  for (path_it = gene.path.begin(); path_it != gene.path.end(); ++path_it){
    free_vertices.erase(std::find(free_vertices.begin(), free_vertices.end(), *path_it));
    seen_vertices.insert(*path_it);
  }
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
      //==============================================================================================================
      double next_vertex_cost = cost_mat[path.back()][next_vertex] + 1;
      double next_vertex_end_cost = path_cost + next_vertex_cost + cost_mat[next_vertex][end_vertex];
      if (next_vertex_end_cost < max_cost || logically_equal(next_vertex_end_cost, max_cost)) {
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
double_t Chromosome::GenerateRandomPath(const Matrix<double_t> &cost_mat,
                                        const std::vector<double_t> &rewards,
                                        const double_t &max_cost,
                                        const uint_fast32_t &start_vertex,
                                        const uint_fast32_t &end_vertex,
                                        std::mt19937 &g,
                                        Path &path) {
  double_t path_cost = 0;
  bool done = false;
  path.push_back(start_vertex);
  while (!done) {
    if (free_vertices.size() == 0) {
      path_cost += cost_mat[path.back()][end_vertex];
      path.push_back(end_vertex);
      done = true;
    }
    uint_fast32_t free_idx = g() % free_vertices.size();
    uint_fast32_t to_insert = free_vertices[free_idx];
    double min_cost_to_end = path_cost + cost_mat[path.back()][to_insert] + cost_mat[to_insert][end_vertex] + 1;
    if (min_cost_to_end < max_cost || logically_equal(min_cost_to_end, max_cost)) {
      path_cost += cost_mat[path.back()][to_insert] + 1;
      path.push_back(to_insert);
      free_vertices.erase(free_vertices.begin() + free_idx);
      seen_vertices.insert(to_insert);
    } else {
      path_cost += cost_mat[path.back()][end_vertex];
      path.push_back(end_vertex);
      done = true;
    }
  }
  return path_cost;
}
void Chromosome::removeCommonVertices(const Gene &gene) {
  Path ::iterator path_it;
  for (path_it = gene.path.begin() + 1; path_it != gene.path.end() - 1; ++path_it){
    std::vector<Gene>::iterator genes_it;
    for (genes_it = genes.begin(); genes_it != genes.end(); ++genes_it){
      genes_it->path.erase(std::find(genes_it->path.begin(), genes_it->path.end(), *path_it));
    }
  }
}
void Chromosome::evaluateGenes(const Matrix<double> &cost_mat,
                               const std::vector<double> &rewards const,
                               const std::vector<uint_fast32_t> free_vertices) {
  std::vector<Gene>::iterator genes_it;
  for (genes_it = genes.begin(); genes_it != genes.end(); ++genes_it){
    genes_it->evaluate_gene(cost_mat, rewards, free_vertices);
  }
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
  std::string method = gen_method;
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
  } else if (method == "RANDOM") {
    for (int i = 0; i < n_agents; ++i) {
      c.genes[i].cost =
          c.GenerateRandomPath(cost_mat, rewards, max_cost_v[i], idx_start, idx_finish, g, c.genes[i].path);
    }
  } else {
    //method == GRASP2
    //c.GenerateGenes(cost_mat, rewards, max_cost_v, idx_start, idx_finish);
    std::cerr << "UNKNOWN GENERATION METHOD USING RANDOM" << std::endl;
    for (int i = 0; i < n_agents; ++i) {
      c.genes[i].cost =
          c.GenerateRandomPath(cost_mat, rewards, max_cost_v[i], idx_start, idx_finish, g, c.genes[i].path);
    }
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
      if (dist < 2) {
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
  Chromosome best_chromosome;
  if (population.size() < tour_size) {
    best_chromosome = population[0];
  } else {
    std::vector<uint> indices(population.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), g);
    double max_fitness = DBL_MIN;

    for (size_t i = 0; i < tour_size; i++) {
      if (population[indices[i]].total_fitness > max_fitness) {
        max_fitness = population[indices[i]].total_fitness;
        best_chromosome = population[indices[i]];
      }
    }
  }
  return best_chromosome;
}

void construct_offspring(Chromosome &offspring,
                         const Chromosome &parent1,
                         const Chromosome &parent2,
                         Matrix<double> &cost_mat,
                         std::vector<double> &max_cost_v,
                         std::vector<double> &rewards,
                         std::mt19937 &g) {

  Chromosome tmp_parent_1(parent1);
  Chromosome tmp_parent_2(parent2);
  uint_fast32_t num_robots = max_cost_v.size();
  std::uniform_real_distribution<> dis(0.0,1.0);
  auto sortRuleLambda = [](const Gene &g1, const Gene &g2) -> bool {
    return g1.fitness < g2.fitness;
  };

  offspring.genes.clear();
  offspring.total_fitness = 0;
  offspring.free_vertices = offspring.all_vertices;
  offspring.seen_vertices.clear();

  // First construct
  for (uint_fast32_t robot = 0; robot < num_robots; ++robot){
    Gene to_insert;
    if (dis(g) < 0.5){
      // Choose best gene from parent 1.
      // Remove that gene from parent 1.
      to_insert = tmp_parent_1.genes[0];
      tmp_parent_1.genes.erase(tmp_parent_1.genes.begin());
      // Insert gene to the offspring.
      offspring.insertGene(to_insert);
      // Remove vertices from genes of parent 2.
      tmp_parent_2.removeCommonVertices(to_insert);
      tmp_parent_2.evaluateGenes(cost_mat, rewards, offspring.free_vertices);
      std::sort(tmp_parent_2.genes.begin(), tmp_parent_2.genes.end(), sortRuleLambda);
    } else {
      // Choose best gene from parent 2.
      // Remove that gene from parent 2.
      to_insert = tmp_parent_2.genes[0];
      tmp_parent_2.genes.erase(tmp_parent_2.genes.begin());
      // Remove these vertices from the free and insert them to the seen.
      offspring.insertGene(to_insert);
      // Remove vertices from genes of parent 2.
      tmp_parent_1.removeCommonVertices(to_insert);
      tmp_parent_1.evaluateGenes(cost_mat, rewards, offspring.free_vertices);
      std::sort(tmp_parent_1.genes.begin(), tmp_parent_1.genes.end(), sortRuleLambda);
    }
  }
  //TODO: Then check if you can improve the offspring genes by inserting free vertices.
}

void cx(Chromosome &c1,
        Chromosome &c2,
        std::vector<std::vector<double> > &cost_mat,
        std::vector<double> &max_cost_v,
        std::vector<double> &rewards,
        std::mt19937 &g) {
  // CX As presented in "A new grouping genetic algorithm approach to the multiple traveling
  // salesperson problem" (or at least inspired by.
  auto sortRuleLambda = [](const Gene &g1, const Gene &g2) -> bool {
    return g1.fitness < g2.fitness;
  };

  std::sort(c1.genes.begin(), c1.genes.end(), sortRuleLambda);
  std::sort(c2.genes.begin(), c2.genes.end(), sortRuleLambda);
  Chromosome parent_1(c1);
  Chromosome parent_2(c2);
  construct_offspring(c1, parent_1, parent_2, cost_mat, max_cost_v, rewards, g);
  construct_offspring(c2, parent_1, parent_2, cost_mat, max_cost_v, rewards, g);
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
                   std::mt19937 &g,
                   size_t population_size = 100,
                   size_t n_generations = 50,
                   double_t cx_rate = 0.5,
                   double_t mutation_rate = 0.75,
                   std::string generation_method = std::string("GRASP"),
                   double_t elite_percent = 0.0) {
  /*
   * Initialise population
   * While gen < max_gen
   *  Select new pop
   *  Cx
   *  Mutate
   * Select fittest
  */

  uint_fast32_t pop_size = population_size;
  int tour_size = 3;
  int max_gen = n_generations;
  uint_fast32_t num_robots = max_cost_v.size();
  Chromosome best;
  if ((max_cost_v[0] - 2) <= cost_mat[idx_start][idx_finish]) {
    for (size_t i = 0; i < num_robots; ++i) {
      best.genes.push_back(Gene());
      best.genes[i].path.push_back(idx_start);
      best.genes[i].path.push_back(idx_finish);
    }
    return best;
  }

  // Initialise population
  std::vector<Chromosome> pop;
  pop.reserve(pop_size);

  for (uint_fast32_t i = 0; i < pop_size; ++i) {
    Chromosome c = generate_chromosome(cost_mat, max_cost_v, rewards, idx_start, idx_finish, g, generation_method);
    c.evaluate_chromosome(cost_mat, rewards, max_cost_v);
    pop.push_back(c);
  }

  for (int gen = 0; gen < max_gen; ++gen) {
    // Select new population
//    std::cout << "Calculating generation " << gen << std::endl;
    std::vector<Chromosome> new_pop;
    new_pop.reserve(pop_size);

    uint_fast32_t num_elites = std::ceil(elite_percent * pop_size);

    if (num_elites > 0) {
      auto sortRuleLambda = [](const Chromosome &c1, const Chromosome &c2) -> bool {
        return c1.total_fitness < c2.total_fitness;
      };
      std::sort(pop.begin(), pop.end(), sortRuleLambda);
      for (uint_fast32_t e = 0; e < num_elites; ++e) {
        new_pop.push_back(pop[e]);
      }
    }

    for (uint_fast32_t i = 0; i < pop_size - num_elites; ++i) {
      new_pop.push_back(tournament_select(pop, tour_size, g));
    }

//     Cx
    for (int i = 0; i < pop_size * cx_rate; ++i) {
      std::vector<size_t> indices = get_population_sample(new_pop.size(), 2, g);
      cx(new_pop[indices[0]], new_pop[indices[1]], cost_mat, max_cost_v, rewards);
    }

    // Mutate
    if (mutation_rate != 0) {
      std::vector<size_t>
          indices = get_population_sample(new_pop.size(), std::floor(mutation_rate * new_pop.size()), g);
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
    }

    pop = new_pop;
  }
  std::sort(pop.begin(), pop.end(), [](Chromosome c1, Chromosome c2) { return c1.total_fitness > c2.total_fitness; });
  best = pop[0];
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
