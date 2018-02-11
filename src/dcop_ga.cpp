//
// Created by nick on 12/01/18.
//

#include "dcop_ga.h"
#include <cassert>

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
  if (path.size() < 3) {
    fitness = 0;
    return;
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
//  std::tie(mutated.path, mutated.angles, mutated.cost) =
//      dubins_two_opt(
//          dubins_cost_mat, std_angles, mutated.nodes, mutated.rho, mutated.path,
//          mutated.angles, mutated.cost);
//  mutated.calculate_cost(dubins_cost_mat);
//  mutated.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
  for (uint_fast32_t iter = 0; iter < 10; ++iter){
    if (std::generate_canonical<double, 10>(g) < 0.9) {
      //Try to add
      if (mutated.cost > 0.95*max_cost) {
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
//        uint_fast32_t prev_vertex = mutated.path[rand_vertex_idx - 1];
//        uint_fast32_t prev_angle = mutated.angles[rand_vertex_idx - 1];
//        uint_fast32_t next_vertex = mutated.path[rand_vertex_idx + 1];
//        uint_fast32_t next_angle = mutated.angles[rand_vertex_idx + 1];
//        double_t cost_removed = mutated.cost
//            - dubins_cost_mat[prev_vertex][rand_vertex][prev_angle][rand_angle]
//            - dubins_cost_mat[rand_vertex][next_vertex][rand_angle][next_angle];
        for (uint_fast32_t i:available_vertices) {
          Chromosome tmp_c(mutated);
          tmp_c.path[rand_vertex_idx] = i;
          uint_fast32_t pa = tmp_c.angles[rand_vertex_idx-1];
          uint_fast32_t na = tmp_c.angles[rand_vertex_idx+1];
          uint8_t angles_num = 5;
          Vector<uint_fast32_t> p_angles, n_angles;
          p_angles.reserve(angles_num);
          n_angles.reserve(angles_num);
          for (int32_t angle_cnt = -2; angle_cnt < 3; ++angle_cnt) {
            if (int32_t(pa) + angle_cnt < 0) {
              p_angles.push_back(std_angles.size() + pa + angle_cnt);
            } else if (pa + angle_cnt >= std_angles.size()) {
              p_angles.push_back(pa + angle_cnt - std_angles.size());
            } else {
              p_angles.push_back(pa + angle_cnt);
            }
            if (int32_t(na) + angle_cnt < 0) {
              n_angles.push_back(std_angles.size() + na + angle_cnt);
            } else if (na + angle_cnt >= std_angles.size()) {
              n_angles.push_back(na + angle_cnt - std_angles.size());
            } else {
              n_angles.push_back(na + angle_cnt);
            }
          }
          uint_fast32_t best_pa, best_a, best_na;
          double_t best_travel_increase = std::numeric_limits<double_t>::max();
          if (rand_vertex_idx == 1) {
            //TODO:Fix case where only one vertext exists in path e.g. 0-1-26 in 5x5
            uint_fast32_t pv = tmp_c.path[rand_vertex_idx-1];
            uint_fast32_t nv = tmp_c.path[rand_vertex_idx+1];
            uint_fast32_t nnv = tmp_c.path[rand_vertex_idx+2];
            uint_fast32_t nna = tmp_c.angles[rand_vertex_idx+2];
//            for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
            for (uint_fast32_t pa:p_angles) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a) {
//                for (uint_fast32_t na = 0; na < std_angles.size(); ++na) {
                for (uint_fast32_t na:n_angles) {
                  double_t travel_increase = 0.0;
                  travel_increase += dubins_cost_mat[pv][i][pa][a];
                  travel_increase += dubins_cost_mat[i][nv][a][na];
                  if (rand_vertex_idx + 2 < tmp_c.path.size()) {
                    travel_increase += dubins_cost_mat[nv][nnv][na][nna];
                  }
                  if (travel_increase < best_travel_increase) {
                    best_travel_increase = travel_increase;
                    best_pa = pa;
                    best_a = a;
                    best_na = na;
                  }
                }
              }
            }
          } else if (rand_vertex_idx == mutated.path.size()-2) {
            uint_fast32_t ppv = tmp_c.path[rand_vertex_idx-2];
            uint_fast32_t ppa = tmp_c.angles[rand_vertex_idx-2];
            uint_fast32_t pv = tmp_c.path[rand_vertex_idx-1];
            uint_fast32_t nv = tmp_c.path[rand_vertex_idx+1];
//            for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
            for (uint_fast32_t pa:p_angles) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a) {
//                for (uint_fast32_t na = 0; na < std_angles.size(); ++na) {
                for (uint_fast32_t na:n_angles) {
                  double_t travel_increase = 0.0;
                  if (rand_vertex_idx - 2 >= 0) {
                    travel_increase += dubins_cost_mat[ppv][pv][ppa][pa];
                  }
                  travel_increase += dubins_cost_mat[pv][i][pa][a];
                  travel_increase += dubins_cost_mat[i][nv][a][na];
                  if (travel_increase < best_travel_increase) {
                    best_travel_increase = travel_increase;
                    best_pa = pa;
                    best_a = a;
                    best_na = na;
                  }
                }
              }
            }
          } else {
            uint_fast32_t ppv = tmp_c.path[rand_vertex_idx-2];
            uint_fast32_t ppa = tmp_c.angles[rand_vertex_idx-2];
            uint_fast32_t pv = tmp_c.path[rand_vertex_idx-1];
            uint_fast32_t nv = tmp_c.path[rand_vertex_idx+1];
            uint_fast32_t nnv = tmp_c.path[rand_vertex_idx+2];
            uint_fast32_t nna = tmp_c.angles[rand_vertex_idx+2];
//            for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
            for (uint_fast32_t pa:p_angles) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a) {
//                for (uint_fast32_t na = 0; na < std_angles.size(); ++na) {
                for (uint_fast32_t na:n_angles) {
                  double_t travel_increase = 0.0;
                  travel_increase += dubins_cost_mat[ppv][pv][ppa][pa];
                  travel_increase += dubins_cost_mat[pv][i][pa][a];
                  travel_increase += dubins_cost_mat[i][nv][a][na];
                  travel_increase += dubins_cost_mat[nv][nnv][na][nna];
                  if (travel_increase < best_travel_increase) {
                    best_travel_increase = travel_increase;
                    best_pa = pa;
                    best_a = a;
                    best_na = na;
                  }
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
          mutated.seen_vertices.erase(rand_vertex);
        }
      } else {
        if (free_vertices.size() == 0)
          continue;
        // Choose the vertex that has the most free neighbours
        // Choose the neighbour that has the most free neighbours
        // Try to add the neighbour before and after the vertex
        // If max cost is not violated
        // Choose addition with minimal cost increase
        std::uniform_int_distribution<> dis(0, mutated.free_vertices.size() - 1);
        size_t vertex_idx = dis(g);
        uint_fast32_t vertex = mutated.free_vertices[vertex_idx];

        double_t best_travel_increase = 0.0;
        double_t best_removed_cost = 0.0;
        double_t best_fitness = 0.0;
        uint_fast32_t best_pa, best_a, best_na;
        size_t ins_pos = 0;
        for (size_t i = 1; i < mutated.path.size(); ++i) {
          double_t travel_increase_pos = std::numeric_limits<double_t>::max();
          uint_fast32_t best_pa_pos, best_a_pos, best_na_pos;
          if (i == 1) {
            uint_fast32_t pv = mutated.path[i-1];
            uint_fast32_t nv = mutated.path[i];
            uint_fast32_t nnv = mutated.path[i+1];
            uint_fast32_t nna = mutated.angles[i+1];
            for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a){
                //TODO: Maybe consider a subset of angles
                for (uint_fast32_t na = 0; na < std_angles.size(); ++na){
                  double_t travel_increase = 0.0;
                  travel_increase += dubins_cost_mat[pv][vertex][pa][a];
                  travel_increase += dubins_cost_mat[vertex][nv][a][na];
                  travel_increase += dubins_cost_mat[nv][nnv][na][nna];
                  if (travel_increase < travel_increase_pos){
                    travel_increase_pos = travel_increase;
                    best_pa_pos = pa;
                    best_a_pos = a;
                    best_na_pos = na;
                  }
                }
              }
            }
          } else if (i == mutated.path.size()-1) {
            uint_fast32_t ppa = mutated.angles[i-2];
            uint_fast32_t ppv = mutated.path[i-2];
            uint_fast32_t pv = mutated.path[i-1];
            uint_fast32_t nv = mutated.path[i];
            for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a) {
                for (uint_fast32_t na = 0; na < std_angles.size(); ++na) {
                  double_t travel_increase = 0.0;
                  if (i - 2 >= 0) {
                    travel_increase += dubins_cost_mat[ppv][pv][ppa][pa];
                  }
                  travel_increase += dubins_cost_mat[pv][vertex][pa][a];
                  travel_increase += dubins_cost_mat[vertex][nv][a][na];
                  if (travel_increase < travel_increase_pos){
                    travel_increase_pos = travel_increase;
                    best_pa_pos = pa;
                    best_a_pos = a;
                    best_na_pos = na;
                  }
                }
              }
            }
          } else {
            uint_fast32_t ppa = mutated.angles[i-2];
            uint_fast32_t nna = mutated.angles[i+1];
            uint_fast32_t ppv = mutated.path[i-2];
            uint_fast32_t pv = mutated.path[i-1];
            uint_fast32_t nv = mutated.path[i];
            uint_fast32_t nnv = mutated.path[i+1];
            for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
              for (uint_fast32_t a = 0; a < std_angles.size(); ++a) {
                for (uint_fast32_t na = 0; na < std_angles.size(); ++na) {
                  double_t travel_increase = 0.0;
                  if (i - 2 >= 0) {
                    travel_increase += dubins_cost_mat[ppv][pv][ppa][pa];
                  }
                  travel_increase += dubins_cost_mat[pv][vertex][pa][a];
                  travel_increase += dubins_cost_mat[vertex][nv][a][na];
                  travel_increase += dubins_cost_mat[nv][nnv][na][nna];
                  if (travel_increase < travel_increase_pos){
                    travel_increase_pos = travel_increase;
                    best_pa_pos = pa;
                    best_a_pos = a;
                    best_na_pos = na;
                  }
                }
              }
            }
          }

          double fit = rewards[vertex]; //TODO: Should we add the extras to the fit?
          if (travel_increase_pos != 0) {
            fit /= travel_increase_pos;
          }
          else {
            fit = std::numeric_limits<double>::infinity();
          }
          if(fit > best_fitness) {
            double_t removed_cost = 0.0;
            if (i == 1) {
              uint_fast32_t pv = mutated.path[i-1];
              uint_fast32_t pa = mutated.angles[i-1];
              uint_fast32_t nv = mutated.path[i];
              uint_fast32_t na = mutated.angles[i];
              uint_fast32_t nnv = mutated.path[i+1];
              uint_fast32_t nna = mutated.angles[i+1];
              removed_cost += dubins_cost_mat[pv][nv][pa][na];
              removed_cost += dubins_cost_mat[nv][nnv][na][nna];
            } else if (i == mutated.path.size() - 1) {
              uint_fast32_t ppv = mutated.path[i-2];
              uint_fast32_t ppa = mutated.angles[i-2];
              uint_fast32_t pv = mutated.path[i-1];
              uint_fast32_t pa = mutated.angles[i-1];
              uint_fast32_t nv = mutated.path[i];
              uint_fast32_t na = mutated.angles[i];
              if (i - 2 >= 0) {
                removed_cost += dubins_cost_mat[ppv][pv][ppa][pa];
              }
              removed_cost += dubins_cost_mat[pv][nv][pa][na];
            } else {
              uint_fast32_t ppv = mutated.path[i-2];
              uint_fast32_t ppa = mutated.angles[i-2];
              uint_fast32_t pv = mutated.path[i-1];
              uint_fast32_t pa = mutated.angles[i-1];
              uint_fast32_t nv = mutated.path[i];
              uint_fast32_t na = mutated.angles[i];
              uint_fast32_t nnv = mutated.path[i+1];
              uint_fast32_t nna = mutated.angles[i+1];
              if (i - 2 >= 0) {
                removed_cost += dubins_cost_mat[ppv][pv][ppa][pa];
              }
              removed_cost += dubins_cost_mat[pv][nv][pa][na];
              removed_cost += dubins_cost_mat[nv][nnv][na][nna];
            }
            if (mutated.cost - removed_cost + travel_increase_pos <= max_cost) {
              best_travel_increase = travel_increase_pos;
              best_removed_cost = removed_cost;
              best_fitness = fit;
              best_pa = best_pa_pos;
              best_a = best_a_pos;
              best_na = best_na_pos;
              ins_pos = i;
            }
          }
        }
        if (ins_pos > 0) {
          /*if (ins_pos == 1) {
            uint_fast32_t pv = mutated.path[ins_pos-1];
            uint_fast32_t pa = mutated.angles[ins_pos-1];
            uint_fast32_t nv = mutated.path[ins_pos];
            uint_fast32_t na = mutated.angles[ins_pos];
            uint_fast32_t nnv = mutated.path[ins_pos+1];
            uint_fast32_t nna = mutated.angles[ins_pos+1];
            mutated.cost -= dubins_cost_mat[pv][nv][pa][na];
            mutated.cost -= dubins_cost_mat[nv][nnv][na][nna];
            mutated.cost += dubins_cost_mat[pv][vertex][best_pa][best_a];
            mutated.cost += dubins_cost_mat[vertex][nv][best_a][best_na];
            mutated.cost += dubins_cost_mat[nv][nnv][best_na][nna];
          } else if (ins_pos == mutated.path.size() - 1) {
            uint_fast32_t ppv = mutated.path[ins_pos-2];
            uint_fast32_t ppa = mutated.angles[ins_pos-2];
            uint_fast32_t pv = mutated.path[ins_pos-1];
            uint_fast32_t pa = mutated.angles[ins_pos-1];
            uint_fast32_t nv = mutated.path[ins_pos];
            uint_fast32_t na = mutated.angles[ins_pos];
            mutated.cost -= dubins_cost_mat[ppv][pv][ppa][pa];
            mutated.cost -= dubins_cost_mat[pv][nv][pa][na];
            mutated.cost += dubins_cost_mat[ppv][pv][ppa][best_pa];
            mutated.cost += dubins_cost_mat[pv][vertex][best_pa][best_a];
            mutated.cost += dubins_cost_mat[vertex][nv][best_a][best_na];
          } else {
            uint_fast32_t ppv = mutated.path[ins_pos-2];
            uint_fast32_t ppa = mutated.angles[ins_pos-2];
            uint_fast32_t pv = mutated.path[ins_pos-1];
            uint_fast32_t pa = mutated.angles[ins_pos-1];
            uint_fast32_t nv = mutated.path[ins_pos];
            uint_fast32_t na = mutated.angles[ins_pos];
            uint_fast32_t nnv = mutated.path[ins_pos+1];
            uint_fast32_t nna = mutated.angles[ins_pos+1];
            mutated.cost -= dubins_cost_mat[ppv][pv][ppa][pa];
            mutated.cost -= dubins_cost_mat[pv][nv][pa][na];
            mutated.cost -= dubins_cost_mat[nv][nnv][na][nna];
            mutated.cost += dubins_cost_mat[ppv][pv][ppa][best_pa];
            mutated.cost += dubins_cost_mat[pv][vertex][best_pa][best_a];
            mutated.cost += dubins_cost_mat[vertex][nv][best_a][best_na];
            mutated.cost += dubins_cost_mat[nv][nnv][best_na][nna];
          }*/
          mutated.cost -= best_removed_cost;
          mutated.cost += best_travel_increase;
          mutated.angles[ins_pos-1] = best_pa;
          mutated.angles[ins_pos] = best_na;
          mutated.path.insert(mutated.path.begin() + ins_pos, vertex);
          mutated.angles.insert(mutated.angles.begin() + ins_pos, best_a);
          mutated.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
          mutated.free_vertices.erase(mutated.free_vertices.begin() + vertex_idx);
          mutated.seen_vertices.insert(vertex);
        }
      }
    } else {
      //Remove
      if (mutated.cost >= 0.5 * max_cost) {
        size_t to_remove = 0;
        double min_loss = std::numeric_limits<double>::infinity();
        uint_fast32_t best_pa, best_na;
        for (size_t i = 1; i < mutated.path.size() - 1; ++i) {
          uint_fast32_t v = mutated.path[i];
          uint_fast32_t a = mutated.angles[i];
          uint_fast32_t pv= mutated.path[i-1];
          uint_fast32_t pa = mutated.angles[i-1];
          uint_fast32_t nv = mutated.path[i+1];
          uint_fast32_t na = mutated.angles[i+1];
          double travel_decrease = dubins_cost_mat[pv][v][pa][a]
              + dubins_cost_mat[v][nv][a][na] - dubins_cost_mat[pv][nv][pa][na];
//          double travel_decrease = dubins_cost_mat[pv][v][pa][a]
//              + dubins_cost_mat[v][nv][a][na];
//          double_t  min_connection = std::numeric_limits<double_t>::max();
//          uint_fast32_t best_connection_pa, best_connection_na;
//          for (uint_fast32_t pa = 0; pa < std_angles.size(); ++pa) {
//            for (uint_fast32_t na = 0; na < std_angles.size(); ++na) {
//              double_t connection = dubins_cost_mat[pv][nv][pa][na];
//              if (connection < min_connection) {
//                min_connection = connection;
//                best_connection_na = na;
//                best_connection_pa = pa;
//              }
//            }
//          }
//          travel_decrease -= dubins_cost_mat[pv][nv][best_connection_pa][best_connection_na];
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
//            best_pa = best_connection_pa;
//            best_na = best_connection_na;
            to_remove = i;
          }
        }
        if (to_remove > 0) {
//          mutated.angles[to_remove-1] = best_pa;
//          mutated.angles[to_remove+1] = best_na;
          mutated.free_vertices.push_back(mutated.path[to_remove]);
          mutated.seen_vertices.erase(mutated.path[to_remove]);
          mutated.path.erase(mutated.path.begin() + to_remove);
          mutated.angles.erase(mutated.angles.begin() + to_remove);
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

void cx(Chromosome &c1,
        Chromosome &c2,
        Matrix<Matrix<double_t>>&dubins_cost_mat,
        Vector<double_t> &std_angles,
        Matrix<double> &cost_mat,
        Vector<double_t> &rewards,
        double_t max_cost,std::mt19937 &g){

  std::vector<int> intersection;
  std::set_intersection(
      c1.seen_vertices.begin(), c1.seen_vertices.end(),
      c2.seen_vertices.begin(), c2.seen_vertices.end(),
      std::back_inserter(intersection));
  if (intersection.empty()) {
    return;
  }
  Chromosome off1;
  Chromosome off2;
  off1.path.reserve(cost_mat.size());
  off1.path.clear();
  off1.angles.reserve(cost_mat.size());
  off1.angles.clear();
  off2.path.reserve(cost_mat.size());
  off2.path.clear();
  off2.angles.reserve(cost_mat.size());
  off2.angles.clear();
  off1.all_vertices = c1.all_vertices;
  off2.all_vertices = c2.all_vertices;
  off1.nodes = c1.nodes;
  off2.nodes = c2.nodes;
  off1.rho = c1.rho;
  off2.rho = c2.rho;
  //1. Choose random intersection point.
  size_t rand_idx = g() % intersection.size();
  //2. Exchange paths.
  //3. Remove duplicates and populate seen/free vertices.
  //4. Connect first part of path with old angles.
  off1.seen_vertices = std::unordered_set<uint_fast32_t> ();
  size_t c1_idx = 0;
  size_t c2_idx = 0;
  while (c1.path[c1_idx] != intersection[rand_idx]) {
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> ret;
    ret = off1.seen_vertices.insert(c1.path[c1_idx]);
    if (ret.second) {
      off1.path.push_back(c1.path[c1_idx]);
    }
    ++c1_idx;
  }
  off1.seen_vertices.insert(c1.path[c1_idx]);
  off1.path.push_back(c1.path[c1_idx]);
  for (size_t i = 0; i < off1.path.size(); ++i) {
    off1.angles.push_back(c1.angles[i]);
  }

  off2.seen_vertices = std::unordered_set<uint_fast32_t> ();
  while (c2.path[c2_idx] != intersection[rand_idx]){
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> ret;
    ret = off2.seen_vertices.insert(c2.path[c2_idx]);
    if (ret.second) {
      off2.path.push_back(c2.path[c2_idx]);
    }
    ++c2_idx;
  }
  off2.seen_vertices.insert(c2.path[c2_idx]);
  off2.path.push_back(c2.path[c2_idx]);
  for (size_t i = 0; i < off2.path.size(); ++i) {
    off2.angles.push_back(c2.angles[i]);
  }

  for (size_t i = c2_idx + 1; i < c2.path.size(); ++i) {
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> ret;
    ret = off1.seen_vertices.insert(c2.path[i]);
    if (ret.second) {
      off1.path.push_back(c2.path[i]);
    }
  }
  for (size_t i = c1_idx + 1; i < c1.path.size(); ++i) {
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> ret;
    ret = off2.seen_vertices.insert(c1.path[i]);
    if (ret.second) {
      off2.path.push_back(c1.path[i]);
    }
  }

  off1.free_vertices.reserve(cost_mat.size());
  off1.seen_vertices.erase(off1.path.back());
  off1.seen_vertices.erase(off1.path.front());
  for (auto v:off1.all_vertices) {
    if (!off1.seen_vertices.count(v)) {
      off1.free_vertices.push_back(v);
    }
  }
  off1.free_vertices.erase(
      std::find(off1.free_vertices.begin(), off1.free_vertices.end(), off1.path.front()));
  off1.free_vertices.erase(
      std::find(off1.free_vertices.begin(), off1.free_vertices.end(), off1.path.back()));

  off2.seen_vertices.erase(off2.path.back());
  off2.seen_vertices.erase(off2.path.front());
  off2.free_vertices.reserve(cost_mat.size());
  for (auto v:off2.all_vertices) {
    if (!off2.seen_vertices.count(v)) {
      off2.free_vertices.push_back(v);
    }
  }
  off2.free_vertices.erase(
      std::find(off2.free_vertices.begin(), off2.free_vertices.end(), off2.path.front()));
  off2.free_vertices.erase(
      std::find(off2.free_vertices.begin(), off2.free_vertices.end(), off2.path.back()));
  //5. Find new angles for second part of path.
  auto idx1 = std::find(off1.path.begin(), off1.path.end(), intersection[rand_idx]);
  auto idx2 = std::find(off2.path.begin(), off2.path.end(), intersection[rand_idx]);
  uint_fast32_t pa;
  uint_fast32_t best_a, best_na;
  best_na = c1.angles.back();
  pa = off1.angles.back();
  for (; idx1 < off1.path.end()-2; ++idx1){
    double_t best_cost = std::numeric_limits<double_t>::max();
    uint_fast32_t a, na, nna;
    for (a = 0; a < std_angles.size(); ++a) {
      for (na = 0; na < std_angles.size(); ++na) {
        double_t cost = 0.0;
        cost += dubins_cost_mat[*idx1][*(idx1+1)][pa][a];
        cost += dubins_cost_mat[*(idx1+1)][*(idx1+2)][a][na];
        if (cost < best_cost) {
          best_cost = cost;
          best_a = a;
          best_na = na;
        }
      }
    }
    pa = best_a;
    off1.angles.push_back(best_a);
  }
  off1.angles.push_back(best_na);

  pa = off2.angles.back();
  best_na = c2.angles.back();
  for (; idx2 < off2.path.end()-2; ++idx2){
    double_t best_cost = std::numeric_limits<double_t>::max();
    uint_fast32_t a, na, nna;
    for (a = 0; a < std_angles.size(); ++a) {
      for (na = 0; na < std_angles.size(); ++na) {
        double_t cost = 0.0;
        cost += dubins_cost_mat[*idx2][*(idx2+1)][pa][a];
        cost += dubins_cost_mat[*(idx2+1)][*(idx2+2)][a][na];
        if (cost < best_cost) {
          best_cost = cost;
          best_a = a;
          best_na = na;
        }
      }
    }
    pa = best_a;
    off2.angles.push_back(best_a);
  }
  off2.angles.push_back(best_na);
  assert(off1.angles.size() == off1.path.size());
  assert(off2.angles.size() == off2.path.size());
  //6. Do 2-opt
  off1.calculate_cost(dubins_cost_mat);
  off2.calculate_cost(dubins_cost_mat);
//  std::tie(off1.path, off1.angles, off1.cost) = dubins_two_opt(
//      dubins_cost_mat, std_angles, off1.nodes, off1.rho, off1.path,
//      off1.angles, off1.cost);
//  std::tie(off2.path, off2.angles, off2.cost) = dubins_two_opt(
//      dubins_cost_mat, std_angles, off2.nodes, off2.rho, off2.path,
//      off2.angles, off2.cost);
  //7. Check feasibility.
  if (off1.cost < max_cost) {
    c1 = off1;
  }
  if (off2.cost < max_cost) {
    c2 = off2;
  }
  //8. Calculate cost and fitness.
  c1.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
  c2.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
//  print_vector(c1.path);
//  print_vector(c2.path);
}

void par_cx(std::vector<size_t> indices,
            std::vector<Chromosome> &pop,
            Matrix<Matrix<double_t>>&dubins_cost_mat,
            Vector<double_t> &std_angles,
            Matrix<double_t> &cost_mat,
            std::vector<double_t> &rewards,
            double_t max_cost) {
  // Get hash of thread id for the seed of the generator.
  std::hash<std::thread::id> hasher;
  static thread_local std::mt19937 g(hasher(std::this_thread::get_id()));
  for (size_t i = 0; i < indices.size(); i = i + 2)
    cx(pop[indices[i]], pop[indices[i + 1]], dubins_cost_mat, std_angles, cost_mat, rewards, max_cost, g);
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
                   uint_fast16_t pop_size,
                   uint_fast8_t num_gen,
                   uint_fast8_t tour_size,
                   double_t cx_rate,
                   double_t mut_rate,
                   double_t elitist_rate,
                   std::mt19937 &g){
  /*
   * Initialise population
   * While gen < max_gen
   *  Select new pop
   *  Cx
   *  Mutate
   * Select fittest
  */

  Chromosome best;
  best.fitness = std::numeric_limits<double_t>::min();
  uint8_t num_best = 0;
//  std::cout << max_cost << " " << cost_mat[idx_start][idx_finish] << std::endl;
  if ((max_cost) <= cost_mat[idx_start][idx_finish]+3) {
    best.path.push_back(idx_start);
    best.path.push_back(idx_finish);
    return best;
  }

  // Initialise population
  std::vector<Chromosome> pop;
  pop.reserve(pop_size);

  for (uint_fast32_t i = 0; i < pop_size; ++i) {
//    std::cout << "Generating chromosome: " << i << std::endl;
    Chromosome c = generate_chromosome(
        nodes, std_angles, rho, max_cost,
        idx_start, idx_finish, dubins_cost_mat, cost_mat, g);
//    std::tie(c.angles, c.cost) = straighten_path(dubins_cost_mat, nodes, rho, c.path, c.angles, c.cost);
//    std::tie(c.path, c.angles, c.cost) =
//        dubins_two_opt(dubins_cost_mat, std_angles, nodes, rho, c.path, c.angles, c.cost);
//    std::tie(c.angles, c.cost) = straighten_path(dubins_cost_mat, nodes, rho, c.path, c.angles, c.cost);
    c.calculate_cost(dubins_cost_mat);
    c.evaluate_chromosome(dubins_cost_mat, cost_mat, rewards);
    pop.push_back(c);
  }

  uint8_t num_stable = 0;
  uint8_t max_stable = 10;

  for (int gen = 0; gen < num_gen && num_stable < max_stable; ++gen) {
    // Select new populations
    std::vector<Chromosome> new_pop;
    new_pop.reserve(pop_size);
    uint_fast32_t num_elites = std::ceil(elitist_rate * pop_size);
    if (num_elites > 0) {
      auto sortRuleLambda = [](const Chromosome &c1, const Chromosome &c2) -> bool {
        return c1.fitness > c2.fitness;
      };
      std::sort(pop.begin(), pop.end(), sortRuleLambda);
      for (uint_fast32_t e = 0; e < num_elites; ++e) {
        new_pop.push_back(pop[e]);
      }
    }

    if (logically_equal(best.fitness, new_pop.front().fitness)) {
      ++num_stable;
    } else {
      best = new_pop.front();
      num_stable = 0;
    }

    for (uint_fast32_t i = 0; i < pop_size - num_elites; ++i) {
      new_pop.push_back(tournament_select(pop, tour_size, g));
    }

//    for (Chromosome c:new_pop)
//      print_vector(c.path);

    // Cx
//    for (int i = 0; i < std::floor(0.9*pop_size); ++i) {
//      std::vector<size_t> indices = get_population_sample(new_pop.size(), 2, g);
//      std::pair<Chromosome, Chromosome> cx_ret =
//          cx(new_pop[indices[0]], new_pop[indices[1]], dubins_cost_mat,
//             std_angles, cost_mat, rewards, max_cost, g);
//      new_pop[indices[0]] = cx_ret.first;
//      new_pop[indices[1]] = cx_ret.second;
//    }
    if (cx_rate != 0) {
      int num_individuals = std::ceil(pop_size * cx_rate);
      if (num_individuals % 2)
        ++num_individuals;
      std::vector<size_t> indices = get_population_sample(new_pop.size(), num_individuals, g);
      for (size_t i = 0; i < indices.size(); i = i + 2)
        cx(new_pop[indices[i]], new_pop[indices[i + 1]], dubins_cost_mat, std_angles, cost_mat, rewards, max_cost, g);
      uint_fast32_t M = 6; //number of cores
      uint_fast32_t chunk_size = indices.size() / M;
      if (chunk_size % 2)
        --chunk_size;

      std::vector<std::future<void> > future_v;
      for (uint_fast32_t thread_count = 0; thread_count < M; ++thread_count) {
        std::vector<size_t>
            tmpv(indices.begin() + thread_count * chunk_size, indices.begin() + (thread_count + 1) * chunk_size);
        future_v.push_back(std::async(std::launch::async,
                                      par_cx,
                                      tmpv,
                                      std::ref(new_pop),
                                      std::ref(dubins_cost_mat),
                                      std::ref(std_angles),
                                      std::ref(cost_mat),
                                      std::ref(rewards),
                                      std::ref(max_cost)));
      }
      if (indices.size() > chunk_size*M){
        std::vector<size_t> tmpv(indices.begin() + (M) * chunk_size, indices.end());
        future_v.push_back(std::async(std::launch::async,
                                      par_cx,
                                      tmpv,
                                      std::ref(new_pop),
                                      std::ref(dubins_cost_mat),
                                      std::ref(std_angles),
                                      std::ref(cost_mat),
                                      std::ref(rewards),
                                      std::ref(max_cost)));
      }
      for (auto &f: future_v) {
        f.get();
      }
    }
    // Mutate
    if (mut_rate != 0) {
      std::vector<size_t> indices =
          get_population_sample(new_pop.size(), std::floor(mut_rate * pop_size), g);
      uint_fast32_t M = 6; //number of cores
      uint_fast32_t chunk_size = indices.size() / M;

//    for (size_t idx : indices) {
//      new_pop[idx].mutate(dubins_cost_mat, cost_mat, std_angles, rewards, max_cost, g);
//    }
      std::vector<std::future<void> > future_v;
      for (uint_fast32_t thread_count = 0; thread_count < M; ++thread_count) {
        //std::launch::deferred|std::launch::async;
        std::vector<size_t> tmpv(indices.begin() + thread_count * chunk_size,
                                 indices.begin()
                                     + (thread_count + 1) * chunk_size);
        future_v.push_back(std::async(std::launch::async,
                                      par_mutate,
                                      tmpv,
                                      std::ref(new_pop),
                                      std::ref(dubins_cost_mat),
                                      std::ref(cost_mat),
                                      std::ref(std_angles),
                                      std::ref(rewards),
                                      std::ref(max_cost)));
      }
      if (indices.size() % M != 0) {
        std::vector<size_t>
            tmpv(indices.begin() + (M) * chunk_size, indices.end());
        future_v.push_back(std::async(std::launch::async,
                                      par_mutate,
                                      tmpv,
                                      std::ref(new_pop),
                                      std::ref(dubins_cost_mat),
                                      std::ref(cost_mat),
                                      std::ref(std_angles),
                                      std::ref(rewards),
                                      std::ref(max_cost)));
      }
      for (auto &f: future_v) {
        f.get();
      }
    }
    pop = new_pop;
  }
  std::sort(pop.begin(), pop.end(), [](Chromosome c1, Chromosome c2) { return c1.fitness > c2.fitness; });
  if (best.fitness < pop[0].fitness)
    best = pop[0];
  return best;
}

}