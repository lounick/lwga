//
// Created by nick on 12/01/18.
//

#include "dcop_ga.h"

namespace dcop_ga {

void Chromosome::calculate_cost() {
  cost = get_dubins_path_cost(nodes, rho, path, angles);
}

void Chromosome::evaluate_chromosome(Matrix<double_t> &cost_mat,
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

void Chromosome::mutate(Matrix<double_t> &cost_mat,
            std::vector<double_t> &rewards,
            double_t max_cost,
            std::mt19937 &g){
  std::tie(path, angles, cost) =
      dubins_two_opt(nodes, rho, path, angles, cost);
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
  double_t q_end[3];
  q_end[0] = nodes->at(idx_finish).first;
  q_end[1] = nodes->at(idx_finish).second;
  q_end[2] = 0.0;
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
      //TODO: fix cost by checking the cheapest angle to insert.
      double_t q_prev[3];
      q_prev[0] = nodes->at(c.path.back()).first;
      q_prev[1] = nodes->at(c.path.back()).second;
      q_prev[2] = c.angles.back();
      double_t q_next[3];
      q_next[0] = nodes->at(next_vertex).first;
      q_next[1] = nodes->at(next_vertex).second;
      q_next[2] = 0;
      std::unique_ptr<DubinsPath> dubins_path = std::make_unique<DubinsPath>();
      double_t min_angle_next, min_angle_end;
      double_t min_distance = std::numeric_limits<double_t>::max();
      for (Vector<double_t>::iterator next_angle_it = std_angles.begin();
          next_angle_it != std_angles.end(); ++next_angle_it) {
        q_next[2] = *next_angle_it;
        for (Vector<double_t>::iterator end_angle_it = std_angles.begin();
            end_angle_it != std_angles.end(); ++end_angle_it){
          double_t dist = 0.0;
          dubins_init(q_prev, q_next, rho, dubins_path.get());
          dist += dubins_path_length(dubins_path.get());
          q_end[2] = *end_angle_it;
          dubins_init(q_next, q_end, rho, dubins_path.get());
          dist += dubins_path_length(dubins_path.get());
          if(dist < min_distance){
            min_distance = dist;
            min_angle_next = q_next[2];
            min_angle_end = q_end[2];
          }
        }
      }


      double next_vertex_end_cost = path_cost + min_distance;
      if (next_vertex_end_cost < max_cost || logically_equal(next_vertex_end_cost, max_cost)) {
        q_next[2] = min_angle_next;
        dubins_init(q_prev, q_next, rho, dubins_path.get());
        path_cost += dubins_path_length(dubins_path.get());
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
  double_t q_back[3];
  q_back[0] = nodes->at(c.path.back()).first;
  q_back[1] = nodes->at(c.path.back()).second;
  q_back[2] = c.angles.back();
  std::unique_ptr<DubinsPath> dubins_path = std::make_unique<DubinsPath>();
  double_t min_angle_end;
  double_t min_distance = std::numeric_limits<double_t>::max();
  for (Vector<double_t>::iterator end_angle_it = std_angles.begin();
       end_angle_it != std_angles.end(); ++end_angle_it) {
    q_end[2] = *end_angle_it;
    dubins_init(q_back, q_end, rho, dubins_path.get());
    double_t dist = dubins_path_length(dubins_path.get());
    if (dist < min_distance){
      min_distance = dist;
      min_angle_end = q_end[2];
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
                Matrix<double> &cost_mat,
                std::vector<double> &rewards,
                double &max_cost){
  // Get hash of thread id for the seed of the generator.
  std::hash<std::thread::id> hasher;
  static thread_local std::mt19937 g(hasher(std::this_thread::get_id()));

  for(size_t idx:indices)
    pop[idx].mutate(cost_mat, rewards, max_cost, g);
}

Chromosome ga_dcop(std::shared_ptr<Vector<Point2D>> nodes,
                   Vector<double_t> std_angles,
                   double_t rho,
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
    Chromosome c = generate_chromosome(
        nodes, std_angles, rho, max_cost,
        idx_start, idx_finish, cost_mat, g);
    std::tie(c.path, c.angles, c.cost) =
        dubins_two_opt(nodes, rho, c.path, c.angles, c.cost);
    c.evaluate_chromosome(cost_mat, rewards);
    pop.push_back(c);
  }

  for (int gen = 0; gen < max_gen; ++gen) {
    // Select new population
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
  return best;
}

}