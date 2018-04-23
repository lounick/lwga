//
// Created by nick on 21/07/17.
//
#include <fstream>
#include <include/top_ga.h>

int main(int argc, char *argv[]) {
  std::string filename("/home/nick/ClionProjects/LWGA/test_instances/Set_100_234/p4.2.a.txt");
  std::ifstream infile(filename);

  std::vector<std::pair<double_t, double_t> > nodes;
  std::vector<std::vector<double_t> > cost_mat;
  std::vector<double_t> rewards;
  std::vector<double_t> max_cost_v;

  std::string cmd;
  size_t num_vertices;
  size_t num_robots;
  double_t max_cost;

  infile >> cmd >> num_vertices;
  infile >> cmd >> num_robots;
  infile >> cmd >> max_cost;

  nodes.reserve(num_vertices);
  cost_mat.reserve(num_vertices);
  for(size_t i = 0; i < num_vertices; ++i){
    cost_mat.push_back(std::vector<double_t>());
  }
  rewards.reserve(num_vertices);
  max_cost_v.reserve(num_robots);

  double_t pos_x;
  double_t pos_y;
  double_t reward;
  while(infile >> pos_x >> pos_y >> reward){
    nodes.push_back(std::make_pair(pos_x,pos_y));
    rewards.push_back(reward);
  }

  for(size_t i = 0; i < num_vertices; ++i){
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < num_vertices; ++j) {
      double_t dist = find_distance(nodes[i],nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  max_cost_v = std::vector<double_t>(num_robots, max_cost);

  size_t nexp = 100;
  std::vector<double_t> fitnesses;
  fitnesses.reserve(nexp);
  std::vector<double_t> times;
  times.reserve(nexp);

  std::random_device rd;
  std::mt19937 g(rd());

  Chromosome best(cost_mat.size(), num_robots, 0, num_vertices-1);
  double best_fit = 0;

  for (size_t exp = 0; exp < nexp; ++exp) {
    auto start = std::chrono::high_resolution_clock::now();
    Chromosome c = ga_top(cost_mat, rewards, max_cost_v, 0, num_vertices - 1, g);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());

    std::unordered_set<uint_fast32_t> seen;
    for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
      for (size_t path_idx = 0; path_idx < c.genes[robot].path.size(); ++path_idx) {
        seen.insert(c.genes[robot].path[path_idx]);
      }
    }

    double fitness = 0;
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

    seen.clear();
    std::pair<std::unordered_set<uint_fast32_t>::iterator, bool> insert_ret;
    for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
      for (size_t i = 1; i < c.genes[robot].path.size() - 1; i++) {
        insert_ret = seen.insert(c.genes[robot].path[i]);
        if (insert_ret.second) {
          fitness += rewards[c.genes[robot].path[i]];
        }
      }
    }
    c.evaluate_chromosome(cost_mat, rewards, max_cost_v);
    std::cout << c.total_fitness << " " << fitness << std::endl;
    fitnesses.push_back(fitness);
    if(fitness > best_fit){
      best = c;
      best_fit = fitness;
    };
  }

  double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
  double fit_var;
  for (int i = 0; i < fitnesses.size(); ++i){
    fit_var += pow((fitnesses[i]-avg_fit),2);
  }
  fit_var /= fitnesses.size();
  double fit_stddev = sqrt(fit_var);
  std::cout << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

  double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
  double time_var;
  for (int i = 0; i < times.size(); ++i){
    time_var += pow((times[i]-avg_time),2);
  }
  time_var /= times.size();
  double time_stddev = sqrt(time_var);
  std::cout << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;

  std::cout << "best_fit: " << best_fit << " max_cost: " << max_cost_v[0] << std::endl;
  for(size_t i = 0; i < num_robots; i++){std::cout << best.genes[i].cost << std::endl;}
  std::cout << "[";
  for(size_t i = 0; i < num_robots; i++){
//    std::cout << best.genes[i].cost << " [";
    std::cout << " [";
    for(size_t j = 0; j < best.genes[i].path.size(); j++){
      if(j < best.genes[i].path.size() - 1)
        std::cout << best.genes[i].path[j] << ", ";
      else {
        if (i < num_robots - 1) {
          std::cout << best.genes[i].path[j] << "], ";
        } else {
          std::cout << best.genes[i].path[j] << "]";
        }
      }
    }
  }
  std::cout << "]" << std::endl;
  return 0;
}