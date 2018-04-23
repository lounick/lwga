//
// Created by nick on 18/03/18.
//

#include <ctop_ga.h>

int main(int argc, char *argv[]) {
  Vector<SamplingArea>
      random_problem = generate_random_problem(100, 100, 9, 0.1, 0);
  Vector<Point2D> problem;
  Vector<double_t> rewards;
  Point2D start = {0, 0};
  Point2D finish = {100, 100};

  problem.push_back(start);
  rewards.push_back(0.0);

  for (SamplingArea &a:random_problem) {
    Point2D area_start = {a.minX() - 1, a.origin().second};
    Vector<double_t> problem_rewards;
    a.nodes(generate_sampling_grid(a.width(),a.length(), area_start));
    problem.insert(std::end(problem), std::begin(*a.nodes()), std::end(*a.nodes()));
    problem_rewards = Vector<double_t>(a.nodes()->size(),1);
    rewards.insert(std::end(rewards), std::begin(problem_rewards), std::end(problem_rewards));
  }

  problem.push_back(finish);
  rewards.push_back(0.0);

/*
  Point2D start = {0, 0};
  Point2D finish = {16, 0};
  SamplingArea s1({8, 5}, 5, 5);
  s1.nodes(generate_sampling_grid(s1.width(),
                                  s1.length(),
                                  {s1.minX() - 1, s1.origin().second}));
//  s1.print(true);
  SamplingArea s2({8, -5}, 5, 5);
  s2.nodes(generate_sampling_grid(s2.width(),
                                  s2.length(),
                                  {s2.minX() - 1, s2.origin().second}));
//  s2.print(true);

  Vector<Point2D> problem;
  Vector<double_t> rewards;
  Vector<double_t> problem_rewards;
  problem.push_back(start);
  rewards.push_back(0.0);
  problem.insert(std::end(problem), std::begin(*s1.nodes()), std::end(*s1.nodes()));
  problem_rewards = Vector<double_t>(s1.nodes()->size(),1);
  rewards.insert(std::end(rewards), std::begin(problem_rewards), std::end(problem_rewards));
  problem.insert(std::end(problem), std::begin(*s2.nodes()), std::end(*s2.nodes()));
  problem_rewards = Vector<double_t>(s2.nodes()->size(),1);
  rewards.insert(std::end(rewards), std::begin(problem_rewards), std::end(problem_rewards));
  problem.push_back(finish);
  rewards.push_back(0.0);
  */


  std::vector<std::vector<double> > cost_mat;
  for (int i = 0; i < problem.size(); i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < problem.size(); i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < problem.size(); j++) {
      double dist = find_distance(problem[i], problem[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  uint_fast32_t num_robots = 2;
//  std::vector<double> max_cost_v(num_robots, 14 + 24);
  std::vector<double> max_cost_v(num_robots, 500);
  std::vector<double> fitnesses;
  std::vector<double> times;

  int nexp = 10;

  std::random_device rd;
  std::mt19937 g(rd());

  ctop_ga::Chromosome best(cost_mat.size(), num_robots, 0, problem.size() - 1);
  double best_fit = 0;

  for (int exp = 0; exp < nexp; exp++) {
    auto start = std::chrono::high_resolution_clock::now();
//    Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 26, g);
    ctop_ga::Chromosome c = ctop_ga::ga_ctop(
//        cost_mat, rewards, max_cost_v, 0, 82, g, 250, 50, 8, 0.9, 0.8,
        cost_mat,
        rewards,
        max_cost_v,
        0,
        problem.size() - 1,
        g,
        250,
        50,
        5,
        0.9,
        0.7,
        "NNGRASP",
        0.03);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    times.push_back(diff.count());

    c.evaluate_chromosome(cost_mat, rewards, max_cost_v);

    std::unordered_set<uint_fast32_t> seen;
    for (uint_fast32_t robot = 0; robot < num_robots; ++robot) {
      for (size_t path_idx = 0; path_idx < c.genes[robot].path.size();
           ++path_idx) {
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
        double extras = 0;
        insert_ret = seen.insert(c.genes[robot].path[i]);
        if (insert_ret.second) {
          for (size_t j = 0; j < free_vertices.size(); j++) {
            if (cost_mat[c.genes[robot].path[i]][free_vertices[j]] < 2) {
              extras += std::exp((log(0.01) / 2)
                                     * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
            }
          }
          fitness += rewards[c.genes[robot].path[i]] + extras;
        }
      }
    }
    std::cout << exp << " " << c.total_fitness << " " << fitness << std::endl;
    fitnesses.push_back(fitness);
    if (fitness > best_fit) {
      best = c;
      best_fit = fitness;
    };
  }
  double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0)
      / fitnesses.size();
  double fit_var;
  for (int i = 0; i < fitnesses.size(); ++i) {
    fit_var += pow((fitnesses[i] - avg_fit), 2);
  }
  fit_var /= fitnesses.size();
  double fit_stddev = sqrt(fit_var);
  std::cout << "Average fitness: " << avg_fit << " Variance: " << fit_var
            << " StdDev: " << fit_stddev << std::endl;

  double avg_time =
      std::accumulate(times.begin(), times.end(), 0.0) / times.size();
  double time_var;
  for (int i = 0; i < times.size(); ++i) {
    time_var += pow((times[i] - avg_time), 2);
  }
  time_var /= times.size();
  double time_stddev = sqrt(time_var);
  std::cout << "Average time: " << avg_time << " Variance: " << time_var
            << " StdDev: " << time_stddev << std::endl;

  std::cout << "best_fit: " << best_fit << " max_cost: " << max_cost_v[0]
            << std::endl;
  for (size_t i = 0; i < num_robots; i++) {
    std::cout << best.genes[i].cost << " " << get_path_cost(best.genes[i].path,
                                                            cost_mat)
              << std::endl;
  }
  std::cout << "[";
  for (size_t i = 0; i < num_robots; i++) {
//    std::cout << best.genes[i].cost << " [";
    std::cout << " [";
    for (size_t j = 0; j < best.genes[i].path.size(); j++) {
      if (j < best.genes[i].path.size() - 1)
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

  std::cout << "Total time: "
            << std::accumulate(times.begin(), times.end(), 0.0) << std::endl;
  return 0;
}