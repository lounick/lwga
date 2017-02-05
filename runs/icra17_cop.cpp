//
// Created by nick on 02/09/16.
//

#include "cop_ga.h"
#include <fstream>

int main(int argc, char* argv[]){
  std::vector<std::pair<double, double> > nodes;
  std::vector<std::vector<double> > cost_mat;
  std::vector<double> max_costs;
  std::vector<double> rewards;
  std::vector<double> fitnesses;
  std::vector<double> times;
  std::random_device rd;
  std::mt19937 g(rd());

  int num_grids = 5;
  int num_exp = 1000;
  double total_cost;

  std::ofstream res_file;
  res_file.open("/home/nick/ClionProjects/LWGA/runs/cop_results.txt");

  nodes.push_back(std::make_pair(0, 0));
  for (int i = 1; i < 6; ++i) {
    for (int j = -2; j < 3; ++j) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(6, 0));

  for (int i = 0; i < 27; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 27; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 27; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  rewards.assign(cost_mat.size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;

  max_costs.clear();
  total_cost = 26+25;
  max_costs.push_back(total_cost);
  max_costs.push_back(total_cost*3/4);
  max_costs.push_back(total_cost/2);
  max_costs.push_back(total_cost/4);
  res_file << "===5x5===" << std::endl;
  for(int mc = 0; mc < max_costs.size(); ++mc) {
    fitnesses.clear();
    times.clear();
    res_file << "===" << max_costs[mc] << "===" << std::endl;
    for (int i = 0; i < num_exp; ++i) {
      auto start = std::chrono::high_resolution_clock::now();
      Chromosome c = ga_cop(cost_mat, rewards, max_costs[mc], 0, 26, g);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      times.push_back(diff.count());
      double fitness = calculate_fitness(c.path, cost_mat, rewards);
      fitnesses.push_back(fitness);
    }
    double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
    double fit_var;
    for (int i = 0; i < fitnesses.size(); ++i){
      fit_var += pow((fitnesses[i]-avg_fit),2);
    }
    fit_var /= fitnesses.size();
    double fit_stddev = sqrt(fit_var);
    res_file << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

    double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double time_var;
    for (int i = 0; i < times.size(); ++i){
      time_var += pow((times[i]-avg_time),2);
    }
    time_var /= times.size();
    double time_stddev = sqrt(time_var);
    res_file << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------------

  nodes.clear();
  nodes.push_back(std::make_pair(0, -0.5));
  for (int i = 1; i < 7; ++i) {
    for (int j = -3; j < 3; ++j) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(7, -0.5));

  cost_mat.clear();
  for (int i = 0; i < 38; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 38; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 38; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  rewards.clear();
  rewards.assign(cost_mat.size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;

  max_costs.clear();
  total_cost = 37.236+36;
  max_costs.push_back(total_cost);
  max_costs.push_back(total_cost*3/4);
  max_costs.push_back(total_cost/2);
  max_costs.push_back(total_cost/4);

  res_file << "===6x6===" << std::endl;
  for(int mc = 0; mc < max_costs.size(); ++mc) {
    fitnesses.clear();
    times.clear();
    res_file << "===" << max_costs[mc] << "===" << std::endl;
    for (int i = 0; i < num_exp; ++i) {
      auto start = std::chrono::high_resolution_clock::now();
      Chromosome c = ga_cop(cost_mat, rewards, max_costs[mc], 0, 37, g);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      times.push_back(diff.count());
      double fitness = calculate_fitness(c.path, cost_mat, rewards);
      fitnesses.push_back(fitness);
    }
    double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
    double fit_var;
    for (int i = 0; i < fitnesses.size(); ++i){
      fit_var += pow((fitnesses[i]-avg_fit),2);
    }
    fit_var /= fitnesses.size();
    double fit_stddev = sqrt(fit_var);
    res_file << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

    double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double time_var;
    for (int i = 0; i < times.size(); ++i){
      time_var += pow((times[i]-avg_time),2);
    }
    time_var /= times.size();
    double time_stddev = sqrt(time_var);
    res_file << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------------

  nodes.clear();
  nodes.push_back(std::make_pair(0, 0));
  for (int i = 1; i < 8; ++i) {
    for (int j = -3; j < 4; ++j) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(8, 0));

  cost_mat.clear();
  for (int i = 0; i < 51; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 51; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 51; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  rewards.clear();
  rewards.assign(cost_mat.size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;

  max_costs.clear();
  total_cost = 50.828+49;
  max_costs.push_back(total_cost);
  max_costs.push_back(total_cost*3/4);
  max_costs.push_back(total_cost/2);
  max_costs.push_back(total_cost/4);

  res_file << "===7x7===" << std::endl;
  for(int mc = 0; mc < max_costs.size(); ++mc) {
    fitnesses.clear();
    times.clear();
    res_file << "===" << max_costs[mc] << "===" << std::endl;
    for (int i = 0; i < num_exp; ++i) {
      auto start = std::chrono::high_resolution_clock::now();
      Chromosome c = ga_cop(cost_mat, rewards, max_costs[mc], 0, 50, g);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      times.push_back(diff.count());
      double fitness = calculate_fitness(c.path, cost_mat, rewards);
      fitnesses.push_back(fitness);
    }
    double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
    double fit_var;
    for (int i = 0; i < fitnesses.size(); ++i){
      fit_var += pow((fitnesses[i]-avg_fit),2);
    }
    fit_var /= fitnesses.size();
    double fit_stddev = sqrt(fit_var);
    res_file << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

    double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double time_var;
    for (int i = 0; i < times.size(); ++i){
      time_var += pow((times[i]-avg_time),2);
    }
    time_var /= times.size();
    double time_stddev = sqrt(time_var);
    res_file << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------------

  nodes.clear();
  nodes.push_back(std::make_pair(0, -0.5));
  for (int i = 1; i < 9; ++i) {
    for (int j = -4; j < 4; ++j) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(9, -0.5));

  cost_mat.clear();
  for (int i = 0; i < 66; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 66; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 66; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  rewards.clear();
  rewards.assign(cost_mat.size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;

  max_costs.clear();
  total_cost = 65.236+64;
  max_costs.push_back(total_cost);
  max_costs.push_back(total_cost*3/4);
  max_costs.push_back(total_cost/2);
  max_costs.push_back(total_cost/4);

  res_file << "===8x8===" << std::endl;
  for(int mc = 0; mc < max_costs.size(); ++mc) {
    fitnesses.clear();
    times.clear();
    res_file << "===" << max_costs[mc] << "===" << std::endl;
    for (int i = 0; i < num_exp; ++i) {
      auto start = std::chrono::high_resolution_clock::now();
      Chromosome c = ga_cop(cost_mat, rewards, max_costs[mc], 0, 65, g);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      times.push_back(diff.count());
      double fitness = calculate_fitness(c.path, cost_mat, rewards);
      fitnesses.push_back(fitness);
    }
    double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
    double fit_var;
    for (int i = 0; i < fitnesses.size(); ++i){
      fit_var += pow((fitnesses[i]-avg_fit),2);
    }
    fit_var /= fitnesses.size();
    double fit_stddev = sqrt(fit_var);
    res_file << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

    double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double time_var;
    for (int i = 0; i < times.size(); ++i){
      time_var += pow((times[i]-avg_time),2);
    }
    time_var /= times.size();
    double time_stddev = sqrt(time_var);
    res_file << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------------

  nodes.clear();
  nodes.push_back(std::make_pair(0, 0));
  for (int i = 1; i < 10; ++i) {
    for (int j = -4; j < 5; ++j) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(10, 0));

  cost_mat.clear();
  for (int i = 0; i < 83; i++) {
    std::vector<double> tmp_vec;
    cost_mat.push_back(tmp_vec);
  }

  for (size_t i = 0; i < 83; i++) {
    cost_mat[i].push_back(0);
    for (size_t j = i + 1; j < 83; j++) {
      double dist = find_distance(nodes[i], nodes[j]);
      cost_mat[i].push_back(dist);
      cost_mat[j].push_back(dist);
    }
  }

  rewards.clear();
  rewards.assign(cost_mat.size(), 1);
  rewards.front() = 0;
  rewards.back() = 0;

  max_costs.clear();
  total_cost = 82.000+81;
  max_costs.push_back(total_cost);
  max_costs.push_back(total_cost*3/4);
  max_costs.push_back(total_cost/2);
  max_costs.push_back(total_cost/4);

  res_file << "===9x9===" << std::endl;
  for(int mc = 0; mc < max_costs.size(); ++mc) {
    fitnesses.clear();
    times.clear();
    res_file << "===" << max_costs[mc] << "===" << std::endl;
    for (int i = 0; i < num_exp; ++i) {
      auto start = std::chrono::high_resolution_clock::now();
      Chromosome c = ga_cop(cost_mat, rewards, max_costs[mc], 0, 82, g);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = end - start;
      times.push_back(diff.count());
      double fitness = calculate_fitness(c.path, cost_mat, rewards);
      fitnesses.push_back(fitness);
    }
    double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
    double fit_var;
    for (int i = 0; i < fitnesses.size(); ++i){
      fit_var += pow((fitnesses[i]-avg_fit),2);
    }
    fit_var /= fitnesses.size();
    double fit_stddev = sqrt(fit_var);
    res_file << "Average fitness: " << avg_fit << " Variance: " << fit_var << " StdDev: " << fit_stddev << std::endl;

    double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double time_var;
    for (int i = 0; i < times.size(); ++i){
      time_var += pow((times[i]-avg_time),2);
    }
    time_var /= times.size();
    double time_stddev = sqrt(time_var);
    res_file << "Average time: " << avg_time << " Variance: " << time_var << " StdDev: " << time_stddev << std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------------
  res_file.close();
  return 0;
}