//
// Created by nick on 23/07/17.
//

#include <fstream>

#include <include/ctop_ga.h>

int main(int argc, char *argv[]) {
  std::vector<std::pair<double_t, double_t> > nodes;
  std::vector<std::vector<double_t> > cost_mat;
  std::vector<double_t> max_cost_v;
  std::vector<uint8_t> num_robots;
  std::vector<std::string> grasp_methods;
  std::vector<double> rewards;
  std::vector<double> fitnesses;
  std::vector<double> mean_fitnesses;
  std::vector<double> stddev_fitnesses;
  std::vector<double> times;
  std::vector<double> mean_times;
  std::vector<double> stdev_times;
  std::random_device rd;
  std::mt19937 g(rd());

  int num_exp = 1000;
  double_t total_cost = 81.0 + 82.0;

  std::ofstream util_file;
  std::ofstream time_file;

  nodes.push_back(std::make_pair(0, 0));
  for (int i = 1; i < 10; i++) {
    for (int j = -4; j < 5; j++) {
      nodes.push_back(std::make_pair(i, j));
    }
  }
  nodes.push_back(std::make_pair(10, 0));

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

  rewards = std::vector<double>(cost_mat.size(), 1);
  rewards[0] = 0;
  rewards[82] = 0;


  grasp_methods = {std::string("GRASP"), std::string("NNGRASP")};

  num_robots = {2, 3, 4, 5, 6};
  util_file.open("/home/nick/ClionProjects/LWGA/runs/ctop_miqp_comparisson_util.txt");
  time_file.open("/home/nick/ClionProjects/LWGA/runs/ctop_miqp_comparisson_time.txt");
  util_file << "\\begin{table}[]\n";
  util_file << "\\centering\n";
  util_file << "\\caption{}\n";
  util_file << "\\label{tab:}\n";
  util_file << "\\begin{tabular}{llllll}\n";
  util_file << "\\hline \\hline\n";
  util_file << "\\multirow{2}{*}{Algorithm} & \\multicolumn{5}{c}{Number of Robots} \\Tstrut \\\\\n";
  util_file << "& 2 & 3 & 4 & 5 &6 \\Bstrut \\\\ \\hline\n";
  util_file << "COP 5\\% & & & & & \\Tstrut \\\\\n";

  time_file << "\\begin{table}[]\n";
  time_file << "\\centering\n";
  time_file << "\\caption{}\n";
  time_file << "\\label{tab:}\n";
  time_file << "\\begin{tabular}{llllll}\n";
  time_file << "\\hline \\hline\n";
  time_file << "\\multirow{2}{*}{Algorithm} & \\multicolumn{5}{c}{Number of Robots} \\Tstrut \\\\\n";
  time_file << "& 2 & 3 & 4 & 5 &6 \\Bstrut \\\\ \\hline\n";
  time_file << "COP 5\\% & & & & & \\Tstrut \\\\\n";
  for (auto &method : grasp_methods) {
    mean_fitnesses.clear();
    stddev_fitnesses.clear();
    mean_times.clear();
    stdev_times.clear();
    for (auto &nr : num_robots) {
      fitnesses.clear();
      times.clear();
      max_cost_v = std::vector<double_t>(nr, double_t(total_cost) / double_t(nr));
      for (size_t i = 0; i < num_exp; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 82, g, 100, 50, 0.75, method);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        times.push_back(diff.count());

        std::unordered_set<uint_fast32_t> seen;
        for (uint_fast32_t robot = 0; robot < nr; ++robot) {
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
        for (uint_fast32_t robot = 0; robot < nr; ++robot) {
          for (size_t i = 1; i < c.genes[robot].path.size() - 1; i++) {
            double extras = 0;
            insert_ret = seen.insert(c.genes[robot].path[i]);
            if (insert_ret.second) {
              for (size_t j = 0; j < free_vertices.size(); j++) {
                if (cost_mat[c.genes[robot].path[i]][free_vertices[j]] < 2) {
                  extras += std::exp(-2 * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
                }
              }
              fitness += rewards[c.genes[robot].path[i]] + extras;
            }
          }
        }
        fitnesses.push_back(fitness);
      }
      double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
      double fit_var;
      for (int i = 0; i < fitnesses.size(); ++i) {
        fit_var += pow((fitnesses[i] - avg_fit), 2);
      }
      fit_var /= fitnesses.size();
      double fit_stddev = sqrt(fit_var);
      mean_fitnesses.push_back(avg_fit);
      stddev_fitnesses.push_back(fit_stddev);
      double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
      double time_var;
      for (int i = 0; i < times.size(); ++i) {
        time_var += pow((times[i] - avg_time), 2);
      }
      time_var /= times.size();
      double time_stddev = sqrt(time_var);
      mean_times.push_back(avg_time);
      stdev_times.push_back(time_stddev);
    }
    util_file << "EV-" << method << " ";
    for (auto &mean_fit: mean_fitnesses) {
      util_file << "& " << mean_fit << " ";
    }
    util_file << "\\\\\n";
    util_file << "St.Dev. ";
    for (auto &stdev_fit: stddev_fitnesses) {
      util_file << "& " << stdev_fit << " ";
    }
    util_file << "\\\\\n";

    time_file << "EV-" << method << " ";
    for (auto &mean_time: mean_times) {
      time_file << "& " << mean_time << " ";
    }
    time_file << "\\\\\n";
    time_file << "St.Dev. ";
    for (auto &stdev_time: stdev_times) {
      time_file << "& " << stdev_time << " ";
    }
    time_file << "\\\\\n";
  }
  util_file << "\\hline\n";
  util_file << "\\end{tabular}\n";
  util_file << "\\end{table}\n";
  util_file << "\n";
  util_file << "Limited energy experiments";
  util_file << "\n";
  time_file << "\\hline\n";
  time_file << "\\end{tabular}\n";
  time_file << "\\end{table}\n";
  time_file << "\n";
  time_file << "Limited energy experiments";
  time_file << "\n";


  util_file << "\\begin{table}[]\n";
  util_file << "\\centering\n";
  util_file << "\\caption{}\n";
  util_file << "\\label{tab:}\n";
  util_file << "\\begin{tabular}{lllll}\n";
  util_file << "\\hline \\hline\n";
  util_file << "\\multirow{2}{*}{Algorithm} & \\multicolumn{4}{c}{Budget} \\Tstrut \\\\\n";
  util_file << "& 100\\% & 75\\%  & 50\\%  & 25\\% \\Bstrut \\\\ \\hline\n";
  util_file << "COP 5\\% & & & & \\Tstrut \\\\\n";

  time_file << "\\begin{table}[]\n";
  time_file << "\\centering\n";
  time_file << "\\caption{}\n";
  time_file << "\\label{tab:}\n";
  time_file << "\\begin{tabular}{lllll}\n";
  time_file << "\\hline \\hline\n";
  time_file << "\\multirow{2}{*}{Algorithm} & \\multicolumn{5}{c}{Budget} \\Tstrut \\\\\n";
  time_file << "& 100\\% & 75\\%  & 50\\%  & 25\\% \\Bstrut \\\\ \\hline\n";
  time_file << "COP 5\\% & & & & \\Tstrut \\\\\n";

  std::vector<double_t> energy_coef = {1, 0.75, 0.5, 0.25};
  size_t nr = 3;
  for (auto &method : grasp_methods) {
    mean_fitnesses.clear();
    stddev_fitnesses.clear();
    mean_times.clear();
    stdev_times.clear();
    for (auto &ec : energy_coef) {
      fitnesses.clear();
      times.clear();
      max_cost_v = std::vector<double_t>(nr, ec * double_t(total_cost) / double_t(nr));
      for (size_t i = 0; i < num_exp; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 82, g, 100, 50, 0.75, method);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        times.push_back(diff.count());

        std::unordered_set<uint_fast32_t> seen;
        for (uint_fast32_t robot = 0; robot < 3; ++robot) {
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
        for (uint_fast32_t robot = 0; robot < nr; ++robot) {
          for (size_t i = 1; i < c.genes[robot].path.size() - 1; i++) {
            double extras = 0;
            insert_ret = seen.insert(c.genes[robot].path[i]);
            if (insert_ret.second) {
              for (size_t j = 0; j < free_vertices.size(); j++) {
                if (cost_mat[c.genes[robot].path[i]][free_vertices[j]] < 2) {
                  extras += std::exp(-2 * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
                }
              }
              fitness += rewards[c.genes[robot].path[i]] + extras;
            }
          }
        }
        fitnesses.push_back(fitness);
      }
      double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
      double fit_var;
      for (int i = 0; i < fitnesses.size(); ++i) {
        fit_var += pow((fitnesses[i] - avg_fit), 2);
      }
      fit_var /= fitnesses.size();
      double fit_stddev = sqrt(fit_var);
      mean_fitnesses.push_back(avg_fit);
      stddev_fitnesses.push_back(fit_stddev);
      double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
      double time_var;
      for (int i = 0; i < times.size(); ++i) {
        time_var += pow((times[i] - avg_time), 2);
      }
      time_var /= times.size();
      double time_stddev = sqrt(time_var);
      mean_times.push_back(avg_time);
      stdev_times.push_back(time_stddev);
    }
    util_file << "EV-" << method << " ";
    for (auto &mean_fit: mean_fitnesses) {
      util_file << "& " << mean_fit << " ";
    }
    util_file << "\\\n";
    util_file << "St.Dev. ";
    for (auto &stdev_fit: stddev_fitnesses) {
      util_file << "& " << stdev_fit << " ";
    }
    util_file << "\\\n";

    time_file << "EV-" << method << " ";
    for (auto &mean_time: mean_times) {
      time_file << "& " << mean_time << " ";
    }
    time_file << "\\\n";
    time_file << "St.Dev. ";
    for (auto &stdev_time: stdev_times) {
      time_file << "& " << stdev_time << " ";
    }
    time_file << "\\\n";
  }
  util_file << "\\hline\n";
  util_file << "\\end{tabular}\n";
  util_file << "\\end{table}\n";
  util_file << "\n";

  time_file << "\\hline\n";
  time_file << "\\end{tabular}\n";
  time_file << "\\end{table}\n";
  time_file << "\n";

  util_file.close();
  time_file.close();
//======================================================================================================================
  util_file.open("/home/nick/ClionProjects/LWGA/runs/ctop_qualitative_comparisson_util.txt");
  time_file.open("/home/nick/ClionProjects/LWGA/runs/ctop_qualitative_comparisson_time.txt");

  util_file << "\\begin{table}[]\n";
  util_file << "\\centering\n";
  util_file << "\\caption{}\n";
  util_file << "\\label{tab:}\n";
  util_file << "\\begin{tabular}{lllll}\n";
  util_file << "\\hline \\hline\n";
  util_file << "\\multirow{2}{*}{Algorithm} & \\multicolumn{4}{c}{Number of Chromosomes} \\Tstrut \\\\\n";
  util_file << "& 1 & 10 & 100 & 1000 \\Bstrut \\\\ \\hline\n";

  time_file << "\\begin{table}[]\n";
  time_file << "\\centering\n";
  time_file << "\\caption{}\n";
  time_file << "\\label{tab:}\n";
  time_file << "\\begin{tabular}{lllll}\n";
  time_file << "\\hline \\hline\n";
  time_file << "\\multirow{2}{*}{Algorithm} & \\multicolumn{5}{c}{Number of Chromosomes} \\Tstrut \\\\\n";
  time_file << "& 1 & 10 & 100 & 1000 \\Bstrut \\\\ \\hline\n";

  max_cost_v = std::vector<double_t>(nr, double_t(total_cost) / double_t(nr));
  std::vector<uint_fast32_t> pop_sizes = {1,10,100,1000};
  for (auto &method : grasp_methods) {
    mean_fitnesses.clear();
    stddev_fitnesses.clear();
    mean_times.clear();
    stdev_times.clear();
    for (auto &pop_size : pop_sizes) {
      fitnesses.clear();
      times.clear();
      for (size_t i = 0; i < num_exp; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 82, g, pop_size, 0, 0, method);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        times.push_back(diff.count());

        std::unordered_set<uint_fast32_t> seen;
        for (uint_fast32_t robot = 0; robot < nr; ++robot) {
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
        for (uint_fast32_t robot = 0; robot < nr; ++robot) {
          for (size_t i = 1; i < c.genes[robot].path.size() - 1; i++) {
            double extras = 0;
            insert_ret = seen.insert(c.genes[robot].path[i]);
            if (insert_ret.second) {
              for (size_t j = 0; j < free_vertices.size(); j++) {
                if (cost_mat[c.genes[robot].path[i]][free_vertices[j]] < 2) {
                  extras += std::exp(-2 * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
                }
              }
              fitness += rewards[c.genes[robot].path[i]] + extras;
            }
          }
        }
        fitnesses.push_back(fitness);
      }
      double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
      double fit_var;
      for (int i = 0; i < fitnesses.size(); ++i) {
        fit_var += pow((fitnesses[i] - avg_fit), 2);
      }
      fit_var /= fitnesses.size();
      double fit_stddev = sqrt(fit_var);
      mean_fitnesses.push_back(avg_fit);
      stddev_fitnesses.push_back(fit_stddev);
      double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
      double time_var;
      for (int i = 0; i < times.size(); ++i) {
        time_var += pow((times[i] - avg_time), 2);
      }
      time_var /= times.size();
      double time_stddev = sqrt(time_var);
      mean_times.push_back(avg_time);
      stdev_times.push_back(time_stddev);
    }
    util_file << "EV-" << method << " ";
    for (auto &mean_fit: mean_fitnesses) {
      util_file << "& " << mean_fit << " ";
    }
    util_file << "\\\\\n";
    util_file << "St.Dev. ";
    for (auto &stdev_fit: stddev_fitnesses) {
      util_file << "& " << stdev_fit << " ";
    }
    util_file << "\\\\\n";

    time_file << "EV-" << method << " ";
    for (auto &mean_time: mean_times) {
      time_file << "& " << mean_time << " ";
    }
    time_file << "\\\\\n";
    time_file << "St.Dev. ";
    for (auto &stdev_time: stdev_times) {
      time_file << "& " << stdev_time << " ";
    }
    time_file << "\\\\\n";
    util_file.flush();
    time_file.flush();
  }
  util_file << "\\hline\n";
  util_file << "\\end{tabular}\n";
  util_file << "\\end{table}\n";
  util_file << "\n";
  util_file << "Local search effect\n";
  util_file << "\n";
  time_file << "\\hline\n";
  time_file << "\\end{tabular}\n";
  time_file << "\\end{table}\n";
  time_file << "\n";
  time_file << "Local search effect\n";
  time_file << "\n";

  std::vector<uint16_t> ngens = {1,10,100,1000};
  std::vector<double_t> mut_rates = {0.25, 0.5, 0.75, 1.0};

  util_file << "\\begin{table}[]\n";
  util_file << "\\centering\n";
  util_file << "\\caption{}\n";
  util_file << "\\label{tab:}\n";
  util_file << "\\begin{tabular}{cllll}\n";
  util_file << "\\hline \\hline\n";
  util_file << "\\multirow{2}{*}{\\# Gen.} &\\multirow{2}{*}{Algortithm} & \\multicolumn{4}{c}{Mutation Rate} \\Tstrut \\\\\n";
  util_file << "& & 25\\% & 50\\% & 75\\% & 100\\% \\Bstrut \\\\ \\hline\n";

  time_file << "\\begin{table}[]\n";
  time_file << "\\centering\n";
  time_file << "\\caption{}\n";
  time_file << "\\label{tab:}\n";
  time_file << "\\begin{tabular}{clllll}\n";
  time_file << "\\hline \\hline\n";
  time_file << "\\multirow{2}{*}{Number of Generations} &\multirow{2}{*}{Algortithm} & \\multicolumn{4}{c}{Mutation Rate} \\Tstrut \\\\\n";
  time_file << "& & 25\\% & 50\\% & 75\\% & 100\\% \\Bstrut \\\\ \\hline\n";

  uint16_t pop_size = 100;
  for (auto &ngen:ngens){
    util_file << "\\multirow{4}{*}{" << ngen << "}";
    time_file << "\\multirow{4}{*}{" << ngen << "}";
    for(auto &method : grasp_methods){
      mean_fitnesses.clear();
      stddev_fitnesses.clear();
      mean_times.clear();
      stdev_times.clear();
      for (auto &mut_rate:mut_rates){
        fitnesses.clear();
        times.clear();
        for (size_t i = 0; i < num_exp; ++i) {
          auto start = std::chrono::high_resolution_clock::now();
          Chromosome c = ga_ctop(cost_mat, rewards, max_cost_v, 0, 82, g, pop_size, ngen, mut_rate, method);
          auto end = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> diff = end - start;
          times.push_back(diff.count());

          std::unordered_set<uint_fast32_t> seen;
          for (uint_fast32_t robot = 0; robot < nr; ++robot) {
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
          for (uint_fast32_t robot = 0; robot < nr; ++robot) {
            for (size_t i = 1; i < c.genes[robot].path.size() - 1; i++) {
              double extras = 0;
              insert_ret = seen.insert(c.genes[robot].path[i]);
              if (insert_ret.second) {
                for (size_t j = 0; j < free_vertices.size(); j++) {
                  if (cost_mat[c.genes[robot].path[i]][free_vertices[j]] < 2) {
                    extras += std::exp(-2 * cost_mat[c.genes[robot].path[i]][free_vertices[j]]);
                  }
                }
                fitness += rewards[c.genes[robot].path[i]] + extras;
              }
            }
          }
          fitnesses.push_back(fitness);
        }
        double avg_fit = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / fitnesses.size();
        double fit_var;
        for (int i = 0; i < fitnesses.size(); ++i) {
          fit_var += pow((fitnesses[i] - avg_fit), 2);
        }
        fit_var /= fitnesses.size();
        double fit_stddev = sqrt(fit_var);
        mean_fitnesses.push_back(avg_fit);
        stddev_fitnesses.push_back(fit_stddev);
        double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        double time_var;
        for (int i = 0; i < times.size(); ++i) {
          time_var += pow((times[i] - avg_time), 2);
        }
        time_var /= times.size();
        double time_stddev = sqrt(time_var);
        mean_times.push_back(avg_time);
        stdev_times.push_back(time_stddev);
      }
      util_file << "&EV-" << method << " ";
      for (auto &mean_fit: mean_fitnesses) {
        util_file << "& " << mean_fit << " ";
      }
      util_file << "\\\\\n";
      util_file << "&St.Dev. ";
      for (auto &stdev_fit: stddev_fitnesses) {
        util_file << "& " << stdev_fit << " ";
      }
      util_file << "\\\\\n";

      time_file << "&EV-" << method << " ";
      for (auto &mean_time: mean_times) {
        time_file << "& " << mean_time << " ";
      }
      time_file << "\\\\\n";
      time_file << "&St.Dev. ";
      for (auto &stdev_time: stdev_times) {
        time_file << "& " << stdev_time << " ";
      }
      time_file << "\\\\\n";
      util_file.flush();
      time_file.flush();
    }
    util_file << "\\hline\n";
    time_file << "\\hline\n";
  }

//  util_file << "\\hline\n";
  util_file << "\\end{tabular}\n";
  util_file << "\\end{table}\n";

  time_file << "\\hline\n";
  time_file << "\\end{tabular}\n";
  time_file << "\\end{table}\n";

  util_file.close();
  time_file.close();

  return 0;
}