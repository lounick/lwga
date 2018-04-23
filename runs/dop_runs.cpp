//
// Created by nick on 12/04/18.
//

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <dcop_ga.h>

void ReadProblem(
    const std::string &filename,
    double_t &budget,
    size_t &paths,
    const std::shared_ptr<Vector<Point2D>> &nodes,
    Vector<double_t> &rewards){
  std::ifstream problem;
  problem.open(filename);
  std::string line;
  line.reserve(100);
  std::stringstream tokens;

  // Read problem spec (budget and tours);
  std::getline(problem, line);
  tokens.str(line);
  tokens >> budget >> paths;

  while (std::getline(problem, line)) {
    tokens.clear();
    tokens.str(line);
    double_t x, y, reward;
    tokens >> x >> y >> reward;
    nodes->push_back(std::make_pair(x, y));
    rewards.push_back(reward);
  }
  problem.close();
}

int main(int argc, char *argv[]) {
  const std::string folder = "/home/nick/Documents/iros18/test-set/";
  const std::vector<std::string> sets = {
//      "Tsiligirides_1/",
//      "Tsiligirides_2/",
      "Tsiligirides_3/",
//      "set_64_1/",
      "set_66_1/"
  };
  const std::map<std::string, std::string> filenames = {
      {"Tsiligirides_1/", "tsiligirides_problem_1_budget_"},
      {"Tsiligirides_2/", "tsiligirides_problem_2_budget_"},
      {"Tsiligirides_3/", "tsiligirides_problem_3_budget_"},
      {"set_64_1/", "set_64_1_"},
      {"set_66_1/", "set_66_1_"}
  };
  const std::map<std::string, std::vector<std::string>> budgets = {
      {"Tsiligirides_1/", {"05","10","15","20","25","30","35","40","46","50","55","60","65","70","73","75","80","85"}},
      {"Tsiligirides_2/", {"15","20","23","25","27","30","32","35","38","40","45"}},
      {"Tsiligirides_3/", {"015","020","025","030","035","040","045","050","055","060","065","070","075","080","085","090","095","100","105","110"}},
//      {"Tsiligirides_3/", {"105"}},
      {"set_64_1/", {"15","20","25","30","35","40","45","50","55","60","65","70","75","80"}},
      {"set_66_1/", {"005","010","015","020","025","030","035","040","045","050","055","060","065","070","075","080","085","090","095","100","105","110","115","120","125","130"}}
  };
  const std::string suffix = ".txt";
  const std::vector<double_t> rhos = {0.0000000001, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3};
  const std::string output_folder = "/home/nick/Documents/iros18/test-set/results/";

  std::random_device rd;
  std::mt19937 g(rd());

  for (const std::string &set : sets){
    std::string output_filename = output_folder + filenames.at(set) + "new" + suffix;
    std::ofstream output;
    output.open(output_filename);
//    std::ostream& output = std::cout;
    for (std::string str_budget : budgets.at(set)) {
      std::string input_filename = folder + set + filenames.at(set) + str_budget + suffix;
      std::cout << input_filename << std::endl;
      std::flush(std::cout);
      double_t budget;
      size_t paths;
      std::shared_ptr<Vector<Point2D>> nodes = std::make_shared<Vector<Point2D>>();
      Vector<double_t> rewards;
      ReadProblem(input_filename, budget, paths, nodes, rewards);

      Vector<double_t> std_angles;
      uint_fast32_t degrees = 15;
      std_angles.reserve(360/degrees);
      for(uint_fast8_t i = 0; i < 360/degrees; ++i){
        std_angles.push_back(i*degrees*M_PI/180.0);
      }

      Matrix <double_t> eucledian_cost_mat;
      for (size_t i = 0; i < nodes->size(); i++) {
        std::vector<double_t > tmp_vec;
        tmp_vec.reserve(nodes->size());
        eucledian_cost_mat.push_back(tmp_vec);
      }

      for (size_t i = 0; i < nodes->size(); i++) {
        eucledian_cost_mat[i].push_back(0);
        for (size_t j = i + 1; j < nodes->size(); j++) {
          double dist = find_distance(nodes->at(i), nodes->at(j));
          eucledian_cost_mat[i].push_back(dist);
          eucledian_cost_mat[j].push_back(dist);
        }
      }

      output << budget << " & ";

      for (double_t rho : rhos) {
        Matrix<Matrix<double_t>> dubins_cost_mat;
        dubins_cost_mat.reserve(nodes->size());
        for (size_t i = 0; i < nodes->size(); ++i) {
          Vector<Matrix<double_t>> tmp_vec;
          tmp_vec.reserve(nodes->size());
          for(size_t j = 0; j < nodes->size(); ++j){
            Matrix<double_t> angle_mat;
            angle_mat.reserve(std_angles.size());
            for (size_t k = 0; k < std_angles.size(); ++k){
              Vector<double_t> angle_vec;
              angle_vec.reserve(std_angles.size());
              angle_mat.push_back(angle_vec);
            }
            tmp_vec.push_back(angle_mat);
          }
          dubins_cost_mat.push_back(tmp_vec);
        }

        for (size_t i = 0; i < nodes->size(); i++) {
          double_t q0[3] = {nodes->at(i).first, nodes->at(i).second, 0};
          std::unique_ptr<DubinsPath> ptp_path = std::make_unique<DubinsPath>();
          for (size_t j = 0; j < nodes->size(); j++) {
            double_t q1[3] = {nodes->at(j).first, nodes->at(j).second, 0};
            for (size_t k = 0; k < std_angles.size(); ++k){
              q0[2] = std_angles[k];
              for (size_t l = 0; l < std_angles.size(); ++l){
                q1[2] = std_angles[l];
                int ret = dubins_init(q0, q1, rho, ptp_path.get());
                if (ret != 0)
                  std::cout << "Dubins ret: " << ret << std::endl;
                dubins_cost_mat[i][j][k].push_back(dubins_path_length(ptp_path.get()));
              }
            }
          }
        }
        int nexp = 10;
        dcop_ga::Chromosome best;
        double_t best_fitness = -1;
        double_t best_time = 1e9;
        for (int exp = 0; exp < nexp; exp++) {
          auto start = std::chrono::high_resolution_clock::now();
          dcop_ga::Chromosome c = dcop_ga::ga_dcop(
              nodes, std_angles, rho, dubins_cost_mat,
              eucledian_cost_mat, rewards, budget, 0, 1, 475, 35, 3, 0.6, 0.8, 0.07, g);

          auto end = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> diff = end - start;

          double fitness = 0;

          for (size_t i = 0; i < c.path.size() ; ++i) {
            fitness += rewards[c.path[i]];
          }

          if (fitness > best_fitness) {
            best = c;
            best_fitness = fitness;
            best_time = diff.count();
          }
        }
        output << best_fitness << " (" << best_time << ") & ";
      }
      output << "\n";
      std::flush(output);
    }
    output.close();
  }

  return 0;
}