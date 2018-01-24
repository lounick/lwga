//
// Created by nick on 24/01/18.
//

#include <dcop_ga.h>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>


namespace dcop_ga {

template<typename T>
struct Limits {
  T min;
  T max;
};

struct Game {
  double_t opponent_rating;
  double_t opponent_rating_interval;
  double_t opponent_rating_deviation;
  Matrix<double_t> results;
};

using Games = std::vector<Game>;

double_t tau = 0.5;

double_t mu(double_t R){
  return (R - 1500)/173.7178;
}

double_t phi(double_t RD){
  return RD/173.7178;
}

double_t g(double_t phi){
  return 1 / std::sqrt(1 + (3 * (std::pow(phi, 2.0) / std::pow(M_PI, 2.0))));
}

double_t E(double_t mu, double_t mui, double phii){
  return 1 / (1 + pow(10, -g(phii) * (mu - mui)));
}

class Configuration {
 public:
  Configuration(uint_fast16_t pop_size,
                uint_fast8_t num_gen,
                uint_fast8_t tour_size,
                double_t cx_rate,
                double_t mut_rate,
                double_t elitist_rate,
                double_t rating_interval,
                Limits<uint_fast16_t> pop_size_lim,
                Limits<uint_fast8_t> num_gen_lim,
                Limits<uint_fast8_t> tour_size_lim,
                Limits<double_t> cx_lim,
                Limits<double_t> mut_lim,
                Limits<double_t> elites_lim)
      : pop_size_(pop_size),
        num_gen_(num_gen),
        tour_size_(tour_size),
        cx_rate_(cx_rate),
        mut_rate_(mut_rate),
        elitist_rate_(elitist_rate),
        rating_interval_(rating_interval),
        pop_size_lim_(pop_size_lim),
        num_gen_lim_(num_gen_lim),
        tour_size_lim_(tour_size_lim),
        cx_lim_(cx_lim),
        mut_lim_(mut_lim),
        elites_lim_(elites_lim) {
    rating_ = 1500;
    rating_deviation_ = 350;
    volatility_ = 0.06;
  }
  Configuration() {
    rating_ = 1500;
    rating_deviation_ = 350;
    volatility_ = 0.06;
  }
  inline double_t rating() const { return rating_; }
  inline void rating(double_t rating) { rating_ = rating; }
  inline double_t rating_deviation() const { return rating_deviation_; }
  inline void rating_deviation(double_t rating_deviation) {
    rating_deviation_ = rating_deviation;
  }
  inline double_t rating_interval() const { return rating_interval_; }
  inline void rating_interval(double_t rating_interval) {
    rating_interval_ = rating_interval;
  }
  inline uint_fast16_t pop_size() const { return pop_size_; }
  inline void pop_size(uint_fast16_t pop_size) {
    if (pop_size > pop_size_lim_.max)
      pop_size_ = pop_size_lim_.max;
    else if (pop_size < pop_size_lim_.min)
      pop_size_ = pop_size_lim_.min;
    else
      pop_size_ = pop_size;
  }
  inline uint_fast8_t num_gen() const { return num_gen_; }
  inline void num_gen(uint_fast8_t num_gen) {
    if(num_gen > num_gen_lim_.max)
      num_gen_ = num_gen_lim_.max;
    else if (num_gen < num_gen_lim_.min)
      num_gen_ = num_gen_lim_.min;
    else
      num_gen_ = num_gen;
  }
  inline uint_fast8_t tour_size() const { return tour_size_; }
  inline void tour_size(uint_fast8_t tour_size) {
    if(tour_size > tour_size_lim_.max)
      tour_size_ = tour_size_lim_.max;
    else if (tour_size < tour_size_lim_.min)
      tour_size_ = tour_size_lim_.min;
    else
      tour_size_ = tour_size;
  }
  inline double_t cx_rate() const { return cx_rate_; }
  inline void cx_rate(double_t cx_rate) {
    if (cx_rate > cx_lim_.max)
      cx_rate_ = cx_lim_.max;
    else if (cx_rate < cx_lim_.min)
      cx_rate_ = cx_lim_.min;
    else
      cx_rate_ = cx_rate;
  }
  inline double_t mut_rate() const { return mut_rate_; }
  inline void mut_rate(double_t mut_rate) {
    if (mut_rate > mut_lim_.max)
      mut_rate_ = mut_lim_.max;
    else if (mut_rate < mut_lim_.min)
      mut_rate_ = mut_lim_.min;
    else
      mut_rate_ = mut_rate;
  }
  inline double_t elitist_rate() const { return elitist_rate_; }
  inline void elitist_rate(double_t elitist_rate) {
    if (elitist_rate > elites_lim_.max)
      elitist_rate_ = elites_lim_.max;
    else if (elitist_rate < elites_lim_.min)
      elitist_rate_ = elites_lim_.min;
    else
      elitist_rate_ = elitist_rate;
  }
  inline double_t results(size_t row, size_t col) const {
    return results_[row][col];
  }
  inline void results(size_t row, size_t col, double_t val) {
    results_[row][col] = val;
  }
  inline double_t time_results(size_t row, size_t col) const {
    return time_results_[row][col];
  }
  inline void time_results(size_t row, size_t col, double_t val) {
    time_results_[row][col] = val;
  }
  void initResults(size_t num_rows, size_t num_cols) {
    results_.clear();
    results_.reserve(num_rows);
    for (size_t row = 0; row < num_rows; ++row) {
      results_.push_back(std::move(std::vector<double_t>(num_cols, 0)));
    }
    time_results_.clear();
    time_results_.reserve(num_rows);
    for (size_t row = 0; row < num_rows; ++row) {
      time_results_.push_back(std::move(std::vector<double_t>(num_cols, 0)));
    }

  }
  void initGames(size_t num_games) {
    games_.clear();
    games_.reserve(num_games);
  }
  void insertGame(Game &&game) {
    games_.push_back(std::move(game));
  }
  void updateRating() {
    double_t mu_ = mu(rating_);
    double_t phi_ = phi(rating_deviation_);
    double_t v = 0.0;
    for (size_t game = 0; game < games_.size(); ++game) {
      double_t mui = mu(games_[game].opponent_rating);
      double_t phii = phi(games_[game].opponent_rating_deviation);
      double_t E_ = E(mu_, mui, phii);
      v += std::pow(g(phii),2)* E_ * (1-E_);
    }
    v = 1/v;
    double_t perfFromGames = 0.0;
    for (size_t game = 0; game < games_.size(); ++game) {
      double_t mui = mu(games_[game].opponent_rating);
      double_t phii = phi(games_[game].opponent_rating_deviation);
      double_t E_ = E(mu_, mui, phii);
      for (size_t prob = 0; prob < games_[game].results.size(); ++prob) {
        for(size_t iter = 0; iter < games_[game].results[prob].size(); ++iter) {
          perfFromGames += g(phii) * (games_[game].results[prob][iter] - E_);
        }
      }
    }
    double_t delta = v * perfFromGames;

    double_t sigma2 = std::pow(volatility_, 2);
    double_t a = std::log(sigma2);
    double_t A = a;
    double_t B;
    double_t delta2 = std::pow(delta, 2);
    double_t phi2 = std::pow(phi_, 2);
    double_t t2 = std::pow(tau, 2);
    auto f = [delta2, phi2, v, a, t2](double_t x) {
      return std::exp(x) * ((delta2 - phi2 - v - std::exp(x)) /
          (2 * std::pow(phi2 + v + std::exp(x),2))) - (x - a)/t2;};
    if (delta2 > phi2 + v)
      B = std::log(delta2 - phi2 - v);
    else {
      uint_fast32_t k = 1;
      while(f(a - k * tau) < 0)
        ++k;
      B = a - k * tau;
    }
    double_t fa = f(A);
    double_t fb = f(B);
    double_t eps = 0.0000001;
    while (std::abs(B - A) > eps){
      double_t C = A + (A - B)*fa/(fb - fa);
      double_t fc = f(C);
      if (fc*fb < 0) {
        A = B;
        fa = fb;
      } else {
        fa = fa/2;
      }
      B = C;
      fb = fc;
    }

    double_t sigma_prime = std::exp(A/2);
    double_t phi_star = std::sqrt(phi2 + std::pow(sigma_prime, 2));
    double_t phi_prime = 1 / std::sqrt(1/std::pow(phi_star,2) + 1/v);
    double_t mu_prime = 0.0;

    mu_prime = mu_ + phi_prime*perfFromGames;
    double_t R_prime = 173.7178 * mu_prime + 1500;
    double_t RD_prime = 173.7178 * phi_prime;
    if (RD_prime < 50)
      RD_prime = 50;
    if (RD_prime > 350)
      RD_prime = 350;
    rating_ = R_prime;
    rating_deviation_ = RD_prime;
    volatility_ = sigma_prime;
  }
  Limits<double_t> getRatingInterval() {
    return {rating_ - rating_interval_ * rating_deviation_,
            rating_ + rating_interval_ * rating_deviation_};
  }
 private:
  double_t rating_;
  double_t rating_deviation_;
  double_t rating_interval_;
  double_t volatility_;
  uint_fast16_t pop_size_;
  uint_fast8_t num_gen_;
  uint_fast8_t tour_size_;
  double_t cx_rate_;
  double_t mut_rate_;
  double_t elitist_rate_;
  Matrix<double_t> results_;
  Matrix<double_t> time_results_;
  Games games_;
  Limits<uint_fast16_t> pop_size_lim_;
  Limits<uint_fast8_t> num_gen_lim_;
  Limits<uint_fast8_t> tour_size_lim_;
  Limits<double_t> cx_lim_;
  Limits<double_t> mut_lim_;
  Limits<double_t> elites_lim_;
};

using Configurations = std::vector<Configuration>;

struct Problem {
  std::shared_ptr<Vector<Point2D>> nodes;
  Matrix<Matrix<double_t>> dubins_cost_mat;
  Matrix<double_t> cost_mat;
  std::vector<double_t> rewards;
  std::vector<double_t> max_cost;
};

using Problems = std::vector<Problem>;

Configurations generatePopulation(uint_fast8_t pop_size,
                                  Limits<uint_fast16_t> pop_size_lim,
                                  Limits<uint_fast8_t> num_generations_lim,
                                  Limits<uint_fast8_t> tour_size_lim,
                                  Limits<double_t> cx_rate_lim,
                                  Limits<double_t> mut_rate_lim,
                                  Limits<double_t> elitist_rate_lim,
                                  double_t rating_interval,
                                  std::mt19937 &g) {
  std::vector<uint_fast16_t> pop_classes;
  pop_classes.reserve(20);
  uint_fast16_t pop_chunk = (pop_size_lim.max - pop_size_lim.min)/19;
  for(uint8_t cnt = 0; cnt < 20; ++cnt){
    pop_classes.push_back(pop_size_lim.min + cnt * pop_chunk);
  }
  std::vector<uint8_t> num_gen_classes;
  num_gen_classes.reserve(10);
  uint8_t num_gen_chunk = (num_generations_lim.max - num_generations_lim.min)/9;
  for(uint8_t cnt = 0; cnt < 10; ++cnt){
    num_gen_classes.push_back(num_generations_lim.min + cnt * num_gen_chunk);
  }
  std::vector<double_t> cx_classes;
  cx_classes.reserve(11);
  double_t cx_chunk = (cx_rate_lim.max - cx_rate_lim.min)/9.0;
  for(uint8_t cnt = 0; cnt < 10; ++cnt){
    cx_classes.push_back(cx_rate_lim.min + cnt * cx_chunk);
  }
  std::vector<double_t> mut_classes;
  mut_classes.reserve(11);
  double_t mut_chunk = (mut_rate_lim.max - mut_rate_lim.min)/9.0;
  for(uint8_t cnt = 0; cnt < 10; ++cnt){
    mut_classes.push_back(mut_rate_lim.min + cnt * mut_chunk);
  }

  std::vector<double_t> elite_classes;
  elite_classes.reserve(20);
  double_t elite_chunk = (elitist_rate_lim.max - elitist_rate_lim.min)/19.0;
  for (uint8_t cnt = 0; cnt < 20; ++cnt){
    elite_classes.push_back(elitist_rate_lim.min + cnt * elite_chunk);
  }

  std::uniform_int_distribution<size_t> dis10(0,9);
  std::uniform_int_distribution<size_t> dis11(0,10);
  std::uniform_int_distribution<size_t> dis20(0,19);

  std::uniform_int_distribution<uint_fast16_t>
      pop_size_dis(pop_size_lim.min, pop_size_lim.max);
  std::uniform_int_distribution<uint_fast8_t>
      num_generations_dis(num_generations_lim.min, num_generations_lim.max);
  std::uniform_int_distribution<uint_fast8_t>
      tour_size_dis(tour_size_lim.min, tour_size_lim.max);
  std::uniform_real_distribution<double_t>
      cx_rate_dis(cx_rate_lim.min, cx_rate_lim.max);
  std::uniform_real_distribution<double_t>
      mut_rate_dis(mut_rate_lim.min, mut_rate_lim.max);
  std::uniform_real_distribution<double_t>
      elitist_rate_dis(elitist_rate_lim.min, elitist_rate_lim.max);

  Configurations pop;
  pop.reserve(pop_size);
  for (uint_fast8_t i = 0; i < pop_size; ++i) {
    pop.emplace_back(pop_classes[dis20(g)],
                     num_gen_classes[dis10(g)],
                     tour_size_dis(g),
                     cx_classes[dis10(g)],
                     mut_classes[dis10(g)],
                     elite_classes[dis20(g)],
                     rating_interval,
                     pop_size_lim,
                     num_generations_lim,
                     tour_size_lim,
                     cx_rate_lim,
                     mut_rate_lim,
                     elitist_rate_lim);
  }
  return pop;
}

Problems generateProblems() {
 //TODO: Generate the problems
}

double_t evaluateChromosome(
    //TODO: Fix me to work with the correct chromosome.
    Chromosome &c, Matrix<double_t> &cost_mat, Vector<double_t> &rewards,
    Vector<double_t> &max_cost_v) {
  c.evaluate_chromosome(cost_mat, rewards, max_cost_v);
  std::unordered_set<uint_fast32_t> seen;
  for (uint_fast32_t robot = 0; robot < max_cost_v.size(); ++robot) {
    for (unsigned long path_idx : c.genes[robot].path) {
      seen.insert(path_idx);
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
  for (uint_fast32_t robot = 0; robot < max_cost_v.size(); ++robot) {
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
  return fitness;
}

void cx(
    const Configuration &parent1, const Configuration &parent2,
    Configuration &child1, Configuration &child2,
    double_t cx_prob, std::mt19937 &g) {
  std::uniform_real_distribution<double_t> prob(0.0, 1.0);
  child1 = parent1;
  child1.rating(1500);
  child1.rating_deviation(350);
  child2 = parent2;
  child2.rating(1500);
  child2.rating_deviation(350);
  if (prob(g) < cx_prob) {
    uint_fast16_t tmp = child1.pop_size();
    child1.pop_size(child2.pop_size());
    child2.pop_size(tmp);
  }
  if (prob(g) < cx_prob) {
    uint_fast8_t tmp = child1.num_gen();
    child1.num_gen(child2.num_gen());
    child2.num_gen(tmp);
  }
  if (prob(g) < cx_prob) {
    uint_fast8_t tmp = child1.tour_size();
    child1.tour_size(child2.tour_size());
    child2.tour_size(tmp);
  }
  if (prob(g) < cx_prob) {
    double_t tmp = child1.cx_rate();
    child1.cx_rate(child2.cx_rate());
    child2.cx_rate(tmp);
  }
  if (prob(g) < cx_prob) {
    double_t tmp = child1.mut_rate();
    child1.mut_rate(child2.mut_rate());
    child2.mut_rate(tmp);
  }
  if (prob(g) < cx_prob) {
    double_t tmp = child1.elitist_rate();
    child1.elitist_rate(child2.elitist_rate());
    child2.elitist_rate(tmp);
  }
}

void cx(
    const Configuration &parent1, const Configuration &parent2,
    Configuration &child, double_t cx_prob, std::mt19937 &g) {
  std::uniform_real_distribution<double_t> prob(0.0, 1.0);
  child = parent1;
  child.rating(1500);
  child.rating_deviation(350);
  if (prob(g) < cx_prob) {
    child.pop_size(parent2.pop_size());
  }
  if (prob(g) < cx_prob) {
    child.num_gen(parent2.num_gen());
  }
  if (prob(g) < cx_prob) {
    child.tour_size(parent2.tour_size());
  }
  if (prob(g) < cx_prob) {
    child.cx_rate(parent2.cx_rate());
  }
  if (prob(g) < cx_prob) {
    child.mut_rate(parent2.mut_rate());
  }
  if (prob(g) < cx_prob) {
    child.elitist_rate(parent2.elitist_rate());
  }
}

void mut(Configuration &c, double_t mut_prob, std::mt19937 &g) {
  uint8_t pop_change = 25;
  uint8_t gen_change = 5;
  uint8_t tour_change = 1;
  double_t cx_change = 0.1;
  double_t mut_change = 0.1;
  double_t elite_change = 0.01;
  std::uniform_real_distribution<double_t> prob(0.0, 1.0);
  if (prob(g) < mut_prob) {
    if (prob(g) < 0.5)
      c.pop_size(c.pop_size() + pop_change);
    else
      c.pop_size(c.pop_size() - pop_change);
  }
  if (prob(g) < mut_prob) {
    if (prob(g) < 0.5)
      c.num_gen(c.num_gen() + gen_change);
    else
      c.num_gen(c.num_gen() - gen_change);
  }
  if (prob(g) < mut_prob) {
    if (prob(g) < 0.5)
      c.tour_size(c.tour_size() + tour_change);
    else
      c.tour_size(c.tour_size() - tour_change);
  }
  if (prob(g) < mut_prob) {
    if (prob(g) < 0.5)
      c.cx_rate(c.cx_rate() + cx_change);
    else
      c.cx_rate(c.cx_rate() - cx_change);
  }
  if (prob(g) < mut_prob) {
    if (prob(g) < 0.5)
      c.mut_rate(c.mut_rate() + mut_change);
    else
      c.mut_rate(c.mut_rate() - mut_change);
  }
  if (prob(g) < mut_prob) {
    if (prob(g) < 0.5)
      c.elitist_rate(c.elitist_rate() + elite_change);
    else
      c.elitist_rate(c.elitist_rate() - elite_change);
  }
}

}

using namespace dcop_ga;

int main(int argc, char *argv[]) {
  uint_fast8_t pop_size = 100;
  Limits<uint_fast16_t> pop_size_lim = {25, 500};
  Limits<uint_fast8_t> num_gen_lim = {5, 50};
  Limits<uint_fast8_t> tour_size_lim = {3, 10};
  Limits<double_t> cx_lim = {0.0, 0.9};
  Limits<double_t> mut_lim = {0.0, 0.9};
  Limits<double_t> elites_lim = {0.01, 0.2};
  double_t rating_interval = 2.0;
  double_t cx_prob = 0.5;
  double_t mut_prob = 0.8;
  std::random_device rd;
  std::mt19937 g(rd());

  Configurations C = generatePopulation(
      pop_size, pop_size_lim, num_gen_lim, tour_size_lim, cx_lim, mut_lim,
      elites_lim, rating_interval, g);
  Problems P = generateProblems();
  uint8_t max_exp = 5;//10;
  uint8_t max_runs = 10;//25;
  double_t MPS = pop_size / 2;
  std::vector<std::string> methods = {"RANDOM", "NNGRASP"};
  for (auto &method: methods) {
    std::cout << method << std::endl;
    for (uint8_t exp = 0; exp < max_exp; ++exp) {
      for (size_t config = 0; config < C.size(); ++config) {
        C[config].initResults(P.size(), max_runs);
        std::cout << "Configuration[" << std::setfill('0') << std::setw(3)
                  << config << "]: [" << unsigned(C[config].num_gen()) << ", "
                  << unsigned(C[config].pop_size()) << ", "
                  << unsigned(C[config].tour_size()) << ", "
                  << C[config].cx_rate() << ", " << C[config].mut_rate() << ", "
                  << C[config].elitist_rate() << "]" << std::endl;
        for (size_t prob = 0; prob < P.size(); ++prob) {
          for (uint8_t run = 0; run < max_runs; ++run) {
//          std::cout << "Running configuration " << std::setfill('0')
//                    << std::setw(3) << config << ", problem "
//                    << std::setfill('0') << std::setw(3) << prob
//                    << ", run " << std::setfill('0') << std::setw(3)
//                    << unsigned(run) << std::endl;
            Chromosome c;
            double_t result;
            auto start = std::chrono::high_resolution_clock::now();
            //TODO: Fix the call to ga_dcop (both here and in the lib)
//            c = ga_dcop(
//                P[prob].cost_mat, P[prob].rewards, P[prob].max_cost_v,
//                0, P[prob].cost_mat.size() - 1, g, C[config].pop_size(),
//                C[config].num_gen(), C[config].tour_size(), C[config].cx_rate(),
//                C[config].mut_rate(), "NNGRASP", C[config].elitist_rate()
//            );
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
            //TODO: Fix call to evaluateChromosome.
//            result = evaluateChromosome(
//                c, P[prob].cost_mat, P[prob].rewards, P[prob].max_cost_v);
            C[config].results(prob, run, result);
            C[config].time_results(prob, run, diff.count());
          }
        }
      }

//    for (size_t config = 0; config < C.size(); ++config) {
//      char comma[3] = {'\0', ' ', '\0'};
//      std::cout << "Config[" << std::setfill('0')
//                << std::setw(3) << config << "]: " << std::endl;
//      for (size_t prob = 0; prob < P.size(); ++prob) {
//        std::cout << "\t" << "Problem[" << std::setfill('0')
//                  << std::setw(3) << prob << "]: [";
//        for (uint8_t run = 0; run < max_runs; ++run) {
//          std::cout << comma << C[config].results(prob, run);
//          comma[0] = ',';
//        }
//        std::cout << "]" << std::endl;
//      }
//    }

      for (size_t config = 0; config < C.size(); ++config) {
        C[config].initGames(C.size() - 1);
      }
      std::cout << "initialised games" << std::endl;
      for (size_t config1 = 0; config1 < C.size(); ++config1) {
        for (size_t config2 = config1 + 1; config2 < C.size(); ++config2) {
          Game g1, g2;
          g1.opponent_rating = C[config2].rating();
          g1.opponent_rating_deviation = C[config2].rating_deviation();
          g1.opponent_rating_interval = C[config2].rating_interval();
          g1.results.reserve(P.size());
          g2.opponent_rating = C[config1].rating();
          g2.opponent_rating_deviation = C[config1].rating_deviation();
          g2.opponent_rating_interval = C[config1].rating_interval();
          g2.results.reserve(P.size());
          for (size_t prob = 0; prob < P.size(); ++prob) {
            Vector<double_t> p_results_1, p_results_2;
            p_results_1.reserve(max_runs);
            p_results_2.reserve(max_runs);
            for (uint8_t run = 0; run < max_runs; ++run) {
              double_t p1_res = C[config1].results(prob, run);
              double_t p2_res = C[config2].results(prob, run);
              double_t p1_time = C[config1].time_results(prob, run);
              double_t p2_time = C[config2].time_results(prob, run);
              double_t timediff = p1_time - p2_time;

//            if(timediff < 0) {
//              p1_res = p1_res + 0.05 * p1_res;
//              p2_res = p2_res - 0.05 * p2_res;
//            }
//            else {
//              p1_res = p1_res - 0.05 * p1_res;
//              p2_res = p2_res + 0.05 * p2_res;
//            }

              if (logically_equal(p1_res, p2_res)) {
                p_results_1.push_back(0.5);
                p_results_2.push_back(0.5);
              } else if (p1_res < p2_res) {
                p_results_1.push_back(0.0);
                p_results_2.push_back(1.0);
              } else {
                p_results_1.push_back(1.0);
                p_results_2.push_back(0.0);
              }
              if (logically_equal(p1_time, p2_time)) {
                p_results_1.back() += 0.0;
                p_results_2.back() += 0.0;
              } else if (p1_time < p2_time) {
                p_results_1.back() += 0.1;
                p_results_2.back() += -0.1;
              } else {
                p_results_2.back() += 0.1;
                p_results_1.back() += -0.1;
              }
            }
            g1.results.push_back(std::move(p_results_1));
            g2.results.push_back(std::move(p_results_2));
          }
          C[config1].insertGame(std::move(g1));
          C[config2].insertGame(std::move(g2));
        }
      }
      std::cout << "Played games" << std::endl;
      for (size_t config = 0; config < C.size(); ++config) {
        C[config].updateRating();
      }
      std::cout << "Updated rating" << std::endl;
      auto sortRuleLambda =
          [](const Configuration &c1, const Configuration &c2) -> bool {
            return c1.rating() > c2.rating();
          };
      std::sort(C.begin(), C.end(), sortRuleLambda);
      for (Configurations::iterator it = C.begin(); it != C.end(); ++it) {
        std::cout << it->rating() << " " << it->getRatingInterval().min << " "
                  << it->getRatingInterval().max << std::endl;
      }
      Configurations parents;
      Configurations newC;
      parents.reserve(C.size());
      newC.reserve(C.size());
//    parents.push_back(*C.begin());
      for (Configurations::iterator it = C.begin(); it != C.begin() + 10;
           ++it) {
        parents.push_back(*it);
      }
      Limits<double_t> best_interval = parents.front().getRatingInterval();
      for (Configurations::iterator it = C.begin() + 1; it != C.end(); ++it) {
        if (!(best_interval.min > it->getRatingInterval().max ||
            parents.size() >= MPS)) {
          parents.push_back(*it);
        }
      }
      for (size_t config = 0; config < 10; ++config) {
        std::cout << "Parent[" << std::setfill('0') << std::setw(3)
                  << config << "]: [" << unsigned(parents[config].num_gen())
                  << ", "
                  << unsigned(parents[config].pop_size()) << ", "
                  << unsigned(parents[config].tour_size()) << ", "
                  << parents[config].cx_rate() << ", "
                  << parents[config].mut_rate() << ", "
                  << parents[config].elitist_rate() << "]" << std::endl;
      }
      newC = parents;
      while (newC.size() < pop_size) {
        std::vector<size_t> parent_indices;
        parent_indices = get_population_sample(parents.size(), 2, g);
        Configuration child;
        cx(C[parent_indices[0]], C[parent_indices[1]], child, cx_prob, g);
        mut(child, mut_prob, g);
        newC.push_back(std::move(child));
      }
      C = newC;
      for (size_t config = 0; config < 10; ++config) {
        std::cout << "Config[" << std::setfill('0') << std::setw(3)
                  << config << "]: [" << unsigned(C[config].num_gen()) << ", "
                  << unsigned(C[config].pop_size()) << ", "
                  << unsigned(C[config].tour_size()) << ", "
                  << C[config].cx_rate() << ", " << C[config].mut_rate() << ", "
                  << C[config].elitist_rate() << "]" << std::endl;
      }
      std::cout << "Generated new population" << std::endl;
    }
    auto sortRuleLambda =
        [](const Configuration &c1, const Configuration &c2) -> bool {
          return c1.rating() > c2.rating();
        };
    std::sort(C.begin(), C.end(), sortRuleLambda);
    std::cout << "Config\t" << "Rating\t" << "Pop. Size\t"
              << "Num. Generations\t" << "Tour size\t" << "CX%\t" << "MUT%\t"
              << "Elite%\t" << std::endl;
    for (size_t config = 0; config < C.size(); ++config) {
      std::cout << std::setfill('0') << std::setw(3) << config << "\t"
                << C[config].rating() << "\t" << unsigned(C[config].pop_size())
                << "\t"
                << unsigned(C[config].num_gen()) << "\t"
                << unsigned(C[config].tour_size()) << "\t"
                << C[config].cx_rate() << "\t" << C[config].mut_rate() << "\t"
                << C[config].elitist_rate() << "\t" << std::endl;
//    char comma[3] = {'\0', ' ', '\0'};
//    std::cout << "Config[" << std::setfill('0')
//              << std::setw(3) << config << "]: " << std::endl;
//    for (size_t prob = 0; prob < P.size(); ++prob) {
//      std::cout << "\t" << "Problem[" << std::setfill('0')
//                << std::setw(3) << prob << "]: [";
//      for (uint8_t run = 0; run < max_runs; ++run) {
//        std::cout << comma << C[config].results(prob, run);
//        comma[0] = ',';
//      }
//      std::cout << "]" << std::endl;
//    }
    }
  }
  return 0;
}