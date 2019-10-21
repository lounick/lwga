#include "pop_ga.h"

int main(int argc, char *argv[]) {
  pop_ga::Properties properties;
  properties.generation_method = pop_ga::GenerationMethod::GRASP;
  properties.start_id = 0;
  properties.end_id = 10;
  properties.maximum_cost = 10;
  properties.cost_per_time_unit = 1;
  properties.grasp_greediness = 0.5;
  properties.grasp_estimated_reward = true;
  properties.max_generations = 50;
  properties.max_stable_generations = 10;
  properties.population_size = 100;
  properties.mutate_add_prob = 0.5;
  properties.num_mutations = 10;
  properties.elite_rate = 0.05;
  properties.tournament_size = 3;
  properties.cx_rate = 0.3;
  properties.mutation_rate = 0.7;
  properties.num_threads = 6;
  return 0;
}
