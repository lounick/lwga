//
// Created by nick on 16/03/18.
//

#include <mdmtsp_solver.h>
#include "gurobi_c++.h"

double_t runMDMTSP(
    const uint_fast8_t num_vehicles, const uint_fast32_t num_vertices,
    const uint_fast8_t num_depots, const Matrix<double_t> &cost_mat,
    const bool fair, Matrix<uint_fast32_t> &tours) {
  int_fast32_t L, K;
  if (fair) {
    K = std::floor((num_vertices - 1) / num_vehicles);
//    std::cout << K << std::endl;
  } else {
    K = 1;
  }

  L = num_vertices;

  GRBEnv *env = nullptr;
  GRBVar **vars = nullptr;
  GRBVar *uvars = nullptr;

  size_t total_vertices = num_depots + num_vertices;
  double_t success = 0.0;
  vars = new GRBVar *[total_vertices];
  for (size_t i = 0; i < total_vertices; ++i) {
    vars[i] = new GRBVar[total_vertices];
  }

  uvars = new GRBVar[total_vertices];

  try {
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    for (size_t i = 0; i < total_vertices; ++i) {
      for (size_t j = 0; j < total_vertices; ++j) {
        if (i == j) {
          vars[i][j] = model.addVar(
              0.0, 1.0, 0.0, GRB_BINARY,
              "x_" + std::to_string(i) + "_" + std::to_string(j));
        } else {
          vars[i][j] = model.addVar(
              0.0, 1.0, cost_mat[i][j], GRB_BINARY,
              "x_" + std::to_string(i) + "_" + std::to_string(j));
        }
      }
      uvars[i] = model.addVar(
          0, total_vertices + 1, 0, GRB_INTEGER, "u_" + std::to_string(i));
    }
    model.update();

    GRBLinExpr expr = 0;

    for (size_t i = 1; i < total_vertices; ++i) {
      expr += vars[0][i];
    }
    model.addConstr(expr == 0, "finish_out");

    expr = 0;
    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr += vars[i][0];
    }
    model.addConstr(expr == num_vehicles, "finish_in");

    for (size_t i = 1; i < num_depots; ++i) {
      expr = 0;
      for (size_t j = num_depots; j < total_vertices; ++j) {
        expr += vars[i][j];
      }
      model.addConstr(expr == 1, "out_depot_" + std::to_string(i));
    }

    for (size_t i = 1; i < num_depots; ++i) {
      expr = 0;
      for (size_t j = num_depots; j < total_vertices; ++j) {
        expr += vars[j][i];
      }
      model.addConstr(expr == 0, "in_depot_" + std::to_string(i));
    }

    // Degree constraints
    for (int j = num_depots; j < total_vertices; ++j) {
      expr = 0;
      for (int i = 1; i < total_vertices; ++i) {
        expr += vars[i][j];
      }
      model.addConstr(expr == 1, "deg_in_" + std::to_string(j));
    }

    for (int i = num_depots; i < total_vertices; ++i) {
      expr = 0;
      for (int j = num_depots; j < total_vertices; ++j) {
        expr += vars[i][j];
      }
      expr += vars[i][0];
      model.addConstr(expr == 1, "deg_out_" + std::to_string(i));
    }

    // Forbid edge from node back to itself
    for (size_t i = 0; i < total_vertices; ++i) {
      vars[i][i].set(GRB_DoubleAttr_UB, 0);
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      for (size_t j = num_depots; j < total_vertices; ++j) {
        expr = 0;
        if (i != j) {
          expr += uvars[i] - uvars[j] + L * vars[i][j] + (L - 2) * vars[j][i];
          model.addConstr(
              expr <= (L - 1),
              "subtour_elim_1_" + std::to_string(i) + std::to_string(j));
        }
      }
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;
      expr += uvars[i];
      GRBLinExpr depot_sum = 0;
      for (size_t j = 1; j < num_depots; ++j) {
        depot_sum += vars[j][i];
      }
      expr += (L - 2) * depot_sum;
      expr -= vars[i][0];
      model.addConstr(expr <= (L - 1), "subtour_elim_2_" + std::to_string(i));
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;
      expr += uvars[i];
      GRBLinExpr depot_sum = 0;
      for (size_t j = 1; j < num_depots; ++j) {
        depot_sum += vars[j][i];
      }
      expr += depot_sum;
      expr += (2 - K) * vars[i][0];
      model.addConstr(expr >= 2, "subtour_elim_3_" + std::to_string(i));
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;
      GRBLinExpr depot_sum = 0;
      for (size_t j = 1; j < num_depots; ++j) {
        depot_sum += vars[j][i];
      }
      expr += depot_sum + vars[i][0];
      model.addConstr(expr <= 1, "subtour_elim_4_" + std::to_string(i));
    }
/*
    //Subroute elimination constraints
    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;

      expr += uvars[i];

      for (size_t k = 0; k < num_depots; ++k) {
        expr += (L - 2) * vars[k][i];
      }

//      for (size_t k = 0; k < num_depots; ++k) {
//        expr += -vars[i][k];
//      }
      expr += -vars[i][0];

      model.addConstr(expr <= (L - 1), "subtour_elim_1_" + std::to_string(i));
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;

      expr += uvars[i];

      for (size_t k = 0; k < num_depots; ++k) {
        expr += vars[k][i];
      }

//      for (size_t k = 0; k < num_depots; ++k) {
//        expr += (2 - K) * vars[i][k];
//      }
      expr += (2 - K) * vars[i][0];
      model.addConstr(expr >= 2, "subtour_elim_2_" + std::to_string(i));
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;
      for (size_t k = 0; k < num_depots; ++k) {
        expr = vars[k][i] + vars[i][0];
        model.addConstr(expr <= 1,
                        "subtour_elim_3_" + std::to_string(i) + "_"
                            + std::to_string(k));
      }
    }

    for (size_t i = num_depots; i < total_vertices; ++i) {
      expr = 0;
      for (int j = num_depots; j < total_vertices; ++j) {
        if (i != j) {
          expr = uvars[i] - uvars[j] + L * vars[i][j] + (L - 2) * vars[j][i];
          model.addConstr(expr <= (L - 1),
                          "subtour_elim_4_" + std::to_string(i) + "_"
                              + std::to_string(j));
        }
      }
    }*/

    model.update();
//    model.write("/tmp/mdmtsp.lp");
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    model.optimize();
    Vector<std::string> solution;
    if (model.get(GRB_IntAttr_SolCount) > 0) {
      double value = model.get(GRB_DoubleAttr_ObjVal);
      for (size_t i = 0; i < total_vertices; ++i) {
        for (size_t j = 0; j < total_vertices; ++j) {
          if (vars[i][j].get(GRB_DoubleAttr_X) > 0.5) {
            solution.push_back(vars[i][j].get(GRB_StringAttr_VarName));
          }
        }
      }
      success = model.get(GRB_DoubleAttr_ObjVal);
    } else {
      success = 0.0;
    }
//    print_vector(solution);
    if (success > 0.0) {
      std::vector<std::string> split_string;
      //numDepots - 1 because no one is starting from the origin.
      for (size_t i = 0; i < num_depots - 1; ++i) {
        bool done = false;
        std::vector<uint_fast32_t> a_tour;
        split_string = split(solution[i], '_');
        int current_verex = atoi(split_string[1].c_str());
        int next_vertex = atoi(split_string[2].c_str());
        a_tour.push_back(current_verex);
        while (!done) {
          //Find the nextWP.
          //Set the nextWP as current and add it to the queue.
          //Set the the new nextWP.
          //If nextWP is 0 then exit the loop adding the nextWP in the vector.
          std::vector<std::string>::iterator tours_it;
          for (const std::string &path:solution) {
            split_string = split(path, '_');
            if (atoi(split_string[1].c_str()) == next_vertex) {
              current_verex = next_vertex;

              a_tour.push_back(current_verex);
              next_vertex = atoi(split_string[2].c_str());
              break;
            }
          }
          if (next_vertex == 0) {
            done = true;
          }
        }
        a_tour.push_back(0);
        tours.push_back(a_tour);
      }
    }
  } catch (GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  }

  catch (...) {
    std::cout << "Error during optimization" <<
              std::endl;
    success = 0.0;
  }
  for (int i = 0; i < total_vertices; i++)
    delete[] vars[i];
  delete[] vars;
  delete env;
  return success;
}