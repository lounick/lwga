//
// Created by nick on 16/03/18.
//

#ifndef LWGA_MDMTSP_SOLVER_H
#define LWGA_MDMTSP_SOLVER_H

#include <ga_types.h>
#include <ga_utils.h>

double_t runMDMTSP(
    const uint_fast8_t num_vehicles, const uint_fast32_t num_vertices,
    const uint_fast8_t num_depots, const Matrix<double_t> &cost_mat,
    const bool fair, Matrix<uint_fast32_t> &tours);

#endif //LWGA_MDMTSP_SOLVER_H
