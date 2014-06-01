#ifndef MAXENTMC_BASIC_ALGORITHM_H_INCLUDED
#define MAXENTMC_BASIC_ALGORITHM_H_INCLUDED

#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include "../user/maxentmc.h"
#include "../user/maxentmc_quad_rectangle_uniform.h"

int maxentmc_basic_algorithm(maxentmc_power_vector_t const v, size_t const * const quad_size, maxentmc_float_t const * const quad_start,
                             maxentmc_float_t const * const quad_end, maxentmc_float_t const tolerance);
/** On input, v contains input contraints. On successful output, v contains computed Lagrange multipliers.
    quad_size, quad_start and quad_end are inputs to maxentmc_quad_rectangle_uniform_ca **/

#endif // MAXENTMC_BASIC_ALGORITHM_H_INCLUDED
