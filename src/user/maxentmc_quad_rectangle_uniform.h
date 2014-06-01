#ifndef TEST_QUAD_HAUSDORFF_UNIFORM_H_INCLUDED
#define TEST_QUAD_HAUSDORFF_UNIFORM_H_INCLUDED

#include <stdio.h>
#include <stdarg.h>
#include "../user/maxentmc.h"

int maxentmc_quadrature_rectangle_uniform(maxentmc_quad_helper_t const quad, ...);

int maxentmc_quadrature_rectangle_uniform_ca(maxentmc_quad_helper_t const quad, size_t const * const num_points,
                                                 maxentmc_float_t const * const start, maxentmc_float_t const * const end);

#endif // TEST_QUAD_HAUSDORFF_UNIFORM_H_INCLUDED
