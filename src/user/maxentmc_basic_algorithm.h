/** This file is part of MaxEntMC, a maximum entropy algorithm with moment constraints. **/
/** Copyright (C) 2014 Rafail V. Abramov.                                               **/
/**                                                                                     **/
/** This program is free software: you can redistribute it and/or modify it under the   **/
/** terms of the GNU General Public License as published by the Free Software           **/
/** Foundation, either version 3 of the License, or (at your option) any later version. **/
/**                                                                                     **/
/** This program is distributed in the hope that it will be useful, but WITHOUT ANY     **/
/** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A     **/
/** PARTICULAR PURPOSE.  See the GNU General Public License for more details.           **/
/**                                                                                     **/
/** You should have received a copy of the GNU General Public License along with this   **/
/** program.  If not, see <http://www.gnu.org/licenses/>.                               **/

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
