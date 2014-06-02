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

#ifndef TEST_QUAD_HAUSDORFF_UNIFORM_H_INCLUDED
#define TEST_QUAD_HAUSDORFF_UNIFORM_H_INCLUDED

#include <stdio.h>
#include <stdarg.h>
#include "../user/maxentmc.h"

int maxentmc_quadrature_rectangle_uniform(maxentmc_quad_helper_t const quad, ...);

int maxentmc_quadrature_rectangle_uniform_ca(maxentmc_quad_helper_t const quad, size_t const * const num_points,
                                                 maxentmc_float_t const * const start, maxentmc_float_t const * const end);

#endif // TEST_QUAD_HAUSDORFF_UNIFORM_H_INCLUDED
