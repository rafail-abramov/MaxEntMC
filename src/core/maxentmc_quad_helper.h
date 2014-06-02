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

#ifndef MAXENTMC_QUAD_HELPER_H_INCLUDED
#define MAXENTMC_QUAD_HELPER_H_INCLUDED

#include "maxentmc_vector.h"

#define MAXENTMC_QUAD_THREAD_HOWMANY_AT_ONCE 8

struct maxentmc_quad_helper_power_list_struct {
    struct maxentmc_power_vector_struct * power_vector;
    struct maxentmc_quad_helper_power_list_struct * next;
};

struct maxentmc_quad_helper_struct {


    maxentmc_index_t dimension, max_power, shift_rotate, locked;
    maxentmc_index_t n_mult, n_mom;
    maxentmc_float_t scale, * shift, * rotate;

    struct maxentmc_power_vector_struct * multipliers, * moments;

    struct maxentmc_quad_helper_power_list_struct * multiplier_list, * moment_list;

};

struct maxentmc_quad_helper_thread_struct {

    struct maxentmc_quad_helper_struct * main_quadrature;
    maxentmc_float_t * moments;
    maxentmc_float_t * scratch;

};

#endif // MAXENTMC_QUAD_HELPER_H_INCLUDED
