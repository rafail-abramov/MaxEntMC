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

#ifndef MAXENTMC_VECTOR_MATRIX_H_INCLUDED
#define MAXENTMC_VECTOR_MATRIX_H_INCLUDED

#include "maxentmc_power.h"
#include "../user/maxentmc.h"

size_t maxentmc_power_vector_size(size_t const s);

struct maxentmc_power_vector_struct * maxentmc_power_vector_alloc_from_power(struct maxentmc_power_struct * const);

int maxentmc_power_vector_init(struct maxentmc_power_struct const * const, struct maxentmc_power_vector_struct * const);

/*
int maxentmc_vector_power_multiply_add(struct maxentmc_power_vector_struct const * const v1, struct maxentmc_power_vector_struct const * const v2, struct maxentmc_power_vector_struct * const v);
*/

#endif // MAXENTMC_VECTOR_MATRIX_H_INCLUDED
