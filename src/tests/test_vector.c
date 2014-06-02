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

#include "test_vector.h"

int test_vector(maxentmc_list_t const list)
{

    maxentmc_power_vector_t v1, v2, v;

    maxentmc_list_create_power_vectors(list,&v1,&v2);

    maxentmc_power_vector_print(v1,stdout);

    maxentmc_power_vector_print(v2,stdout);

    v = maxentmc_power_vector_product_alloc(v1,v2);

    gsl_vector_set_zero(&v->gsl_vec);

    maxentmc_power_vector_print(v,stdout);

    maxentmc_power_vector_free(v);
    maxentmc_power_vector_free(v1);
    maxentmc_power_vector_free(v2);

    return 0;

}

