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

#include <math.h>
#include "test_gradient_hessian.h"

#define MAXENTMC_DIM 2
#define QUAD_SIZE 100
#define QUAD_AMP 8.0

int test_gradient_hessian(void)
{
    maxentmc_float_t const pi = 4.0*atan(1.0);

    maxentmc_list_t list = maxentmc_list_alloc(2,2,MAXENTMC_LIST_ORDERED,MAXENTMC_LIST_ASCEND);

    maxentmc_index_t pow1, pow2;

    maxentmc_float_t x, y;

    pow1=0; pow2 = 0; x = -2.0*log(sqrt(2.0*pi)); y = 1;
    maxentmc_list_insert(list,pow1,pow2,x,y);

    pow1=1; pow2 = 0; x = 0; y = 0;
    maxentmc_list_insert(list,pow1,pow2,x,y);

    pow1=0; pow2 = 1; x = 0; y = 0;
    maxentmc_list_insert(list,pow1,pow2,x,y);

    pow1=2; pow2 = 0; x = -0.5; y = 1;
    maxentmc_list_insert(list,pow1,pow2,x,y);

    pow1=0; pow2 = 2; x = -0.5; y = 1;
    maxentmc_list_insert(list,pow1,pow2,x,y);

    pow1=1; pow2 = 1; x = 0; y = 0;
    maxentmc_list_insert(list,pow1,pow2,x,y);

    maxentmc_power_vector_t multipliers, constraints;

    maxentmc_list_create_power_vectors(list,&multipliers,&constraints);

    maxentmc_power_vector_t moments = maxentmc_power_vector_product_alloc(multipliers,multipliers);

    maxentmc_LGH_t LGH = maxentmc_LGH_alloc(constraints);

    maxentmc_LGH_add_power_vector(LGH,moments);

    maxentmc_quad_helper_t quad = maxentmc_quad_helper_alloc(MAXENTMC_DIM);

    maxentmc_quad_helper_set_multipliers(quad,multipliers);

    maxentmc_quad_helper_set_moments(quad,moments);

    size_t quad_num_points[MAXENTMC_DIM], i;
    maxentmc_float_t quad_start[MAXENTMC_DIM], quad_end[MAXENTMC_DIM];

    for(i=0;i<MAXENTMC_DIM;++i){
        quad_num_points[i] = 2*QUAD_SIZE;
        quad_start[i] = -QUAD_AMP;
        quad_end[i] = QUAD_AMP;
    }

    maxentmc_quadrature_rectangle_uniform_ca(quad, quad_num_points, quad_start, quad_end);

    maxentmc_quad_helper_get_moments(quad,moments);

    puts("Moments:");
    maxentmc_power_vector_print(moments,stdout);

    maxentmc_float_t lagrangian;

    maxentmc_LGH_compute_lagrangian(LGH,moments,constraints,multipliers,&lagrangian);

    printf("Lagrangian = %g\n",lagrangian);

    maxentmc_gsl_vector_t * gradient = gsl_vector_alloc(constraints->gsl_vec.size);

    maxentmc_LGH_compute_gradient(LGH,moments,constraints,gradient);

//    puts("Gradient:");
//    maxentmc_vector_print(gradient,stdout);

    maxentmc_gsl_matrix_t * hessian = gsl_matrix_alloc(constraints->gsl_vec.size,constraints->gsl_vec.size);

    maxentmc_LGH_compute_hessian(LGH,moments,hessian);

//    puts("Hessian:");
//    maxentmc_matrix_print(hessian,stdout);

    return 0;

}
