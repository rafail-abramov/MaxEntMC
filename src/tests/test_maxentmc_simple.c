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

#include "test_maxentmc_simple.h"
#include "test_read_old_data.h"

#define MAXENTMC_DIM 2
#define QUAD_AMP 5
#define QUAD_SIZE (80*QUAD_AMP)
#define MAXENTMC_TOL 1E-06

#if (MAXENTMC_DIM == 2)
#define MAXENTMC_TEST_VAR "Mono_2D"
#endif

#if (MAXENTMC_DIM == 3)
#define MAXENTMC_TEST_VAR "Mono_3D"
#endif

#if (MAXENTMC_DIM == 4)
#define MAXENTMC_TEST_VAR "Mono_4D"
#endif

int test_maxentmc_simple(void)
{

    gsl_set_error_handler_off();

#if (MAXENTMC_DIM == 1)

    /** 1D test **/

    maxentmc_index_t pow;

    maxentmc_float_t x;

    maxentmc_list_t list = maxentmc_list_alloc(MAXENTMC_DIM,1,MAXENTMC_LIST_ORDERED,MAXENTMC_LIST_ASCEND);

    pow=0; x = 1;
    maxentmc_list_insert(list,pow,x);

    pow=1; x = 0;
    maxentmc_list_insert(list,pow,x);

    pow=2; x = 1;
    maxentmc_list_insert(list,pow,x);

    pow=3; x = 0.137193470359717;
    maxentmc_list_insert(list,pow,x);

    pow=4; x = 2.51091277331384;
    maxentmc_list_insert(list,pow,x);

    pow=5; x = 1.37475554727826;
    maxentmc_list_insert(list,pow,x);

    pow=6; x = 9.5608266563126;
    maxentmc_list_insert(list,pow,x);

    pow=7; x = 10.4940641636803;
    maxentmc_list_insert(list,pow,x);

    pow=8; x = 48.3127755980057;
    maxentmc_list_insert(list,pow,x);

    pow=9; x = 79.6716333150277;
    maxentmc_list_insert(list,pow,x);

    pow=10; x = 300.747351222544;
    maxentmc_list_insert(list,pow,x);

    pow=11; x = 636.826473835387;
    maxentmc_list_insert(list,pow,x);

    pow=12; x = 2190.07613030075;
    maxentmc_list_insert(list,pow,x);

    pow=13; x = 5433.37958129264;
    maxentmc_list_insert(list,pow,x);

    pow=14; x = 17984.9455750663;
    maxentmc_list_insert(list,pow,x);

    pow=15; x = 49502.7898338205;
    maxentmc_list_insert(list,pow,x);

    pow=16; x = 162276.494119514;
    maxentmc_list_insert(list,pow,x);

    maxentmc_power_vector_t constraints;
    maxentmc_list_create_power_vectors(list,&constraints);

#elif (MAXENTMC_DIM <= 4)

    char filename[256];

    sprintf(filename,"data/data_%uD/constraints_dim%u_pow%u_1.dat",MAXENTMC_DIM,MAXENTMC_DIM,2*(6-MAXENTMC_DIM));

    FILE * in = fopen(filename,"r");
    if(in == NULL){
        fputs("Could not open constraints file\n",stderr);
        return -1;
    }

    maxentmc_power_vector_t constraints = maxentmc_power_vector_fread_power(in);
    if(constraints == NULL){
        fputs("Could not read powers\n",stderr);
        return -1;
    }

    if(maxentmc_power_vector_fread_values(constraints,in)){
        fputs("Could not read values\n",stderr);
        return -1;
    }
    fclose(in);

#endif

#if (MAXENTMC_DIM > 4)
#error Dimension not supported
#endif

    /** Print the constraints **/

    puts("Constraints");
    maxentmc_power_vector_print(constraints,stdout);

    maxentmc_power_vector_t multipliers = maxentmc_power_vector_alloc(constraints);
    gsl_vector_memcpy(&multipliers->gsl_vec,&constraints->gsl_vec);

    size_t quad_size[MAXENTMC_DIM], i;
    maxentmc_float_t quad_start[MAXENTMC_DIM], quad_end[MAXENTMC_DIM];

    for(i=0;i<MAXENTMC_DIM;++i){
        quad_size[i] = QUAD_SIZE;
        quad_start[i] = -QUAD_AMP;
        quad_end[i] = QUAD_AMP;
    }

    if(maxentmc_basic_algorithm(multipliers, quad_size, quad_start, quad_end, MAXENTMC_TOL))
        puts("Some error happened");
    else{

        puts("Resulting multipliers:");
        maxentmc_power_vector_print(multipliers,stdout);

    }

    return 0;

}
