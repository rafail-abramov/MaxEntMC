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
#include <gsl/gsl_vector.h>
#include "test_quad.h"

/** NOTE! For the moments of the gaussian density with zero mean and unit variance,
    the rectangle rule yields exact results with QUAD_SIZE=10 and QUAD_AMP=8.0 !!! **/

/*
#define QUAD_SIZE 10
#define QUAD_AMP 8.0
#define MAX_MOMENT 16
*/

#define MULT_AMP 2.0
#define MULT_POINTS 10

#define QUAD_SIZE 30
#define QUAD_AMP 4.0
#define MAXENTMC_DIM 2
#define MAXENTMC_POW 6

#define POWER_COUNTER_INIT(dim,counter,flag)    \
  maxentmc_index_t (counter)[(dim)];			    \
  {								                \
    maxentmc_index_t _i;						    \
    for(_i=0;_i<(dim);++_i)					    \
      (counter)[_i]=0;						    \
  }								                \
  size_t (flag) = 0

#define POWER_COUNTER_RESET(dim,counter,flag)   \
  {								                \
    maxentmc_index_t _i;						    \
    for(_i=0;_i<(dim);++_i)					    \
      (counter)[_i]=0;						    \
  }								                \
  (flag) = 0

#define POWER_COUNTER_INCREMENT(dim,pow,counter,flag)   \
  if((counter)[(dim)-1]<(pow)){						        \
    ++((counter)[0]);							            \
    ++(flag);								                \
    maxentmc_index_t _i=0;							        \
    while((flag)>(pow)){						            \
      (flag) -= (counter)[_i]-1;					        \
      (counter)[_i]=0;							            \
      ++((counter)[++_i]);						            \
    }									                    \
  }									                        \
  else (flag) = 0

/*
#define QUAD_COUNTER_INIT(dim,counter,flag) \
  size_t (counter)[(dim)];			        \
  {								            \
    maxentmc_index_t _i;				        \
    for(_i=0;_i<(dim);++_i)			        \
      (counter)[_i]=0;					    \
  }								            \
  maxentmc_index_t (flag) = 1

#define QUAD_COUNTER_RESET(dim,counter,flag)    \
  {								                \
    maxentmc_index_t _i;			                \
    for(_i=0;_i<(dim);++_i)			            \
      (counter)[_i]=0;				            \
  }                                             \
  (flag) = 1
*/
/** MUST CHECK THE INCREMENT **/
/*
#define QUAD_COUNTER_INCREMENT(dim,size,counter,flag)   \
    ++((counter)[0]);							        \
    if((counter)[(dim)-1]<(size)){				        \
        maxentmc_index_t _i=1;							\
        while((_i<(dim)) && ((counter)[_i-1]>=(size))){ \
            (counter)[_i-1] = 0;                        \
            ++((counter)[_i++]);                        \
        }                                               \
    }                                                   \
    else (flag) = 0
*/

#define CORNER_MOMENT 1
#define EACH_POWER_EVEN 2

#define TEST_QUAD_MAX(a,b)  ((a)>(b))?(a):(b)

int test_quadrature(void)
{

    maxentmc_float_t x=0.0;

    maxentmc_list_t list = maxentmc_list_alloc(MAXENTMC_DIM,1,MAXENTMC_LIST_ORDERED,MAXENTMC_LIST_ASCEND);

    POWER_COUNTER_INIT(MAXENTMC_DIM,pow,flag);

    maxentmc_list_insert_ca(list,pow,&x);

    maxentmc_power_vector_t moments;

    maxentmc_list_create_power_vectors(list,&moments);

//    maxentmc_vector_print(moments,stdout);

    maxentmc_list_clear(list);

    do{

        maxentmc_list_insert_ca(list,pow,&x);

        POWER_COUNTER_INCREMENT(MAXENTMC_DIM,MAXENTMC_POW,pow,flag);


    }while(flag);

    maxentmc_power_vector_t multipliers;

    maxentmc_list_create_power_vectors(list,&multipliers);

//    maxentmc_vector_print(multipliers,stdout);

    size_t i;

    puts("% Line 1: values of multipliers");
    for(i=0;i<(2*MULT_POINTS);++i){
        x = (0.5+i)*MULT_AMP/MULT_POINTS-MULT_AMP;
        printf(" %g",x);
    }
    puts("");

    maxentmc_quad_helper_t quad = maxentmc_quad_helper_alloc(MAXENTMC_DIM);

    size_t quad_num_points[MAXENTMC_DIM];
    maxentmc_float_t quad_start[MAXENTMC_DIM], quad_end[MAXENTMC_DIM];

    for(i=0;i<MAXENTMC_DIM;++i){
        quad_num_points[i] = 2*QUAD_SIZE;
        quad_start[i] = -QUAD_AMP;
        quad_end[i] = QUAD_AMP;
    }

    i = 0;
    flag = 0;

    do{

        size_t j;

        for(j=0;j<MAXENTMC_DIM;++j)
            pow[j] = 0;
        x=-1.0;
        for(j=0;j<MAXENTMC_DIM;++j){
            pow[j] = MAXENTMC_POW;
            size_t pos=0;
            if(maxentmc_power_vector_find_element_ca(multipliers,pow,&pos)){
                MAXENTMC_MESSAGE(stderr,"error: could not find the power");
                return -1;
            }
            gsl_vector_set(&multipliers->gsl_vec,pos,x);
            pow[j] = 0;
        }

        size_t special_case = 0, total_power = 0, max_power = 0, all_even = 1;
        maxentmc_power_vector_get_powers_ca(multipliers,i,pow);
        for(j=0;j<MAXENTMC_DIM;++j){
            total_power += pow[j];
            if(pow[j]&1)
                all_even = 0;
            max_power = TEST_QUAD_MAX(max_power,pow[j]);
        }
        if(all_even)
            special_case = EACH_POWER_EVEN;
        if(max_power == MAXENTMC_POW)
            special_case = CORNER_MOMENT;

        printf("%% Line %zu: moments [%u",i+2,pow[0]);
        for(j=1;j<MAXENTMC_DIM;++j)
            printf(" %u",pow[j]);
        puts("]");

        for(j=0;j<(2*MULT_POINTS);++j){

            x = (0.5+j)*MULT_AMP/MULT_POINTS-MULT_AMP;

            switch(special_case){
            case CORNER_MOMENT:
                x = -exp(x);
                break;
            case EACH_POWER_EVEN:
                x += 0.4*exp(-x);
                break;
            }

            gsl_vector_set(&multipliers->gsl_vec,i,x);

            //maxentmc_vector_print(multipliers,stdout);

            maxentmc_quad_helper_set_multipliers(quad,multipliers);

            maxentmc_quad_helper_set_moments(quad,moments);

            maxentmc_quadrature_rectangle_uniform_ca(quad, quad_num_points, quad_start, quad_end);

            maxentmc_quad_helper_get_moments(quad,moments);

            if(total_power==0)
                printf(" %g",moments->gsl_vec.data[0]-x);
            else
                printf(" %g",moments->gsl_vec.data[0]);

        }

        puts("");

        x=0.0;
        gsl_vector_set(&multipliers->gsl_vec,i,x);
/*
        for(j=0;j<MAXENTMC_DIM;++j)
            pow[j] = 0;
        x=-1.0;
        for(j=0;j<MAXENTMC_DIM;++j){
            pow[j] = MAXENTMC_POW;
            size_t pos=0;
            if(maxentmc_power_find(multipliers->powers,pow,&pos)){
                MAXENTMC_MESSAGE(stderr,"error: could not find the power");
                return -1;
            }
            maxentmc_vector_set(multipliers,pos,x);
            pow[j] = 0;
        }
*/
        ++i;
/*
        for(j=0,flag=0;j<MAXENTMC_DIM;++j)
            flag += multipliers->powers->power[i][j];
*/

    }while(i<multipliers->gsl_vec.size /*flag<MAXENTMC_POW*/);


    return 0;

}

