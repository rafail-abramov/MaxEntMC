#include <math.h>
#include "test_quad_gauss_1D.h"

/** NOTE! For the moments of the gaussian density with zero mean and unit variance,
    the rectangle rule yields exact results with QUAD_SIZE=10 and QUAD_AMP=8.0 !!! **/

#define QUAD_SIZE 10
#define QUAD_AMP 8.0
#define MAXENTMC_POW 16

int test_quad_gauss_1D(void)
{

    maxentmc_index_t p;
    maxentmc_float_t x;

    maxentmc_list_t list = maxentmc_list_alloc(1,1,MAXENTMC_LIST_ORDERED,MAXENTMC_LIST_DESCEND);

    maxentmc_float_t const pi = 4.0*atan(1.0);

    p = 0;
    x = -log(sqrt(2.0*pi));
    maxentmc_list_insert(list,p,x);

    p = 2;
    x = -0.5;
    maxentmc_list_insert(list,p,x);

    maxentmc_list_print(list,stdout);

    maxentmc_power_vector_t multipliers;

    maxentmc_list_create_power_vectors(list,&multipliers);

    maxentmc_power_vector_print(multipliers,stdout);

    maxentmc_list_clear(list);

    x = 0;

    for(p=0;p<=MAXENTMC_POW;++p)
        maxentmc_list_insert(list,p,x);

    maxentmc_power_vector_t moments;

    maxentmc_list_create_power_vectors(list,&moments);

    maxentmc_list_free(list);

    maxentmc_quad_helper_t quad = maxentmc_quad_helper_alloc(1);

    maxentmc_quad_helper_set_multipliers(quad,multipliers);

    maxentmc_quad_helper_set_moments(quad,moments);

    maxentmc_quadrature_rectangle_uniform(quad, 2*QUAD_SIZE, -QUAD_AMP, QUAD_AMP);

    maxentmc_quad_helper_get_moments(quad,moments);

    maxentmc_quad_helper_free(quad);

    maxentmc_power_vector_print(moments,stdout);

    maxentmc_power_vector_free(moments);

    maxentmc_power_vector_free(multipliers);

    return 0;

}

