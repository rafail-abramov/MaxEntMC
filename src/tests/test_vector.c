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

    /* maxentmc_power_print(A->powers,stdout);

    maxentmc_power_t * pb = maxentmc_power_alloc_product(A->powers,A->powers);

    maxentmc_power_print_product(pb,stdout);

    maxentmc_matrix_t * B = maxentmc_matrix_alloc(pb,A->size1,0);

    gsl_matrix_set_zero((gsl_matrix *)B);

    maxentmc_matrix_power_multiply_add(A,A,B);

    maxentmc_matrix_print(A,stdout);

    maxentmc_matrix_print(B,stdout);
    */

    return 0;

}

