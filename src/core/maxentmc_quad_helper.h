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
