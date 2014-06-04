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

#ifndef MAXENTMC_H_INCLUDED
#define MAXENTMC_H_INCLUDED

#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>

/********* Basic type definitions ********/

typedef uint8_t maxentmc_index_t;

#ifdef MAXENTMC_SINGLE_PRECISION
typedef float maxentmc_float_t;
#else
typedef double maxentmc_float_t;
#endif

#ifndef __GSL_BLOCK_H__
struct gsl_block_struct;
#endif

/***** Vector and matrix definitions *****/

#ifndef __GSL_VECTOR_H__
typedef struct {
    size_t size, stride;
    maxentmc_float_t * data;
    struct gsl_block_struct * block;
    int owner;
} maxentmc_gsl_vector_t;
#else
#ifdef MAXENTMC_SINGLE_PRECISION
typedef gsl_vector_float maxentmc_gsl_vector_t;
#else
typedef gsl_vector maxentmc_gsl_vector_t;
#endif
#endif

#ifndef __GSL_MATRIX_H__
typedef struct {
    size_t size1, size2, tda;
    maxentmc_float_t * data;
    struct gsl_block_struct * block;
    int owner;
} maxentmc_gsl_matrix_t;
#else
#ifdef MAXENTMC_SINGLE_PRECISION
typedef gsl_matrix_float maxentmc_gsl_matrix_t;
#else
typedef gsl_matrix maxentmc_gsl_matrix_t;
#endif
#endif

/**** Power type declaration and maxent power-vector type definition ****/

struct maxentmc_power_struct;
struct maxentmc_power_vector_struct {
    maxentmc_gsl_vector_t gsl_vec;
    struct maxentmc_power_struct * powers;
};
typedef struct maxentmc_power_vector_struct * maxentmc_power_vector_t;

/** Vector functions **/

struct maxentmc_power_vector_struct * maxentmc_power_vector_alloc(struct maxentmc_power_vector_struct const *);

struct maxentmc_power_vector_struct * maxentmc_power_vector_product_alloc(struct maxentmc_power_vector_struct const *,
                                                                          struct maxentmc_power_vector_struct const *);

void maxentmc_power_vector_free(struct maxentmc_power_vector_struct * const);

int maxentmc_power_vector_find_element(struct maxentmc_power_vector_struct const * v, ...);

int maxentmc_power_vector_find_element_ca(struct maxentmc_power_vector_struct const * v, maxentmc_index_t const * p, size_t * pos);

int maxentmc_power_vector_get_powers(struct maxentmc_power_vector_struct const * v, size_t pos, ...);

int maxentmc_power_vector_get_powers_ca(struct maxentmc_power_vector_struct const * v, size_t pos, maxentmc_index_t * p);

maxentmc_index_t maxentmc_power_vector_get_dimension(struct maxentmc_power_vector_struct const *);

int maxentmc_power_vector_get_max_power(struct maxentmc_power_vector_struct const *, maxentmc_index_t * max_power);

int maxentmc_power_vector_print(struct maxentmc_power_vector_struct const *, FILE * out);

struct maxentmc_power_vector_struct * maxentmc_power_vector_fread_power(FILE * in);

int maxentmc_power_vector_fwrite_power(struct maxentmc_power_vector_struct const *, FILE * out);

int maxentmc_power_vector_fread_values(struct maxentmc_power_vector_struct *, FILE * in);

int maxentmc_power_vector_fwrite_values(struct maxentmc_power_vector_struct const *, FILE * out);

int maxentmc_power_vector_compute_polynomial(struct maxentmc_power_vector_struct const * v, ...);

int maxentmc_power_vector_compute_polynomial_ca(struct maxentmc_power_vector_struct const * v,
                                                maxentmc_float_t const * x, maxentmc_float_t * result);

/** List structure declaration **/

enum MAXENTMC_LIST_POWER_ORDER {MAXENTMC_LIST_FORWARD, MAXENTMC_LIST_BACKWARD, MAXENTMC_LIST_ORDERED};
enum MAXENTMC_LIST_POWER_INSERTION_ORDER {MAXENTMC_LIST_ASCEND, MAXENTMC_LIST_DESCEND};

struct maxentmc_list_struct;
typedef struct maxentmc_list_struct * maxentmc_list_t;

/** List functions **/

struct maxentmc_list_struct * maxentmc_list_alloc(maxentmc_index_t dimension, size_t data_size,
                                                  enum MAXENTMC_LIST_POWER_ORDER order,
                                                  enum MAXENTMC_LIST_POWER_INSERTION_ORDER insertion_order);

int maxentmc_list_clear(struct maxentmc_list_struct *);

void maxentmc_list_free(struct maxentmc_list_struct *);

int maxentmc_list_insert(struct maxentmc_list_struct * const, ...);

int maxentmc_list_insert_ca(struct maxentmc_list_struct *, maxentmc_index_t const * powers, maxentmc_float_t const * data);

int maxentmc_list_delete(struct maxentmc_list_struct *, ...);

int maxentmc_list_delete_ca(struct maxentmc_list_struct *, maxentmc_index_t const * powers);

int maxentmc_list_print(struct maxentmc_list_struct const *, FILE * out);

/** Create vectors from a list **/

struct maxentmc_power_vector_struct ** maxentmc_list_create_power_vector_array(struct maxentmc_list_struct const *);

int maxentmc_list_create_power_vectors(struct maxentmc_list_struct const *, ...);

/** Quadrature helper structures and functions **/

struct maxentmc_quad_helper_struct;
typedef struct maxentmc_quad_helper_struct * maxentmc_quad_helper_t;
struct maxentmc_quad_helper_thread_struct;
typedef struct maxentmc_quad_helper_thread_struct * maxentmc_quad_helper_thread_t;

struct maxentmc_quad_helper_struct * maxentmc_quad_helper_alloc(maxentmc_index_t dimension);

void maxentmc_quad_helper_free(struct maxentmc_quad_helper_struct * q);

maxentmc_index_t maxentmc_quad_helper_get_dimension(struct maxentmc_quad_helper_struct const * q);

int maxentmc_quad_helper_set_shift_rotation(struct maxentmc_quad_helper_struct * q, struct maxentmc_power_vector_struct const * constraints);

int maxentmc_quad_helper_set_multipliers(struct maxentmc_quad_helper_struct * q, struct maxentmc_power_vector_struct const * multipliers);

int maxentmc_quad_helper_set_moments(struct maxentmc_quad_helper_struct * q, struct maxentmc_power_vector_struct const * moments);

int maxentmc_quad_helper_get_moments(struct maxentmc_quad_helper_struct * q, struct maxentmc_power_vector_struct * moments);

struct maxentmc_quad_helper_thread_struct * maxentmc_quad_helper_thread_alloc(struct maxentmc_quad_helper_struct *);

int maxentmc_quad_helper_thread_merge(struct maxentmc_quad_helper_thread_struct *);

int maxentmc_quad_helper_thread_compute(struct maxentmc_quad_helper_thread_struct *, ...);

int maxentmc_quad_helper_thread_compute_1(struct maxentmc_quad_helper_thread_struct *,
                                          maxentmc_float_t const * x1, maxentmc_float_t w1);
/** Length of x is [dimension] **/

int maxentmc_quad_helper_thread_compute_2(struct maxentmc_quad_helper_thread_struct *,
                                          maxentmc_float_t const * x1, maxentmc_float_t w1,
                                          maxentmc_float_t const * x2, maxentmc_float_t w2);

int maxentmc_quad_helper_thread_compute_3(struct maxentmc_quad_helper_thread_struct *,
                                          maxentmc_float_t const * x1, maxentmc_float_t w1,
                                          maxentmc_float_t const * x2, maxentmc_float_t w2,
                                          maxentmc_float_t const * x3, maxentmc_float_t w3);

int maxentmc_quad_helper_thread_compute_4(struct maxentmc_quad_helper_thread_struct *,
                                          maxentmc_float_t const * x1, maxentmc_float_t w1,
                                          maxentmc_float_t const * x2, maxentmc_float_t w2,
                                          maxentmc_float_t const * x3, maxentmc_float_t w3,
                                          maxentmc_float_t const * x4, maxentmc_float_t w4);


/** Lagrangian, gradient and Hessian structures and functions **/

struct maxentmc_LGH_struct;
typedef struct maxentmc_LGH_struct * maxentmc_LGH_t;

/** Allocation/deallocation routines **/

struct maxentmc_LGH_struct * maxentmc_LGH_alloc(struct maxentmc_power_vector_struct const * constraints);

void maxentmc_LGH_free(struct maxentmc_LGH_struct *);


/** List addition routines **/

int maxentmc_LGH_add_power_vector(struct maxentmc_LGH_struct * d,
                                  struct maxentmc_power_vector_struct const * p);

/** Computation routines **/

int maxentmc_LGH_compute_lagrangian(struct maxentmc_LGH_struct const * d,
                                    struct maxentmc_power_vector_struct const * moments,
                                    struct maxentmc_power_vector_struct const * constraints,
                                    struct maxentmc_power_vector_struct const * multipliers,
                                    maxentmc_float_t * L);

int maxentmc_LGH_compute_gradient(struct maxentmc_LGH_struct const * d,
                                  struct maxentmc_power_vector_struct const * moments,
                                  struct maxentmc_power_vector_struct const * constraints,
                                  maxentmc_gsl_vector_t * G);

int maxentmc_LGH_compute_hessian(struct maxentmc_LGH_struct const * d,
                                 struct maxentmc_power_vector_struct const * moments,
                                 maxentmc_gsl_matrix_t * H);


#endif // MAXENTMC_H_INCLUDED
