#ifndef MAXENTMC_POWER_H_INCLUDED
#define MAXENTMC_POWER_H_INCLUDED

#include <stdio.h>
#include "../user/maxentmc.h"
#include "maxentmc_defs.h"

#define MAXENTMC_POWER_COMPLETE 1
#define MAXENTMC_POWER_ORDERED 2

#define MAXENTMC_POWER_SET_PROPERTY(_prop,_pow) ((_pow)->properties)|=(_prop)
#define MAXENTMC_POWER_CLEAR_PROPERTY(_prop,_pow) ((_pow)->properties)&=~(_prop)
#define MAXENTMC_POWER_CHECK_PROPERTY(_prop,_pow) (((_pow)->properties)&(_prop))

#define MAXENTMC_POWER_STREAM_HEADER "MAXENTMCPOWER"
#define MAXENTMC_POWER_STREAM_HEADER_SIZE 13

size_t maxentmc_bincoeff(size_t const n, size_t const k);

int maxentmc_power_comparison(maxentmc_index_t const dim, maxentmc_index_t const * const p1, maxentmc_index_t const * const p2);

/** Returns -1 if p1 is "less" than p2, 1 if p1 is "greater" than p2, and 0 if "equal" **/

/** The data structure and routines for the product of polynomials **/

/***** Power structure definition ****/

struct maxentmc_product_struct;

struct maxentmc_power_struct {

    unsigned char properties; /** Stores up to eight flags **/
    maxentmc_index_t dimension, max_power, num_refs; /** num_refs stores the number of references from vectors. When zero, the power is freed **/
    size_t size;
    maxentmc_index_t * max_power_per_dimension; /** [dimension] **/
    maxentmc_index_t ** power; /** [size][dimension] **/
    struct maxentmc_product_struct * product;

};

int maxentmc_power_find(struct maxentmc_power_struct const * const d, maxentmc_index_t const * const p, size_t * pos);

int maxentmc_power_get_power(struct maxentmc_power_struct const * const d, size_t const pos, maxentmc_index_t * const p);

struct maxentmc_power_struct * maxentmc_power_alloc(maxentmc_index_t const dimension, size_t const size);

struct maxentmc_power_struct * maxentmc_power_alloc_product(struct maxentmc_power_struct const * const p1, struct maxentmc_power_struct const * const p2);

void maxentmc_power_free(struct maxentmc_power_struct * const); /** Decrements num_refs and uses free() to deallocate **/

int maxentmc_power_update_max_power(struct maxentmc_power_struct * const p);

int maxentmc_power_print(struct maxentmc_power_struct const * const, FILE * const out);

struct maxentmc_power_struct * maxentmc_power_fread(FILE * const in);

int maxentmc_power_fwrite(struct maxentmc_power_struct const * const, FILE * const out);

/** Structures for product **/

struct maxentmc_product_sublink_struct {
    size_t i1, i2;
};

struct maxentmc_product_link_struct {
    size_t np;
    struct maxentmc_product_sublink_struct * p;
};

struct maxentmc_product_struct{

    struct maxentmc_power_struct * p1, * p2;
    struct maxentmc_product_link_struct * entry;
    /** product[size], for each entry, np contains the length of p12, where p[i].i1 and p[i].i2 contain the indices of p1 and p2 to multiply **/

};

int maxentmc_power_print_product(struct maxentmc_power_struct const * const d, FILE * const out);


int maxentmc_power_multiply_add(size_t const N,
                                struct maxentmc_power_struct const * const p1, maxentmc_index_t const x1_power_dim, maxentmc_float_t const * const * const x1,
                                struct maxentmc_power_struct const * const p2, maxentmc_index_t const x2_power_dim, maxentmc_float_t const * const * const x2,
                                struct maxentmc_power_struct const * const p, maxentmc_index_t const x_power_dim, maxentmc_float_t * const * const x);

#endif // MAXENTMC_POWER_H_INCLUDED
