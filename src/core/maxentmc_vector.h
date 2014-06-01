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
