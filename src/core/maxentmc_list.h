#ifndef MAXENTMC_LIST_H_INCLUDED
#define MAXENTMC_LIST_H_INCLUDED

#include "../user/maxentmc.h"
#include "maxentmc_vector.h"

struct maxentmc_list_struct {

    maxentmc_index_t dimension;
    enum MAXENTMC_LIST_POWER_ORDER order;
    enum MAXENTMC_LIST_POWER_INSERTION_ORDER insertion_order;
    size_t data_size, list_size;
    struct maxentmc_list_link_struct * link;

};

struct maxentmc_power_struct * maxentmc_list_create_power(struct maxentmc_list_struct const * const);

#endif // MAXENTMC_LIST_H_INCLUDED
