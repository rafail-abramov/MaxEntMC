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
