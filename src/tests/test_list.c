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

#include "test_list.h"

#define DIMENSION 2
#define DATA_SIZE 2
#define SIZE 100

maxentmc_list_t test_list(void)
{

    maxentmc_index_t power[DIMENSION];
    maxentmc_float_t data[DATA_SIZE], data1, data2;

    maxentmc_list_t list = maxentmc_list_alloc(DIMENSION,DATA_SIZE,MAXENTMC_LIST_BACKWARD/*MAXENTMC_LIST_ORDERED*/,MAXENTMC_LIST_ASCEND/*MAXENTMC_LIST_DESCEND*/);

    power[0] = 2; power[1] = 0;
    data1 = 0.5; data2 = -0.5;
    maxentmc_list_insert(list,power[0],power[1],data1,data2);

    power[0] = 0; power[1] = 1;
    data[0] = -0.3; data[1] = 1;
    maxentmc_list_insert_ca(list,power,data);

    power[0] = 1; power[1] = 0;
    data[0] = 1; data[1] = -2;
    maxentmc_list_insert_ca(list,power,data);

    power[0] = 0; power[1] = 0;
    data[0] = 0; data[1] = -1;
    maxentmc_list_insert_ca(list,power,data);

    power[0] = 0; power[1] = 2;
    data[0] = -0.3; data[1] = 1;
    maxentmc_list_insert_ca(list,power,data);

    power[0] = 1; power[1] = 1;
    data[0] = -0.5; data[1] = -1;
    maxentmc_list_insert_ca(list,power,data);

    maxentmc_list_print(list,stdout);

/*
    maxentmc_list_delete_ca(list,power);
    maxentmc_list_print(list,stdout);

    power[0] = 0; power[1] = 0;
    maxentmc_list_delete_ca(list,power);

    maxentmc_list_print(list,stdout);

    power[0] = 0; power[1] = 1;
    maxentmc_list_delete_ca(list,power);

    maxentmc_list_print(list,stdout);

    power[0] = 2; power[1] = 0;
    maxentmc_list_delete_ca(list,power);

    maxentmc_list_print(list,stdout);
*/

    return list;

}
