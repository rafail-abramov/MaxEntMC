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

#ifndef MAXENTMC_GRADIENT_HESSIAN_H_INCLUDED
#define MAXENTMC_GRADIENT_HESSIAN_H_INCLUDED

#include "maxentmc_vector.h"

struct maxentmc_LGH_power_list_struct;

struct maxentmc_LGH_struct {

    /** This structure contains all necessary data (including constraints)
        to compute the Lagrangian, its gradient and Hessian from the
        quadrature data **/

    /** We need to store the constraints here **/

    struct maxentmc_power_struct * powers;
    maxentmc_index_t have_L;
    size_t L_element;

    /** In addition to the constraint_power, we are going to have a list
        of powers from which the extraction of the Lagrangian, gradient
        and Hessian is possible. This is done to speed up computations,
        so that it would not be necessary to recompute quadratures. **/

    struct maxentmc_LGH_power_list_struct * LGH_powers;

};

#endif // MAXENTMC_GRADIENT_HESSIAN_H_INCLUDED
