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
