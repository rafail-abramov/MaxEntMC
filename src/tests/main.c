#include <stdio.h>

#include "test_list.h"
#include "test_vector.h"
#include "test_quad_gauss_1D.h"
#include "test_quad.h"
#include "test_gradient_hessian.h"
#include "test_maxentmc_simple.h"
#include "test_convert_constraints.h"

int main(int argc, char ** argv)
{
/*
    maxentmc_list_t list = test_list();

    test_vector(list);
*/
/*
    test_quad_gauss_1D();
*/
/*
    test_quadrature();
*/
/*
    test_gradient_hessian();
*/

    test_maxentmc_simple();

/*
    if(argc<3){
        fputs("Not enough command line arguments\n",stderr);
        return -1;
    }

    if(argc<4)
        test_convert_constraints(argv[1],argv[2],NULL);
    else
        test_convert_constraints(argv[1],argv[2],argv[3]);
*/
    return 0;

}
