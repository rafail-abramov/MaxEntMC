#include "maxentmc_quad_rectangle_uniform.h"

int maxentmc_quadrature_rectangle_uniform(maxentmc_quad_helper_t const quad, ...)
{
    if(quad == NULL){
        fputs("maxentmc_quad_hausdorff_uniform: NULL pointer is given as quadrature helper structure",stderr);
        return -1;
    }

    maxentmc_index_t const dim = maxentmc_quad_helper_get_dimension(quad);

    size_t num_points[dim];
    maxentmc_float_t start[dim], end[dim];

    va_list ap;
    va_start(ap,quad);
    maxentmc_index_t i;
    for(i=0;i<dim;++i){
        num_points[i] = va_arg(ap,size_t);
        start[i] = va_arg(ap,maxentmc_float_t);
        end[i] = va_arg(ap,maxentmc_float_t);
    }
    va_end(ap);
    return maxentmc_quadrature_rectangle_uniform_ca(quad, num_points, start, end);

}

int maxentmc_quadrature_rectangle_uniform_ca(maxentmc_quad_helper_t const quad, size_t const * const num_points,
                                                 maxentmc_float_t const * const start, maxentmc_float_t const * const end)
{

    if(quad == NULL){
        fputs("maxentmc_quad_hausdorff_uniform: NULL pointer is given as quadrature helper structure",stderr);
        return -1;
    }

    maxentmc_index_t const dim = maxentmc_quad_helper_get_dimension(quad);

    size_t quad_point[dim];

    maxentmc_float_t abscissa[4][dim], dx[dim], weight = 1.0;

    maxentmc_index_t i;

    for(i=0;i<dim;++i){
        quad_point[i] = 0;
        dx[i] = (end[i] - start[i])/num_points[i];
        weight *= dx[i];
    }

    struct maxentmc_quad_helper_thread_struct * quad_thread = maxentmc_quad_helper_thread_alloc(quad);

    while(quad_point[dim-1]<num_points[dim-1]){

        for(i=1;i<dim;++i){
            abscissa[0][i] = start[i]+(0.5+quad_point[i])*dx[i];
            abscissa[1][i] = abscissa[0][i];
            abscissa[2][i] = abscissa[0][i];
            abscissa[3][i] = abscissa[0][i];
        }

        while((quad_point[0]+3) < num_points[0]){

            abscissa[0][0] = start[0]+(0.5+quad_point[0])*dx[0];
            abscissa[1][0] = start[0]+(1.5+quad_point[0])*dx[0];
            abscissa[2][0] = start[0]+(2.5+quad_point[0])*dx[0];
            abscissa[3][0] = start[0]+(3.5+quad_point[0])*dx[0];
            maxentmc_quad_helper_thread_compute_4(quad_thread,abscissa[0],weight,abscissa[1],weight,abscissa[2],weight,abscissa[3],weight);
            (quad_point[0]) += 4;

        }

        i = num_points[0] - quad_point[0];

        switch(i){
            case 3:
                abscissa[0][0] = start[0]+(0.5+quad_point[0])*dx[0];
                abscissa[1][0] = start[0]+(1.5+quad_point[0])*dx[0];
                abscissa[2][0] = start[0]+(2.5+quad_point[0])*dx[0];
                maxentmc_quad_helper_thread_compute_3(quad_thread,abscissa[0],weight,abscissa[1],weight,abscissa[2],weight);
                (quad_point[0]) += 3;
                break;
            case 2:
                abscissa[0][0] = start[0]+(0.5+quad_point[0])*dx[0];
                abscissa[1][0] = start[0]+(1.5+quad_point[0])*dx[0];
                maxentmc_quad_helper_thread_compute_2(quad_thread,abscissa[0],weight,abscissa[1],weight);
                (quad_point[0]) += 2;
                break;
            case 1:
                abscissa[0][0] = start[0]+(0.5+quad_point[0])*dx[0];
                maxentmc_quad_helper_thread_compute_1(quad_thread,abscissa[0],weight);
                (quad_point[0]) += 1;
                break;
        }

        i=1;

        while((i<dim) && (quad_point[i-1]==num_points[i-1])){

            quad_point[i-1] = 0;
            ++(quad_point[i++]);

        }

    }

    maxentmc_quad_helper_thread_merge(quad_thread);

    return 0;

}
