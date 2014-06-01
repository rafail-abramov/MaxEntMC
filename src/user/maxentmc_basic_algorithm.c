#include "maxentmc_basic_algorithm.h"

int maxentmc_basic_algorithm(maxentmc_power_vector_t const constraints, size_t const * const quad_size, maxentmc_float_t const * const quad_start,
                             maxentmc_float_t const * const quad_end, maxentmc_float_t const tolerance)
{

    /** First, check that the constraints are valid **/
    if(constraints == NULL){
        fputs(" MaxEntMC basic algorithm error: provided constraint vector is NULL\n",stderr);
        return -1;
    }



    /** Determine the dimension of the problem **/
    maxentmc_index_t const dimension = maxentmc_power_vector_get_dimension(constraints);



    /** Now we need to create a starting multiplier vector. For simplicity, we will start with Gaussian density, zero mean, identity covariance **/

    /** Allocate the array of powers**/
    maxentmc_index_t powers[dimension], i;

    /** Allocate the multiplier vector and zero it (using compatibility with GSL) **/
    maxentmc_power_vector_t multipliers = maxentmc_power_vector_alloc(constraints);
    gsl_vector_set_zero(&multipliers->gsl_vec);

    /** Find the position of element with zero total power and set it to what it should be in the Gaussian distribution **/
    for(i=0;i<dimension;++i)
        powers[i] = 0;
    size_t pos;
    if(maxentmc_power_vector_find_element_ca(multipliers,powers,&pos))
        return -1;
    gsl_vector_set(&multipliers->gsl_vec,pos,-log(sqrt(8.0*atan(1.0)))*dimension); /** Use the formula pi = 4 atan 1 **/

    /** Set the corner multipliers of power 2 to -1/2 **/
    for(i=0;i<dimension;++i){
        powers[i] = 2;
        if(maxentmc_power_vector_find_element_ca(multipliers,powers,&pos))
            return -1;
        gsl_vector_set(&multipliers->gsl_vec,pos,-0.5);
        powers[i] = 0;
    }

    size_t const size = constraints->gsl_vec.size; /** This is the total number of constraints in the problem **/

    maxentmc_power_vector_t moments_grad = maxentmc_power_vector_alloc(constraints); /** This is used to hold moments for gradient computation **/

    maxentmc_power_vector_t moments_hess = maxentmc_power_vector_product_alloc(constraints,constraints); /** This is used to hold moments for hessian computation **/

    maxentmc_quad_helper_t quad = maxentmc_quad_helper_alloc(dimension); /** This is quadrature helper structure **/

    maxentmc_LGH_t LGH = maxentmc_LGH_alloc(moments_grad); /** This is the object for computing the lagrangian, gradient and hessian from moments.
                                                         Allocated from any vector with constraint powers (gradient moments have suitable powers,
                                                         constraints vector could have been used too) **/

    maxentmc_LGH_add_power_vector(LGH,moments_hess);  /** Add the vector for hessian moments, to be able to extract the hessian **/

    /** Allocate the gradient, hessian and auxiliary data structures **/

    maxentmc_power_vector_t temp_multipliers = maxentmc_power_vector_alloc(constraints);
    gsl_vector * gradient = gsl_vector_alloc(size);
    gsl_matrix * hessian = gsl_matrix_alloc(size,size);
    gsl_vector * temp_gradient = gsl_vector_alloc(size);
    gsl_vector * step = gsl_vector_alloc(size);

    /** This is diagnostics (can be commented out) **/
    maxentmc_float_t lagrangian;
    maxentmc_gsl_matrix_t * eigvec = gsl_matrix_alloc(size,size);
    maxentmc_gsl_vector_t * eigval = gsl_vector_alloc(size);
    gsl_eigen_symm_workspace * eigen_workspace = gsl_eigen_symm_alloc(size);
    /** End diagnostics **/


    int do_it = 1; /** This is flag variable used to determine whether iterations ended or not **/
    int error_flag = 0; /** This is what is returned by this function. If non-zero, indicates error **/
    size_t num_iter=0; /** This is iteration counter **/

    do{

        /** This is start of the iterations **/

        /** First, compute the gradient and Hessian from the current set of Lagrange multipliers **/

        maxentmc_quad_helper_set_multipliers(quad,multipliers); /** Setting Lagrange multipliers for quadrature **/
        maxentmc_quad_helper_set_moments(quad,moments_hess);    /** Setting the moments for quadrature (currently hessian moments, since we will need the hessian at this stage **/
        maxentmc_quadrature_rectangle_uniform_ca(quad, quad_size, quad_start, quad_end); /** Use rectangular uniform quadrature **/
        maxentmc_quad_helper_get_moments(quad,moments_hess); /** Extract computed moments **/
        maxentmc_LGH_compute_gradient(LGH,moments_hess,constraints,gradient); /** Compute the gradient vector from the moments **/
        maxentmc_LGH_compute_hessian(LGH,moments_hess,hessian); /** Compute the hessian matrix from the same moments **/
        maxentmc_float_t gnorm = gsl_blas_dnrm2(gradient);   /** Compute the square norm of the gradient (if small enough, the iterations will be stopped) **/

        /** Diagnostic info (can be commented out) **/
        maxentmc_LGH_compute_lagrangian(LGH,moments_hess,constraints,multipliers,&lagrangian);
        printf("----------- Iteration %zu -----------\nValue of Lagrangian %g\n",++num_iter,lagrangian);
        printf("Norm of gradient %g\n",gnorm);
        gsl_matrix_memcpy(eigvec,hessian);
        gsl_eigen_symm(eigvec,eigval,eigen_workspace);
        printf("Hessian condition number %g\n",eigval->data[0]/eigval->data[size-1]);
        /** End diagnostic info **/

        if(isnan(gnorm) || isinf(gnorm) || (gnorm<tolerance)){
            do_it = 0; /** The iterations are stopped **/
            if(isnan(gnorm) || isinf(gnorm))
                error_flag = -1; /** Something failed, will return error **/
        }
        else{

            /** Do stepping here **/

            /** First, determine the step through Cholesky decomposition **/

            if(gsl_linalg_cholesky_decomp(hessian)){

                do_it = 0;
                error_flag = -1;
                puts("Cholesky decomposition failed, convergence failed");

            }
            else{

                gsl_linalg_cholesky_solve(hessian,gradient,step);

                /** Here do simple line search (halving the distance if line minimum is overshot) **/

                int do_line_search = 1;
                size_t num_line_search=0;
                maxentmc_float_t step_scale = -1.0;

                do{

                    /** Here we compute the temporary set of multipliers from the step **/
                    gsl_vector_memcpy(&temp_multipliers->gsl_vec,&multipliers->gsl_vec);
                    gsl_blas_daxpy(step_scale,step,&temp_multipliers->gsl_vec);
                    /** Done computing temporary multipliers **/

                    maxentmc_quad_helper_set_multipliers(quad,temp_multipliers); /** Set the temporary multipliers in the quadrature **/
                    maxentmc_quad_helper_set_moments(quad,moments_grad);  /** Here we do not need hessian, so set gradient moments (faster computation) **/
                    maxentmc_quadrature_rectangle_uniform_ca(quad, quad_size, quad_start, quad_end); /** Compute quadrature **/
                    maxentmc_quad_helper_get_moments(quad,moments_grad); /** Extract moments **/
                    maxentmc_LGH_compute_gradient(LGH,moments_grad,constraints,temp_gradient); /** Compute the temporary gradient **/

                    maxentmc_float_t gdot;
                    gsl_blas_ddot(step,temp_gradient,&gdot);
                    if(isnan(gdot) || isinf(gdot) || (gdot<0)){

                        /** Here one needs to define a guard against too many line search iterations, but we are too lazy for it (this is an example algorithm anyway) **/

                        step_scale *= 0.5;
                        ++num_line_search;
                    }
                    else{
                        /** Line search successful, copy the temporary multipliers into the main multipliers **/
                        gsl_vector_memcpy(&multipliers->gsl_vec,&temp_multipliers->gsl_vec);
                        do_line_search = 0;
                    }

                }while(do_line_search);

                /** Diagnostic info (can be commented out) **/

                printf("Line search rescalings %zu\n",num_line_search);

                /** End diagnostic info **/

            }

        }

    }while(do_it);

    if(!error_flag){

        /** Computations are successful, copy the computed multipliers into the constraint vector **/
        gsl_vector_memcpy(&constraints->gsl_vec,&multipliers->gsl_vec);

    }

    /** Release everything allocated for the computation **/

    gsl_matrix_free(eigvec);
    gsl_vector_free(eigval);
    gsl_eigen_symm_free(eigen_workspace);

    maxentmc_power_vector_free(temp_multipliers);
    gsl_vector_free(gradient);
    gsl_matrix_free(hessian);
    gsl_vector_free(temp_gradient);
    gsl_vector_free(step);

    maxentmc_power_vector_free(moments_grad);
    maxentmc_power_vector_free(moments_hess);
    maxentmc_quad_helper_free(quad);
    maxentmc_LGH_free(LGH);

    maxentmc_power_vector_free(multipliers);

    return error_flag;

}
