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

#include <string.h>
#include <pthread.h>
#include <math.h>

#include "maxentmc_quad_helper.h"

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

#define QUAD_MAX(a,b)  ((a)>(b))?(a):(b)

#define MAXENTMC_QUAD_HELPER_HEADER_SIZE MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_float_t),sizeof(struct maxentmc_quad_helper_struct))

#define MAXENTMC_QUAD_HELPER_SIZE(_s_) MAXENTMC_ALIGNED_SIZE(MAXENTMC_FLOAT_ALIGNMENT,MAXENTMC_QUAD_HELPER_HEADER_SIZE+sizeof(maxentmc_float_t)*(_s_)*((_s_)+1))

#define MAXENTMC_QUAD_THREAD_HEADER_SIZE MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,sizeof(struct maxentmc_quad_helper_thread_struct))

#define MAXENTMC_QUAD_THREAD_MOMENT_SIZE(_s_) MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,MAXENTMC_QUAD_THREAD_HEADER_SIZE+sizeof(maxentmc_float_t)*(_s_))

#define MAXENTMC_QUAD_THREAD_FULL_SIZE(_s_) MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,MAXENTMC_QUAD_THREAD_MOMENT_SIZE(_s_)+sizeof(maxentmc_float_t)*(_s_)*MAXENTMC_QUAD_THREAD_HOWMANY_AT_ONCE)

static pthread_mutex_t maxentmc_quad_helper_global_lock = PTHREAD_MUTEX_INITIALIZER;

#ifdef MAXENTMC_SINGLE_PRECISION
void ssyev_(char const * const, char const * const, unsigned int const * const, maxentmc_float_t * const,
            unsigned int const * const, maxentmc_float_t * const, maxentmc_float_t * const,
            int const * const, int * const);
#else
void dsyev_(char const * const, char const * const, unsigned int const * const, maxentmc_float_t * const,
            unsigned int const * const, maxentmc_float_t * const, maxentmc_float_t * const,
            int const * const, int * const);
#endif

struct maxentmc_quad_helper_struct * maxentmc_quad_helper_alloc(maxentmc_index_t const dim)
{
    if(dim == 0){
        MAXENTMC_MESSAGE(stderr,"error: zero dimension");
        return NULL;
    }

    int status;

    struct maxentmc_quad_helper_struct * q;

    MAXENTMC_ALLOC(q,MAXENTMC_QUAD_HELPER_SIZE(dim),status);

    if(status)
        return NULL;

    q->dimension = dim;

    q->max_power = 0;

    q->locked = 0;

    q->n_mult = 0;

    q->n_mom = 0;

    q->shift = MAXENTMC_INCREMENT_POINTER(q,MAXENTMC_QUAD_HELPER_HEADER_SIZE);

    q->rotate = q->shift + dim;

    q->multipliers = NULL;

    q->moments = NULL;

    q->moment_list = NULL;

    q->multiplier_list = NULL;

    maxentmc_quad_helper_set_shift_rotation(q,NULL);

    return q;
}

void maxentmc_quad_helper_free(struct maxentmc_quad_helper_struct * const q)
{
    if(q){
        pthread_mutex_lock(&maxentmc_quad_helper_global_lock);
        if(q->locked){
            MAXENTMC_MESSAGE(stderr,"error: quadrature helper is locked");
        }
        else{
            struct maxentmc_quad_helper_power_list_struct * temp = q->multiplier_list;
            while(temp){
                maxentmc_power_vector_free(temp->power_vector);
                struct maxentmc_quad_helper_power_list_struct * temp2 = temp;
                temp = temp2->next;
                free(temp2);
            }
            temp = q->moment_list;
            while(temp){
                maxentmc_power_vector_free(temp->power_vector);
                struct maxentmc_quad_helper_power_list_struct * temp2 = temp;
                temp = temp2->next;
                free(temp2);
            }
            free(q);
        }
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
    }
}

maxentmc_index_t maxentmc_quad_helper_get_dimension(struct maxentmc_quad_helper_struct const * const q)
{
    if(q)
        return q->dimension;
    else{
        MAXENTMC_MESSAGE(stderr,"error: NULL pointer provided");
        return 0;
    }
}

int maxentmc_quad_helper_set_shift_rotation(struct maxentmc_quad_helper_struct * const q, struct maxentmc_power_vector_struct * const constraints)
{
    MAXENTMC_CHECK_NULL(q);
    pthread_mutex_lock(&maxentmc_quad_helper_global_lock);
    if(q->locked){
        MAXENTMC_MESSAGE(stderr,"error: quadrature helper is locked");
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        return -1;
    }

    if(constraints){

        if(constraints->powers->dimension != q->dimension){
            MAXENTMC_MESSAGE(stderr,"error: dimensions of quadrature helper and constraints do not match");
            pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
            return -1;
        }

        /** Cycle through constraints, find the mean and covariance **/

        maxentmc_index_t const dim = q->dimension;
        maxentmc_float_t mean[dim], cov[dim*dim];
        size_t const size = constraints->gsl_vec.size;
        size_t i, c_s = 0, c_d = 0;

        do{
            maxentmc_index_t const * const p = constraints->powers->power[c_s];

            /** Compute the total power, see if it is the mean or the covariance **/

            size_t total_power = 0;
            for(i=0;i<dim;++i)
                total_power += p[i];

            if(total_power == 1){

                /** This is the mean power **/

                /** See where it is and copy the data into the mean **/

                i=0;
                while(!(p[i])) ++i;
                mean[i] = constraints->gsl_vec.data[c_s];

                ++c_d;
            }

            if(total_power == 2){

                /** This is the covariance power **/

                /** See if it is diagonal or off-diagonal element **/

                total_power = 0;

                for(i=0;i<dim;++i)
                    total_power = QUAD_MAX(total_power,p[i]);

                if(total_power == 2){
                    /** Diagonal element **/
                    i=0;
                    while(!(p[i])) ++i;
                    cov[i*(dim+1)] = constraints->gsl_vec.data[c_s];
                    ++c_d;
                }
                else{
                    /** Off-diagonal element **/
                    i=0;
                    while(!(p[i])) ++i;
                    size_t j = i+1;
                    while(!(p[j])) ++j;
                    cov[i*dim+j] = constraints->gsl_vec.data[c_s];
                    cov[j*dim+i] = cov[i*dim+j];
                    c_d += 2;
                }

            }

            ++c_s;

        }while((c_s<size) && (c_d<dim*(dim+1)));

        if(c_d<dim*(dim+1)){
            MAXENTMC_MESSAGE(stderr,"warning: could not extract shift-rotation data");
            pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
            return -1;
        }

        for(i=0;i<dim;++i){
            size_t j;
            for(j=0;j<dim;++j)
                cov[i*dim+j] -= mean[i]*mean[j];
        }


        /** DEBUG **/
/*
        puts("Quadrature mean and covariance");
        for(i=0;i<dim;++i){
            printf("%g |",mean[i]);
            size_t j;
            for(j=0;j<dim;++j)
                printf(" %g",cov[j*dim+i]);
            puts("");
        }
*/
        /** END DEBUG **/
        /** Compute the shift **/

        for(i=0;i<dim;++i)
            q->shift[i] = mean[i];

        /** Now need to compute the rotation **/

        /** This is Fortran garbage **/

        char const jobz = 'V';
        char const uplo = 'U';
        unsigned int const n = dim;
        unsigned int const lda = dim;
        int lwork = -1, info;

#ifdef MAXENTMC_SINGLE_PRECISION
        ssyev_(&jobz, &uplo, &n, NULL, &lda, NULL, mean, &lwork, &info);
        if(info){
            MAXENTMC_MESSAGE(stderr,"error: ssyev failed on query stage");
            pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
            return -1;
        }
#else
        dsyev_(&jobz, &uplo, &n, NULL, &lda, NULL, mean, &lwork, &info);
        if(info){
            MAXENTMC_MESSAGE(stderr,"error: dsyev failed on query stage");
            pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
            return -1;
        }
#endif

        lwork = mean[0]+0.1;
        maxentmc_float_t work[lwork];

#ifdef MAXENTMC_SINGLE_PRECISION
        ssyev_(&jobz, &uplo, &n, cov, &lda, mean, work, &lwork, &info);
        if(info){
            MAXENTMC_MESSAGE(stderr,"error: ssyev failed on computation stage");
            pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
            return -1;
        }
#else
        dsyev_(&jobz, &uplo, &n, cov, &lda, mean, work, &lwork, &info);
        if(info){
            MAXENTMC_MESSAGE(stderr,"error: dsyev failed on computation stage");
            pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
            return -1;
        }
#endif

        /** DEBUG **/
/*
        puts("Quadrature eigenvalues and eigenvectors");
        for(i=0;i<dim;++i){
            printf("%g |",mean[i]);
            size_t j;
            for(j=0;j<dim;++j)
                printf(" %g",cov[j*dim+i]);
            puts("");
        }
*/
        /** END DEBUG **/


        /** Check if eigenvalues are nonnegative **/

        for(i=0;i<dim;++i)
            if(mean[i]<=0){
                MAXENTMC_MESSAGE(stderr,"error: an eigenvalue is not positive");
                pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
                return -1;
            }

        q->scale = 1.0;
        for(i=0;i<dim;++i){
            maxentmc_float_t const e_temp = sqrt(mean[i]);
            q->scale *= e_temp;
            size_t j;
            for(j=0;j<dim;++j)
                q->rotate[j*dim+i] = cov[i*dim+j]*e_temp; /** Transposing on the fly because Fortran **/
        }

        q->shift_rotate = 1;

    }
    else{

        q->shift_rotate = 0;
        q->scale = 0;
        maxentmc_index_t const dim = q->dimension;
        memset(q->shift,0,sizeof(maxentmc_float_t)*dim);
        memset(q->rotate,0,sizeof(maxentmc_float_t)*dim*dim);

    }

    pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);

    return 0;
}

static struct maxentmc_power_vector_struct * maxentmc_quad_helper_find_power_vector(struct maxentmc_quad_helper_struct * const q, struct maxentmc_power_vector_struct const * const power_vector, maxentmc_index_t const which)
{

    struct maxentmc_quad_helper_power_list_struct * temp = (which)?q->moment_list:q->multiplier_list;

    if(temp){

        do{
            if(temp->power_vector->powers == power_vector->powers)
                return temp->power_vector;
            else{
                if(temp->next)
                    temp = temp->next;
                else{
                    temp->next = malloc(sizeof(struct maxentmc_quad_helper_power_list_struct));
                    temp = temp->next;
                    temp->power_vector = maxentmc_power_vector_alloc(power_vector);
                    temp->next = NULL;
                    if(which)
                        ++(q->n_mom);
                    else
                        ++(q->n_mult);
                    return temp->power_vector;
                }
            }

        }while(temp);

    }
    else{

        if(which){
            q->moment_list = malloc(sizeof(struct maxentmc_quad_helper_power_list_struct));
            temp = q->moment_list;
            ++(q->n_mom);
        }
        else{
            q->multiplier_list = malloc(sizeof(struct maxentmc_quad_helper_power_list_struct));
            temp = q->multiplier_list;
            ++(q->n_mult);
        }
        temp->power_vector = maxentmc_power_vector_alloc(power_vector);
        temp->next = NULL;
        return temp->power_vector;
    }

    return NULL;
}

int maxentmc_quad_helper_set_multipliers(struct maxentmc_quad_helper_struct * const q, struct maxentmc_power_vector_struct const * const power_vector)
{
    MAXENTMC_CHECK_NULL(q);
    MAXENTMC_CHECK_NULL(power_vector);
    pthread_mutex_lock(&maxentmc_quad_helper_global_lock);
    if(q->locked){
        MAXENTMC_MESSAGE(stderr,"error: quadrature helper is locked");
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        return -1;
    }

    if(q->dimension != power_vector->powers->dimension){
        MAXENTMC_MESSAGE(stderr,"error: dimensions do not match");
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        return -1;
    }

    struct maxentmc_power_vector_struct * temp_v = maxentmc_quad_helper_find_power_vector(q,power_vector,0);
    if(temp_v == NULL){
        MAXENTMC_MESSAGE(stderr,"error: maxentmc_quad_helper_find_power_vector returned NULL for some reason");
        return -1;
    }

/*
    memcpy(temp_v->gsl_vec.data, power_vector->gsl_vec.data, sizeof(maxentmc_float_t)*temp_v->gsl_vec.size);
*/
    size_t const stride = power_vector->gsl_vec.stride;
    size_t i;
    for(i=0;i<temp_v->gsl_vec.size;++i)
        temp_v->gsl_vec.data[i] = power_vector->gsl_vec.data[i*stride];

    q->multipliers = temp_v;

    if(q->moments && q->multipliers)
        q->max_power = QUAD_MAX(q->moments->powers->max_power,q->multipliers->powers->max_power);
    else{
        if(q->moments)
            q->max_power = q->moments->powers->max_power;
        else
            q->max_power = q->multipliers->powers->max_power;
    }

    pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);

    return 0;


}

int maxentmc_quad_helper_set_moments(struct maxentmc_quad_helper_struct * const q, struct maxentmc_power_vector_struct const * const power_vector)
{
    MAXENTMC_CHECK_NULL(q);
    MAXENTMC_CHECK_NULL(power_vector);
    pthread_mutex_lock(&maxentmc_quad_helper_global_lock);
    if(q->locked){
        MAXENTMC_MESSAGE(stderr,"error: quadrature helper is locked");
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        return -1;
    }

    if(q->dimension != power_vector->powers->dimension){
        MAXENTMC_MESSAGE(stderr,"error: dimensions do not match");
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        return -1;
    }

    struct maxentmc_power_vector_struct * temp_v = maxentmc_quad_helper_find_power_vector(q,power_vector,1);
    if(temp_v == NULL){
        MAXENTMC_MESSAGE(stderr,"error: maxentmc_quad_helper_find_power_vector returned NULL for some reason");
        return -1;
    }

    memset(temp_v->gsl_vec.data, 0, sizeof(maxentmc_float_t)*temp_v->gsl_vec.size);

    q->moments = temp_v;
    q->locked = 1;

    if(q->moments && q->multipliers)
        q->max_power = QUAD_MAX(q->moments->powers->max_power,q->multipliers->powers->max_power);
    else{
        if(q->moments)
            q->max_power = q->moments->powers->max_power;
        else
            q->max_power = q->multipliers->powers->max_power;
    }
    pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);

    return 0;
}

struct maxentmc_quad_helper_thread_struct * maxentmc_quad_helper_thread_alloc(struct maxentmc_quad_helper_struct * const q)
{

    MAXENTMC_CHECK_NULL_PT(q);

    pthread_mutex_lock(&maxentmc_quad_helper_global_lock);

    if(q->multipliers == NULL){
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        MAXENTMC_MESSAGE(stderr,"error: multipliers is NULL");
        return NULL;
    }

    if(q->moments == NULL){
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        MAXENTMC_MESSAGE(stderr,"error: moments is NULL");
        return NULL;
    }

    if(!q->locked){
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        MAXENTMC_MESSAGE(stderr,"error: quad helper not locked");
        return NULL;
    }

    int status;

    size_t size = sizeof(maxentmc_float_t)*q->moments->gsl_vec.size;

    struct maxentmc_quad_helper_thread_struct * qt;

    MAXENTMC_ALLOC(qt,MAXENTMC_QUAD_THREAD_FULL_SIZE(size),status);

    if(status)
        return NULL;

    qt->main_quadrature = q;

    qt->moments = MAXENTMC_INCREMENT_POINTER(qt,MAXENTMC_QUAD_THREAD_HEADER_SIZE);

    qt->scratch = MAXENTMC_INCREMENT_POINTER(qt,MAXENTMC_QUAD_THREAD_MOMENT_SIZE(size));

    memset(qt->moments,0,size);

    /** DEBUG **/
    /*
    MAXENTMC_MESSAGE_VARARG(stdout,"n_mult = %u, n_mom = %u",q->n_mult,q->n_mom);
    */
    /** END DEBUG **/

    pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);

    return qt;
}

int maxentmc_quad_helper_thread_compute(struct maxentmc_quad_helper_thread_struct * const qt, ...)
{
    MAXENTMC_CHECK_NULL(qt);
    maxentmc_index_t const dim = qt->main_quadrature->multipliers->powers->dimension;
    maxentmc_float_t x[dim], w;
    maxentmc_index_t i;
    va_list ap;
    va_start(ap,qt);
    for(i=0;i<dim;++i)
        x[i] = va_arg(ap,maxentmc_float_t);
    w = va_arg(ap,maxentmc_float_t);
    va_end(ap);
    return maxentmc_quad_helper_thread_compute_1(qt,x,w);
}

#define MAXENTMC_QUADRATURE_THREAD_COMPUTE(_N_)                                                 \
                                                                                                \
    struct maxentmc_power_vector_struct const * const multipliers = q->multipliers;             \
    struct maxentmc_power_vector_struct const * const moments = q->moments;                     \
                                                                                                \
    size_t const d_size = multipliers->gsl_vec.size;                                            \
    maxentmc_index_t const * const * const __restrict d_powers =                                \
        (maxentmc_index_t const * const * const)multipliers->powers->power;                     \
    maxentmc_float_t const * const __restrict d_data = multipliers->gsl_vec.data;               \
    maxentmc_float_t const scale = q->scale;                                                    \
                                                                                                \
    size_t m_size = moments->gsl_vec.size;                                                      \
    maxentmc_index_t const * const * const __restrict m_powers =                                \
        (maxentmc_index_t const * const * const)moments->powers->power;                         \
    maxentmc_float_t * const __restrict m_data = qt->moments;                                   \
                                                                                                \
    for(i=2;i<p1;++i){                                                                          \
        maxentmc_index_t j;                                                                     \
        for(j=0;j<dim;++j){                                                                     \
            maxentmc_float_t * __restrict _x1_ = x_pow[j*p1+1];                                 \
            maxentmc_float_t * __restrict _x_ = x_pow[j*p1+i];                                  \
            maxentmc_float_t * __restrict _x_prev_ = x_pow[j*p1+i-1];                           \
            maxentmc_index_t k;                                                                 \
            for(k=0;k<(_N_);++k)                                                                \
                _x_[k] = _x_prev_[k] * _x1_[k];                                                 \
        }                                                                                       \
    }                                                                                           \
                                                                                                \
    maxentmc_float_t rho[(_N_)];                                                                \
                                                                                                \
    memset(rho,0,sizeof(maxentmc_float_t)*(_N_));                                               \
                                                                                                \
    if(d_powers == m_powers){                                                                   \
                                                                                                \
        maxentmc_float_t * const temp_moments = qt->scratch;                                    \
                                                                                                \
        for(i=0;i<d_size;++i){                                                                  \
            maxentmc_index_t const * const __restrict d_p = d_powers[i];                        \
            maxentmc_float_t * const __restrict temp_m = temp_moments + i*(_N_);                \
            maxentmc_index_t j;                                                                 \
            for(j=0;j<(_N_);++j)                                                                \
                temp_m[j] = 1.0;                                                                \
            for(j=0;j<dim;++j){                                                                 \
                maxentmc_float_t const * const __restrict _x_ = x_pow[j*p1+d_p[j]];             \
                maxentmc_index_t k;                                                             \
                for(k=0;k<(_N_);++k)                                                            \
                    temp_m[k] *= _x_[k];                                                        \
            }                                                                                   \
            for(j=0;j<(_N_);++j)                                                                \
                rho[j] += d_data[i] * temp_m[j];                                                \
        }                                                                                       \
                                                                                                \
        for(i=0;i<(_N_);++i)                                                                    \
            rho[i] = exp(rho[i]) * w[i];                                                        \
                                                                                                \
        if(shift_rotate)                                                                        \
            for(i=0;i<(_N_);++i)                                                                \
                rho[i] *= scale;                                                                \
                                                                                                \
        for(i=0;i<m_size;++i){                                                                  \
            maxentmc_index_t j;                                                                 \
            maxentmc_float_t const * const __restrict temp_m = temp_moments + i*(_N_);          \
            for(j=0;j<(_N_);++j)                                                                \
                m_data[i] += temp_m[j]*rho[j];                                                  \
        }                                                                                       \
                                                                                                \
    }                                                                                           \
    else{                                                                                       \
                                                                                                \
        maxentmc_float_t * const __restrict temp_m = qt->scratch;                               \
                                                                                                \
        for(i=0;i<d_size;++i){                                                                  \
            maxentmc_index_t const * const __restrict d_p = d_powers[i];                        \
            maxentmc_index_t j;                                                                 \
            for(j=0;j<(_N_);++j)                                                                \
                temp_m[j] = 1.0;                                                                \
            for(j=0;j<dim;++j){                                                                 \
                maxentmc_float_t const * const __restrict _x_ = x_pow[j*p1+d_p[j]];             \
                maxentmc_index_t k;                                                             \
                for(k=0;k<(_N_);++k)                                                            \
                    temp_m[k] *= _x_[k];                                                        \
            }                                                                                   \
            for(j=0;j<(_N_);++j)                                                                \
                rho[j] += d_data[i]*temp_m[j];                                                  \
        }                                                                                       \
                                                                                                \
        for(i=0;i<(_N_);++i)                                                                    \
            rho[i] = exp(rho[i]) * w[i];                                                        \
                                                                                                \
        if(shift_rotate)                                                                        \
            for(i=0;i<(_N_);++i)                                                                \
                rho[i] *= scale;                                                                \
                                                                                                \
        for(i=0;i<m_size;++i){                                                                  \
            maxentmc_index_t const * const __restrict m_p = m_powers[i];                        \
            maxentmc_index_t j;                                                                 \
            for(j=0;j<(_N_);++j)                                                                \
                temp_m[j] = rho[j];                                                             \
            for(j=0;j<dim;++j){                                                                 \
                maxentmc_float_t const * const __restrict _x_ = x_pow[j*p1+m_p[j]];             \
                maxentmc_index_t k;                                                             \
                for(k=0;k<(_N_);++k)                                                            \
                    temp_m[k] *= _x_[k];                                                        \
            }                                                                                   \
            for(j=0;j<(_N_);++j)                                                                \
                m_data[i] += temp_m[j];                                                         \
        }                                                                                       \
                                                                                                \
    }


int maxentmc_quad_helper_thread_compute_1(struct maxentmc_quad_helper_thread_struct * const qt,
                                          maxentmc_float_t const * const x1, maxentmc_float_t const w1)
{
    MAXENTMC_CHECK_NULL(qt);
    MAXENTMC_CHECK_NULL(x1);

    struct maxentmc_quad_helper_struct const * const q = qt->main_quadrature;
    maxentmc_index_t const dim = q->multipliers->powers->dimension;
    maxentmc_index_t const p1 = q->max_power+1;
    maxentmc_index_t const shift_rotate = q->shift_rotate;
    maxentmc_float_t const * const __restrict shift = q->shift;
    maxentmc_float_t const * const __restrict rotate = q->rotate;

    maxentmc_float_t x_pow[dim*p1][1];

    maxentmc_float_t w[1];

    size_t i;

    for(i=0;i<dim;++i)
        x_pow[i*p1][0] = 1.0;

    if(shift_rotate){

        for(i=0;i<dim;++i){
            x_pow[i*p1+1][0] = shift[i];
            maxentmc_index_t j;
            for(j=0;j<dim;++j)
                x_pow[i*p1+1][0] += rotate[i*dim+j]*x1[j];
        }

    }
    else
        for(i=0;i<dim;++i)
            x_pow[i*p1+1][0] = x1[i];

    w[0] = w1;

    MAXENTMC_QUADRATURE_THREAD_COMPUTE(1);

    return 0;
}

int maxentmc_quad_helper_thread_compute_2(struct maxentmc_quad_helper_thread_struct * const qt,
                                          maxentmc_float_t const * const x1, maxentmc_float_t const w1,
                                          maxentmc_float_t const * const x2, maxentmc_float_t const w2)
{
    MAXENTMC_CHECK_NULL(qt);
    MAXENTMC_CHECK_NULL(x1);
    MAXENTMC_CHECK_NULL(x2);

    struct maxentmc_quad_helper_struct const * const q = qt->main_quadrature;
    maxentmc_index_t const dim = q->multipliers->powers->dimension;
    maxentmc_index_t const p1 = q->max_power+1;
    maxentmc_index_t const shift_rotate = q->shift_rotate;
    maxentmc_float_t const * const __restrict shift = q->shift;
    maxentmc_float_t const * const __restrict rotate = q->rotate;

    maxentmc_float_t x_pow[dim*p1][2];
    maxentmc_float_t w[2];

    size_t i;

    for(i=0;i<dim;++i){
        maxentmc_index_t j;
        for(j=0;j<2;++j)
            x_pow[i*p1][j] = 1.0;
    }

    if(shift_rotate){

        for(i=0;i<dim;++i){
            maxentmc_index_t j;
            for(j=0;j<2;++j)
                x_pow[i*p1+1][j] = shift[i];
            for(j=0;j<dim;++j){
                x_pow[i*p1+1][0] += rotate[i*dim+j]*x1[j];
                x_pow[i*p1+1][1] += rotate[i*dim+j]*x2[j];
            }
        }

    }
    else
        for(i=0;i<dim;++i){
            x_pow[i*p1+1][0] = x1[i];
            x_pow[i*p1+1][1] = x2[i];
        }


    w[0] = w1;
    w[1] = w2;

    MAXENTMC_QUADRATURE_THREAD_COMPUTE(2);

    return 0;
}

int maxentmc_quad_helper_thread_compute_3(struct maxentmc_quad_helper_thread_struct * const qt,
                                          maxentmc_float_t const * const x1, maxentmc_float_t const w1,
                                          maxentmc_float_t const * const x2, maxentmc_float_t const w2,
                                          maxentmc_float_t const * const x3, maxentmc_float_t const w3)
{
    MAXENTMC_CHECK_NULL(qt);
    MAXENTMC_CHECK_NULL(x1);
    MAXENTMC_CHECK_NULL(x2);
    MAXENTMC_CHECK_NULL(x3);

    struct maxentmc_quad_helper_struct const * const q = qt->main_quadrature;
    maxentmc_index_t const dim = q->multipliers->powers->dimension;
    maxentmc_index_t const p1 = q->max_power+1;
    maxentmc_index_t const shift_rotate = q->shift_rotate;
    maxentmc_float_t const * const __restrict shift = q->shift;
    maxentmc_float_t const * const __restrict rotate = q->rotate;

    maxentmc_float_t x_pow[dim*p1][3];
    maxentmc_float_t w[3];

    size_t i;

    for(i=0;i<dim;++i){
        maxentmc_index_t j;
        for(j=0;j<3;++j)
            x_pow[i*p1][j] = 1.0;
    }


    if(shift_rotate){

        for(i=0;i<dim;++i){
            maxentmc_index_t j;
            for(j=0;j<3;++j)
                x_pow[i*p1+1][j] = shift[i];
            for(j=0;j<dim;++j){
                x_pow[i*p1+1][0] += rotate[i*dim+j]*x1[j];
                x_pow[i*p1+1][1] += rotate[i*dim+j]*x2[j];
                x_pow[i*p1+1][2] += rotate[i*dim+j]*x3[j];
            }
        }

    }
    else
        for(i=0;i<dim;++i){
            x_pow[i*p1+1][0] = x1[i];
            x_pow[i*p1+1][1] = x2[i];
            x_pow[i*p1+1][2] = x3[i];
        }

    w[0] = w1;
    w[1] = w2;
    w[2] = w3;

    MAXENTMC_QUADRATURE_THREAD_COMPUTE(3);

    return 0;
}

int maxentmc_quad_helper_thread_compute_4(struct maxentmc_quad_helper_thread_struct * const qt,
                                          maxentmc_float_t const * const x1, maxentmc_float_t const w1,
                                          maxentmc_float_t const * const x2, maxentmc_float_t const w2,
                                          maxentmc_float_t const * const x3, maxentmc_float_t const w3,
                                          maxentmc_float_t const * const x4, maxentmc_float_t const w4)
{
    MAXENTMC_CHECK_NULL(qt);
    MAXENTMC_CHECK_NULL(x1);
    MAXENTMC_CHECK_NULL(x2);
    MAXENTMC_CHECK_NULL(x3);
    MAXENTMC_CHECK_NULL(x4);

    struct maxentmc_quad_helper_struct const * const q = qt->main_quadrature;
    maxentmc_index_t const dim = q->multipliers->powers->dimension;
    maxentmc_index_t const p1 = q->max_power+1;
    maxentmc_index_t const shift_rotate = q->shift_rotate;
    maxentmc_float_t const * const __restrict shift = q->shift;
    maxentmc_float_t const * const __restrict rotate = q->rotate;

    maxentmc_float_t x_pow[dim*p1][4];
    maxentmc_float_t w[4];

    size_t i;

    for(i=0;i<dim;++i){
        maxentmc_index_t j;
        for(j=0;j<4;++j)
            x_pow[i*p1][j] = 1.0;
    }

    if(shift_rotate){

        for(i=0;i<dim;++i){
            maxentmc_index_t j;
            for(j=0;j<4;++j)
                x_pow[i*p1+1][j] = shift[i];
            for(j=0;j<dim;++j){
                x_pow[i*p1+1][0] += rotate[i*dim+j]*x1[j];
                x_pow[i*p1+1][1] += rotate[i*dim+j]*x2[j];
                x_pow[i*p1+1][2] += rotate[i*dim+j]*x3[j];
                x_pow[i*p1+1][3] += rotate[i*dim+j]*x4[j];
            }
        }

    }
    else
        for(i=0;i<dim;++i){
            x_pow[i*p1+1][0] = x1[i];
            x_pow[i*p1+1][1] = x2[i];
            x_pow[i*p1+1][2] = x3[i];
            x_pow[i*p1+1][3] = x4[i];
        }

    w[0] = w1;
    w[1] = w2;
    w[2] = w3;
    w[3] = w4;

    MAXENTMC_QUADRATURE_THREAD_COMPUTE(4);

    return 0;
}


int maxentmc_quad_helper_thread_merge(struct maxentmc_quad_helper_thread_struct * const qt)
{
    MAXENTMC_CHECK_NULL(qt);

    pthread_mutex_lock(&maxentmc_quad_helper_global_lock);
    if(!qt->main_quadrature->locked){
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        MAXENTMC_MESSAGE(stderr,"error: quad helper not locked");
        return -1;
    }

    size_t i;

    for(i=0;i<qt->main_quadrature->moments->gsl_vec.size;++i)
        qt->main_quadrature->moments->gsl_vec.data[i] += qt->moments[i];

    pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);

    free(qt);

    return 0;
}

int maxentmc_quad_helper_get_moments(struct maxentmc_quad_helper_struct * const q, struct maxentmc_power_vector_struct * const moments)
{
    MAXENTMC_CHECK_NULL(q);
    MAXENTMC_CHECK_NULL(moments);
    pthread_mutex_lock(&maxentmc_quad_helper_global_lock);
    if(!q->locked){
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        MAXENTMC_MESSAGE(stderr,"error: quad helper not locked");
        return -1;
    }

    if(q->moments->powers != moments->powers){
        pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);
        MAXENTMC_MESSAGE(stderr,"error: powers do not match");
        return -1;
    }

/*
    memcpy(moments->gsl_vec.data,q->moments->gsl_vec.data,sizeof(maxentmc_float_t)*moments->gsl_vec.size);
*/
    size_t const stride = moments->gsl_vec.stride;
    size_t i;
    for(i=0;i<moments->gsl_vec.size;++i)
        moments->gsl_vec.data[i*stride] = q->moments->gsl_vec.data[i];

    q->locked = 0;

    pthread_mutex_unlock(&maxentmc_quad_helper_global_lock);

    return 0;
}
