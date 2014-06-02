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
#include <stdint.h>
#include <netinet/in.h>
#include "maxentmc_vector.h"

/** Auxiliary routines **/

/** Both ntoh and hton functions assume IEEE 754 floating point **/

static maxentmc_float_t maxentmc_float_t_ntoh(maxentmc_float_t x)
{
#if __BYTE_ORDER == __LITTLE_ENDIAN
    if(sizeof(maxentmc_float_t) == 8)
        return x;
    else if(sizeof(maxentmc_float_t) == 16){
        union {maxentmc_float_t x; uint16_t i;} v;
        v.x = x;
        v.i = __bswap_16(v.i);
        return v.x;
    }
    else if(sizeof(maxentmc_float_t) == 32){
        union {maxentmc_float_t x; uint32_t i;} v;
        v.x = x;
        v.i = __bswap_32(v.i);
        return v.x;
    }
    else if(sizeof(maxentmc_float_t) == 64){
        union {maxentmc_float_t x; uint64_t i;} v;
        v.x = x;
        v.i = __bswap_64(v.i);
        return v.x;
    }
    else{
        MAXENTMC_MESSAGE(stderr,"error: cannot byte swap this type");
        return x;
    }
#endif

#if __BYTE_ORDER == __BIG_ENDIAN
            return x;
#endif

#if __BYTE_ORDER == __PDP_ENDIAN
#error PDP byte order not supported
#endif

}

static maxentmc_float_t maxentmc_float_t_hton(maxentmc_float_t x)
{
#if __BYTE_ORDER == __LITTLE_ENDIAN
    if(sizeof(maxentmc_float_t) == 8)
        return x;
    else if(sizeof(maxentmc_float_t) == 16){
        union {maxentmc_float_t x; uint16_t i;} v;
        v.x = x;
        v.i = __bswap_16(v.i);
        return v.x;
    }
    else if(sizeof(maxentmc_float_t) == 32){
        union {maxentmc_float_t x; uint32_t i;} v;
        v.x = x;
        v.i = __bswap_32(v.i);
        return v.x;
    }
    else if(sizeof(maxentmc_float_t) == 64){
        union {maxentmc_float_t x; uint64_t i;} v;
        v.x = x;
        v.i = __bswap_64(v.i);
        return v.x;
    }
    else{
        MAXENTMC_MESSAGE(stderr,"error: cannot byte swap this type");
        return x;
    }
#endif

#if __BYTE_ORDER == __BIG_ENDIAN
            return x;
#endif

#if __BYTE_ORDER == __PDP_ENDIAN
#error PDP byte order not supported
#endif

}

/** Block structure definition **/

#ifndef __GSL_BLOCK_H__
struct gsl_block_struct {

    size_t size;
    maxentmc_float_t * data;

};
typedef struct gsl_block_struct maxentmc_gsl_block_t;
#else
#ifdef MAXENTMC_SINGLE_PRECISION
typedef gsl_block_float maxentmc_gsl_block_t;
#else
typedef gsl_block maxentmc_gsl_block_t;
#endif
#endif

/** Vector routines **/

#define MAXENTMC_VECTOR_SUB_HEADER_SIZE     MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_gsl_block_t),sizeof(struct maxentmc_power_vector_struct))
#define MAXENTMC_VECTOR_HEADER_SIZE         MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,MAXENTMC_VECTOR_SUB_HEADER_SIZE+sizeof(maxentmc_gsl_block_t))
#define MAXENTMC_VECTOR_SIZE(_s_)           MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,MAXENTMC_VECTOR_HEADER_SIZE+sizeof(maxentmc_float_t)*(_s_))

size_t maxentmc_power_vector_size(size_t const s)
{
    return MAXENTMC_VECTOR_SIZE(s);
}

int maxentmc_power_vector_init(struct maxentmc_power_struct const * const p, struct maxentmc_power_vector_struct * const v)
{
    MAXENTMC_CHECK_NULL(p);
    MAXENTMC_CHECK_NULL(v);

    v->gsl_vec.size = p->size;
    v->gsl_vec.stride = 1;
    v->gsl_vec.data = MAXENTMC_INCREMENT_POINTER(v,MAXENTMC_VECTOR_HEADER_SIZE);
    v->gsl_vec.block = MAXENTMC_INCREMENT_POINTER(v,MAXENTMC_VECTOR_SUB_HEADER_SIZE);
    v->gsl_vec.block->size = v->gsl_vec.size;
    v->gsl_vec.block->data = v->gsl_vec.data;
    v->gsl_vec.owner = 0;
    v->powers = (struct maxentmc_power_struct *)p;

    return 0;
}

maxentmc_index_t maxentmc_power_vector_get_dimension(struct maxentmc_power_vector_struct const * const v)
{
    return (v)?v->powers->dimension:0;
}

int maxentmc_power_vector_get_max_power(struct maxentmc_power_vector_struct const * const v, maxentmc_index_t * const max_power)
{
    MAXENTMC_CHECK_NULL(v);
    MAXENTMC_CHECK_NULL(max_power);
    *max_power = v->powers->max_power;
    return 0;
}

struct maxentmc_power_vector_struct * maxentmc_power_vector_alloc_from_power(struct maxentmc_power_struct * const dp)
{
    MAXENTMC_CHECK_NULL_PT(dp);

    struct maxentmc_power_vector_struct * d;
    int status;

    MAXENTMC_ALLOC(d,MAXENTMC_VECTOR_SIZE(dp->size),status);
    if(status)
        return NULL;

    if(maxentmc_power_vector_init(dp,d)){
        free(d);
        return NULL;
    }
    ++(d->powers->num_refs);

    return d;
}

struct maxentmc_power_vector_struct * maxentmc_power_vector_alloc(struct maxentmc_power_vector_struct const * const d)
{
    MAXENTMC_CHECK_NULL_PT(d);
    return maxentmc_power_vector_alloc_from_power(d->powers);
}

struct maxentmc_power_vector_struct * maxentmc_power_vector_product_alloc(struct maxentmc_power_vector_struct const * const d1, struct maxentmc_power_vector_struct const * const d2)
{
    MAXENTMC_CHECK_NULL_PT(d1);
    MAXENTMC_CHECK_NULL_PT(d2);

    struct maxentmc_power_struct * dp = maxentmc_power_alloc_product(d1->powers,d2->powers);

    MAXENTMC_CHECK_NULL_PT(dp);

    struct maxentmc_power_vector_struct * d = maxentmc_power_vector_alloc_from_power(dp);

    if(d==NULL)
        maxentmc_power_free(dp);

    return d;
}

void maxentmc_power_vector_free(struct maxentmc_power_vector_struct * const d)
{
    if(d){
        maxentmc_power_free(d->powers);
        free(d);
    }
}

int maxentmc_power_vector_print(struct maxentmc_power_vector_struct const * const d, FILE * const out)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(out);

    fprintf(out,"dimension = %u, size = %zu, max_power = %u,",d->powers->dimension,d->powers->size,d->powers->max_power);

    if(MAXENTMC_POWER_CHECK_PROPERTY(MAXENTMC_POWER_COMPLETE,d->powers))
        fputs(" complete,",out);
    else
        fputs(" incomplete,",out);
    if(MAXENTMC_POWER_CHECK_PROPERTY(MAXENTMC_POWER_ORDERED,d->powers))
        fputs(" ordered\n",out);
    else
        fputs(" unordered\n",out);

    size_t i;
    for(i=0;i<d->powers->size;++i){
        maxentmc_index_t j;
        for(j=0;j<d->powers->dimension;++j)
            fprintf(out,"%u ",d->powers->power[i][j]);
        fprintf(out,"| %g\n",d->gsl_vec.data[i]);
    }

    /** DEBUG **/
/*
    maxentmc_power_print_product(d->powers,out);
*/

    return 0;
}

int maxentmc_power_vector_find_element(struct maxentmc_power_vector_struct const * const v, ...)
{
    MAXENTMC_CHECK_NULL(v);
    maxentmc_index_t const dim = v->powers->dimension;
    maxentmc_index_t p[dim], i;
    size_t * pos;
    va_list ap;
    va_start(ap,v);
    for(i=0;i<dim;++i)
        p[i] = va_arg(ap,int);
    pos = va_arg(ap,size_t *);
    va_end(ap);
    return maxentmc_power_find(v->powers,p,pos);
}

int maxentmc_power_vector_find_element_ca(struct maxentmc_power_vector_struct const * const v, maxentmc_index_t const * const p, size_t * const pos)
{
    MAXENTMC_CHECK_NULL(v);
    return maxentmc_power_find(v->powers,p,pos);
}

int maxentmc_power_vector_get_powers_ca(struct maxentmc_power_vector_struct const * const d, size_t const pos, maxentmc_index_t * const p)
{
    MAXENTMC_CHECK_NULL(d);
    return maxentmc_power_get_power(d->powers,pos,p);
}

int maxentmc_power_vector_get_powers(struct maxentmc_power_vector_struct const * const d, size_t const pos, ...)
{
    MAXENTMC_CHECK_NULL(d);
    maxentmc_index_t const dim = d->powers->dimension;
    maxentmc_index_t p[dim];
    int status = maxentmc_power_get_power(d->powers,pos,p);
    if(!status){
        va_list ap;
        va_start(ap,pos);
        maxentmc_index_t i;
        for(i=0;i<dim;++i)
            *va_arg(ap,maxentmc_index_t *) = p[i];
        va_end(ap);
    }
    return status;
}

/*
int maxentmc_power_vector_power_multiply_add(maxentmc_vector_t const * const v1, maxentmc_vector_t const * const v2, maxentmc_vector_t * const v)
{

    MAXENTMC_CHECK_NULL(v1);
    MAXENTMC_CHECK_NULL(v2);
    MAXENTMC_CHECK_NULL(v);

    return maxentmc_power_multiply_add(1,v1->powers,0,(maxentmc_float_t const * const * const)(&(v1->data)),
                                       v2->powers,0,(maxentmc_float_t const * const * const)(&(v2->data)),
                                       v->powers,0,(maxentmc_float_t * const * const)(&(v->data)));

}
*/

struct maxentmc_power_vector_struct * maxentmc_power_vector_fread_power(FILE * const in)
{
    return maxentmc_power_vector_alloc_from_power(maxentmc_power_fread(in));
}

int maxentmc_power_vector_fwrite_power(struct maxentmc_power_vector_struct const * const v, FILE * const out)
{
    return maxentmc_power_fwrite(v->powers,out);
}

int maxentmc_power_vector_fread_values(struct maxentmc_power_vector_struct * const v, FILE * const in)
{
    MAXENTMC_CHECK_NULL(v);
    MAXENTMC_CHECK_NULL(in);

    size_t const size = v->gsl_vec.size, stride = v->gsl_vec.stride;
    maxentmc_float_t temp[size];

    MAXENTMC_FREAD(temp,sizeof(maxentmc_float_t)*size,in);

    size_t i;
    for(i=0;i<size;++i)
        v->gsl_vec.data[i*stride] = maxentmc_float_t_ntoh(temp[i]);

    return 0;
}

int maxentmc_power_vector_fwrite_values(struct maxentmc_power_vector_struct const * const v, FILE * const out)
{
    MAXENTMC_CHECK_NULL(v);
    MAXENTMC_CHECK_NULL(out);

    size_t const size = v->gsl_vec.size, stride = v->gsl_vec.stride;
    maxentmc_float_t temp[size];

    size_t i;
    for(i=0;i<size;++i)
        temp[i] = maxentmc_float_t_hton(v->gsl_vec.data[i*stride]);

    MAXENTMC_FWRITE(temp,sizeof(maxentmc_float_t)*size,out);

    return 0;
}

int maxentmc_power_vector_compute_polynomial(struct maxentmc_power_vector_struct const * const v, ...)
{
    MAXENTMC_CHECK_NULL(v);
    maxentmc_index_t const dim = v->powers->dimension;
    maxentmc_float_t x[dim], *result;
    va_list ap;
    va_start(ap,v);
    maxentmc_index_t i;
    for(i=0;i<dim;++i)
        x[i] = va_arg(ap,maxentmc_float_t);
    result = va_arg(ap,maxentmc_float_t *);
    va_end(ap);
    return maxentmc_power_vector_compute_polynomial_ca(v,x,result);
}

int maxentmc_power_vector_compute_polynomial_ca(struct maxentmc_power_vector_struct const * const v,
                                                          maxentmc_float_t const * const x,
                                                          maxentmc_float_t * const result)
{
    MAXENTMC_CHECK_NULL(v);
    MAXENTMC_CHECK_NULL(x);
    MAXENTMC_CHECK_NULL(result);

    maxentmc_index_t const dim = v->powers->dimension;
    maxentmc_index_t const p1 = v->powers->max_power+1;
    size_t const size = v->gsl_vec.size;
    maxentmc_index_t const * const * const __restrict powers =
        (maxentmc_index_t const * const * const)v->powers->power;
    maxentmc_float_t const * const __restrict data = v->gsl_vec.data;

    maxentmc_float_t x_pow[dim*p1];

    size_t i;

    for(i=0;i<dim;++i){
        x_pow[i*p1] = 1.0;
        x_pow[i*p1+1] = x[i];
    }

    for(i=2;i<p1;++i){
        maxentmc_index_t j;
        for(j=0;j<dim;++j)
            x_pow[j*p1+i] = x_pow[j*p1+i-1] * x_pow[j*p1+1];
    }

    maxentmc_float_t res = 0.0;

    for(i=0;i<size;++i){
        maxentmc_index_t const * const p = powers[i];
        maxentmc_float_t temp_moments = 1.0;
        maxentmc_index_t j;
        for(j=0;j<dim;++j)
            temp_moments *= x_pow[j*p1+p[j]];
        res += data[i]*temp_moments;
    }
    return 0;
}
