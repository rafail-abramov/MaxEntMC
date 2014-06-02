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
#include <netinet/in.h>
#include "maxentmc_power.h"

#define MAXENTMC_POWER_SUB_HEADER_SIZE MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_index_t *),sizeof(struct maxentmc_power_struct))
#define MAXENTMC_POWER_HEADER_SIZE(s) MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_index_t),MAXENTMC_POWER_SUB_HEADER_SIZE+sizeof(maxentmc_index_t *)*(s))
#define MAXENTMC_POWER_MAX_POWER_SIZE(s,d) MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_index_t),MAXENTMC_POWER_HEADER_SIZE(s)+sizeof(maxentmc_index_t)*(d))
#define MAXENTMC_POWER_SIZE(s,d) MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,MAXENTMC_POWER_MAX_POWER_SIZE(s,d)+sizeof(maxentmc_index_t)*(s)*(d))

#define MAXENTMC_MAX(a,b)  ((a)>(b))?(a):(b)

/** Auxiliary routines **/

size_t maxentmc_bincoeff(size_t const n, size_t const k)
{
    size_t i, b[n+1];
    b[0] = 1;
    for (i=1; i<=n; ++i){
        b[i] = 1;
        size_t j;
        for(j=i-1; j; --j)
            b[j] += b[j-1];
    }
    return b[k];
}

int maxentmc_power_comparison(maxentmc_index_t const dimension, maxentmc_index_t const * const p1, maxentmc_index_t const * const p2)
{

    /** First see how the sums are related **/

    size_t i, sum_p1=0, sum_p2=0;
    maxentmc_float_t m1 = 0.0, m2 = 0.0, v1 = 0.0, v2 = 0.0;
    for(i=0;i<dimension;++i){
        sum_p1 += p1[i];
        sum_p2 += p2[i];
        maxentmc_float_t x = (maxentmc_float_t)(p1[i]);
        m1 += x*i;
        v1 += (x*i)*i;
        x = (maxentmc_float_t)(p2[i]);
        m2 += x*i;
        v2 += (x*i)*i;
    }

    /** Sum test **/

    if(sum_p1 < sum_p2)
        return -1;
    if(sum_p1 > sum_p2)
        return 1;

    m1 /= sum_p1;
    m2 /= sum_p1;
    v1 /= sum_p1;
    v2 /= sum_p1;
    v1 -= m1*m1;
    v2 -= m2*m2;

    /** Variance test **/

    if(v1 < v2)
        return -1;
    if(v1 > v2)
        return 1;

    /** Mean test **/

    if(m1 < m2)
        return -1;
    if(m1 > m2)
        return 1;

    /** Guaranteed test **/

    i=dimension-1;
    size_t i1 = p1[i], i2 = p2[i];
    while(i){
        --i;
        i1 *= sum_p1+1;
        i2 *= sum_p1+1;
        i1 += p1[i];
        i2 += p2[i];
    }

    if(i1 < i2)
        return -1;
    if(i1 > i2)
        return 1;

    /** The powers are equal **/

    return 0;

}

static maxentmc_index_t maxentmc_index_t_ntoh(maxentmc_index_t x)
{
    switch(sizeof(maxentmc_index_t)){
        case sizeof(uint8_t):
            return x;
            break;
        case sizeof(uint16_t):
            return ntohs(x);
            break;
        case sizeof(uint32_t):
            return ntohl(x);
            break;
        case sizeof(uint64_t):
#if __BYTE_ORDER == __LITTLE_ENDIAN
            return __bswap_64(x);
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
            return x;
#endif
#if __BYTE_ORDER == __PDP_ENDIAN
#error PDP byte order not supported
#endif
            break;
        default:
            MAXENTMC_MESSAGE(stderr,"error: cannot byte swap this type");
            return 0;
    }
}

static size_t maxentmc_size_t_ntoh(size_t x)
{
    switch(sizeof(size_t)){
        case sizeof(uint8_t):
            return x;
            break;
        case sizeof(uint16_t):
            return ntohs(x);
            break;
        case sizeof(uint32_t):
            return ntohl(x);
            break;
        case sizeof(uint64_t):
#if __BYTE_ORDER == __LITTLE_ENDIAN
            return __bswap_64(x);
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
            return x;
#endif
#if __BYTE_ORDER == __PDP_ENDIAN
#error PDP byte order not supported
#endif
            break;
        default:
            MAXENTMC_MESSAGE(stderr,"error: cannot byte swap this type");
            return 0;
    }
}

static maxentmc_index_t maxentmc_index_t_hton(maxentmc_index_t x)
{
    switch(sizeof(maxentmc_index_t)){
        case sizeof(uint8_t):
            return x;
            break;
        case sizeof(uint16_t):
            return htons(x);
            break;
        case sizeof(uint32_t):
            return htonl(x);
            break;
        case sizeof(uint64_t):
#if __BYTE_ORDER == __LITTLE_ENDIAN
            return __bswap_64(x);
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
            return x;
#endif
#if __BYTE_ORDER == __PDP_ENDIAN
#error PDP byte order not supported
#endif
            break;
        default:
            MAXENTMC_MESSAGE(stderr,"error: cannot byte swap this type");
            return 0;
    }
}

static size_t maxentmc_size_t_hton(size_t x)
{
    switch(sizeof(size_t)){
        case sizeof(uint8_t):
            return x;
            break;
        case sizeof(uint16_t):
            return htons(x);
            break;
        case sizeof(uint32_t):
            return htonl(x);
            break;
        case sizeof(uint64_t):
#if __BYTE_ORDER == __LITTLE_ENDIAN
            return __bswap_64(x);
#endif
#if __BYTE_ORDER == __BIG_ENDIAN
            return x;
#endif
#if __BYTE_ORDER == __PDP_ENDIAN
#error PDP byte order not supported
#endif
            break;
        default:
            MAXENTMC_MESSAGE(stderr,"error: cannot byte swap this type");
            return 0;
    }
}

/** Power routines **/

int maxentmc_power_find(struct maxentmc_power_struct const * const d, maxentmc_index_t const * const p, size_t * const pos)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(p);
    MAXENTMC_CHECK_NULL(pos);

    size_t const sdim = sizeof(maxentmc_index_t)*d->dimension;
    maxentmc_index_t const * const * const pow = (maxentmc_index_t const * const * const)(d->power);
    size_t const size = d->size;
    size_t i=0;

    while((i<size) && memcmp(pow[i],p,sdim))
        ++i;

    if(i==size)
        return -1;

    *pos = i;

    return 0;
}

int maxentmc_power_get_power(struct maxentmc_power_struct const * const d, size_t const pos, maxentmc_index_t * const p)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(p);

    if(pos>=d->size)
        return -1;

    memcpy(p,d->power[pos],sizeof(maxentmc_index_t)*d->dimension);

    return 0;
}

struct maxentmc_power_struct * maxentmc_power_alloc(maxentmc_index_t const dimension, size_t const size)
{

    if(!dimension){
        MAXENTMC_MESSAGE(stderr,"error: dimension is zero");
        return NULL;
    }
    if(!size){
        MAXENTMC_MESSAGE(stderr,"error: size is zero");
        return NULL;
    }

    struct maxentmc_power_struct * d;
    int status;

    MAXENTMC_ALLOC(d,MAXENTMC_POWER_SIZE(size,dimension),status);
    if(status)
        return NULL;

    /** Partition the space **/

    d->properties = 0;
    d->dimension = dimension;
    d->max_power = 0;
    d->num_refs = 0;
    d->size = size;
    d->power = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_SUB_HEADER_SIZE);
    d->max_power_per_dimension = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_HEADER_SIZE(size));
    d->power[0] = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_MAX_POWER_SIZE(size,dimension));
    size_t i;
    for(i=1;i<size;++i)
        d->power[i] = d->power[i-1]+dimension;
    for(i=0;i<dimension;++i)
        d->max_power_per_dimension[i] = 0;
    d->product = NULL;

    /** done **/

    return d;

}

void maxentmc_power_free(struct maxentmc_power_struct * const d)
{
    if(d){
        if(d->num_refs > 1)
            --(d->num_refs);
        else
            free(d);
    }
}

int maxentmc_power_update_max_power(struct maxentmc_power_struct * const p)
{
    MAXENTMC_CHECK_NULL(p);

    maxentmc_index_t const dimension = p->dimension;
    maxentmc_index_t * const max_power_per_dim = p->max_power_per_dimension;
    size_t const size = p->size;
    size_t i;

    for(i=0;i<dimension;++i)
        max_power_per_dim[i] = 0;

    for(i=0;i<size;++i){
        maxentmc_index_t k;
        for(k=0;k<dimension;++k)
            max_power_per_dim[k] = MAXENTMC_MAX(max_power_per_dim[k],p->power[i][k]);
    }

    p->max_power = 0;
    for(i=0;i<dimension;++i)
        p->max_power = MAXENTMC_MAX(p->max_power,max_power_per_dim[i]);

    if(maxentmc_bincoeff(dimension+p->max_power,p->max_power) == size)
        MAXENTMC_POWER_SET_PROPERTY(MAXENTMC_POWER_COMPLETE,p);
    else
        MAXENTMC_POWER_CLEAR_PROPERTY(MAXENTMC_POWER_COMPLETE,p);

    return 0;
}

int maxentmc_power_print(struct maxentmc_power_struct const * const d, FILE * const out)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(out);

    fprintf(out,"dimension = %u, size = %zu, max_power = %u, num_refs = %u,",d->dimension,d->size,d->max_power,d->num_refs);

    if(MAXENTMC_POWER_CHECK_PROPERTY(MAXENTMC_POWER_COMPLETE,d))
        fputs(" complete,",out);
    else
        fputs(" incomplete,",out);
    if(MAXENTMC_POWER_CHECK_PROPERTY(MAXENTMC_POWER_ORDERED,d))
        fputs(" ordered\n",out);
    else
        fputs(" unordered\n",out);

    size_t i;
    for(i=0;i<d->size;++i){
        maxentmc_index_t j;
        for(j=0;j<d->dimension;++j)
            fprintf(out," %u",d->power[i][j]);
        fputs("\n",out);
    }

    return 0;

}

struct maxentmc_power_struct * maxentmc_power_fread(FILE * const in)
{
    MAXENTMC_CHECK_NULL_PT(in);

    /** First, read the header **/

    char header[MAXENTMC_POWER_STREAM_HEADER_SIZE];
    MAXENTMC_FREAD_PT(header,MAXENTMC_POWER_STREAM_HEADER_SIZE,in);
    if(memcmp(header,MAXENTMC_POWER_STREAM_HEADER,MAXENTMC_POWER_STREAM_HEADER_SIZE)){
        MAXENTMC_MESSAGE(stderr,"error: stream did not contain a proper header");
        return NULL;
    }

    /** Now read properties, dimension and size **/

    uint8_t properties;
    MAXENTMC_FREAD_PT(&properties,sizeof(uint8_t),in);

    maxentmc_index_t dimension;
    MAXENTMC_FREAD_PT(&dimension,sizeof(maxentmc_index_t),in);
    dimension = maxentmc_index_t_ntoh(dimension);

    size_t size;
    MAXENTMC_FREAD_PT(&size,sizeof(size_t),in);
    size = maxentmc_size_t_ntoh(size);

    /** Read powers **/

    maxentmc_index_t p[size*dimension];
    MAXENTMC_FREAD_PT(p,sizeof(maxentmc_index_t)*size*dimension,in);

    /** Everything is loaded successfully, now need to create the power structure **/

    struct maxentmc_power_struct * power = maxentmc_power_alloc(dimension,size);
    power->properties = properties;
    size_t i;
    for(i=0;i<size;++i){
        maxentmc_index_t j;
        for(j=0;j<dimension;++j)
            power->power[i][j] = maxentmc_index_t_ntoh(p[i*dimension+j]);
    }

    maxentmc_power_update_max_power(power);

    return power;
}

int maxentmc_power_fwrite(struct maxentmc_power_struct const * const power, FILE * const out)
{
    MAXENTMC_CHECK_NULL(power);
    MAXENTMC_CHECK_NULL(out);

    MAXENTMC_FWRITE(MAXENTMC_POWER_STREAM_HEADER,MAXENTMC_POWER_STREAM_HEADER_SIZE,out);

    MAXENTMC_FWRITE(&(power->properties),sizeof(uint8_t),out);

    maxentmc_index_t const dimension = power->dimension;
    maxentmc_index_t const n_dimension = maxentmc_index_t_hton(dimension);
    MAXENTMC_FWRITE(&n_dimension,sizeof(maxentmc_index_t),out);

    size_t const size = power->size;
    size_t const n_size = maxentmc_size_t_hton(size);
    MAXENTMC_FWRITE(&n_size,sizeof(size_t),out);

    maxentmc_index_t p[size*dimension];
    size_t i;
    for(i=0;i<size;++i){
        maxentmc_index_t j;
            for(j=0;j<dimension;++j)
                p[i*dimension+j] = maxentmc_index_t_hton(power->power[i][j]);
    }

    MAXENTMC_FWRITE(p,sizeof(maxentmc_index_t)*size*dimension,out);

    return 0;
}






/** The functions for the product of two polynomials **/

struct maxentmc_product_list_link_struct {

    size_t i1, i2;
    struct maxentmc_product_list_link_struct * next;

};

struct maxentmc_product_sparse_array_struct {

    maxentmc_index_t * power;
    struct maxentmc_product_list_link_struct * product_data;
    size_t product_data_size;

};

struct maxentmc_product_ordered_list_struct {

    maxentmc_index_t * power;
    struct maxentmc_product_list_link_struct * product_data;
    size_t product_data_size;
    struct maxentmc_product_ordered_list_struct * next;

};

struct maxentmc_power_struct * maxentmc_power_alloc_product(struct maxentmc_power_struct const * const p1, struct maxentmc_power_struct const * const p2)
{

    if((p1==NULL) || (p2==NULL)){
        MAXENTMC_MESSAGE(stderr,"error: NULL pointer provided");
        return NULL;
    }
    if((p1)->dimension != (p2)->dimension){
        MAXENTMC_MESSAGE(stderr,"error: dimension does not match");
        return NULL;
    }

    maxentmc_index_t const dimension = p1->dimension;
    maxentmc_index_t max_power_per_dim[dimension];
    size_t i;

    for(i=0;i<dimension;++i)
        max_power_per_dim[i] = p1->max_power_per_dimension[i]+p2->max_power_per_dimension[i];

    size_t prod_power_size = 1;
    for(i=0;i<dimension;++i)
        prod_power_size *= max_power_per_dim[i]+1;

    struct maxentmc_product_sparse_array_struct ** temp_product_powers = malloc(sizeof(struct maxentmc_product_sparse_array_struct *)*prod_power_size);
    memset(temp_product_powers,0,sizeof(struct maxentmc_product_sparse_array_struct *)*prod_power_size);

    for(i=0;i<p1->size;++i){
        size_t j;
        for(j=0;j<p2->size;++j){
            maxentmc_index_t k, power_per_dim[dimension];
            for(k=0;k<dimension;++k)
                power_per_dim[k] = p1->power[i][k] + p2->power[j][k];
            size_t c = power_per_dim[0];
            for(k=1;k<dimension;++k){
                c *= max_power_per_dim[k]+1;
                c += power_per_dim[k];
            }

            /** Now add the list at c **/

            if(temp_product_powers[c] == NULL){
                temp_product_powers[c] = malloc(sizeof(struct maxentmc_product_sparse_array_struct));
                temp_product_powers[c]->power = malloc(sizeof(maxentmc_index_t)*dimension);
                for(k=0;k<dimension;++k)
                    temp_product_powers[c]->power[k] = power_per_dim[k];
                temp_product_powers[c]->product_data = NULL;
                temp_product_powers[c]->product_data_size = 0;
            }

            struct maxentmc_product_list_link_struct * tp = malloc(sizeof(struct maxentmc_product_list_link_struct));

            tp->i1 = i;
            tp->i2 = j;

            tp->next = temp_product_powers[c]->product_data;
            temp_product_powers[c]->product_data = tp;
            ++(temp_product_powers[c]->product_data_size);
        }
    }

    /** The sparse list is built, now need to build the ordered list and destroy the sparse list in the process **/

    size_t size = 0, num_links = 0;

    struct maxentmc_product_ordered_list_struct * ordered_power_list = NULL;

    for(i=0;i<prod_power_size;++i)
        if(temp_product_powers[i]){
            ++size;
            num_links += temp_product_powers[i]->product_data_size;

            maxentmc_index_t * tpow = temp_product_powers[i]->power;
            struct maxentmc_product_list_link_struct * tpd = temp_product_powers[i]->product_data;
            size_t tpds = temp_product_powers[i]->product_data_size;
            free(temp_product_powers[i]);

            /** Here insert the data from the temp_product_powers into the ordered_power_list **/

            struct maxentmc_product_ordered_list_struct * topl = ordered_power_list;

            if(topl){

                /** Find the right position, malloc it **/

                while(topl->next && maxentmc_power_comparison(dimension, tpow, topl->next->power)>0)
                    topl = topl->next;

                struct maxentmc_product_ordered_list_struct * topl2 = malloc(sizeof(struct maxentmc_product_ordered_list_struct));
                topl2->next = topl->next;
                topl->next = topl2;
                topl = topl2;

            }
            else{
                /** This is the first link **/

                ordered_power_list = malloc(sizeof(struct maxentmc_product_ordered_list_struct));
                ordered_power_list->next = NULL;
                topl = ordered_power_list;

            }

            topl->power = tpow;
            topl->product_data = tpd;
            topl->product_data_size = tpds;

        }

    free(temp_product_powers);

    /** Ordered list is built, now need to build power structure **/

    /** Need to determine the size of the product polynomial **/

    struct maxentmc_power_struct * d;

    size_t const product_sub_header_length = MAXENTMC_ALIGNED_SIZE(sizeof(struct maxentmc_product_link_struct),
                                                                 MAXENTMC_POWER_SIZE(size,dimension)+sizeof(struct maxentmc_product_struct));
    size_t const product_header_length = MAXENTMC_ALIGNED_SIZE(sizeof(struct maxentmc_product_sublink_struct),
                                                             product_sub_header_length+sizeof(struct maxentmc_product_link_struct)*size);
    size_t const product_length = MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,product_header_length
                                                      +sizeof(struct maxentmc_product_sublink_struct)*num_links);

    int status;

    MAXENTMC_ALLOC(d,product_length,status);
    if(status)
        return NULL;

    /** Partition the space **/

    d->properties = 0;
    d->dimension = dimension;
    d->max_power = 0;
    d->num_refs = 0;
    d->size = size;
    d->power = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_SUB_HEADER_SIZE);
    d->max_power_per_dimension = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_HEADER_SIZE(size));
    d->power[0] = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_MAX_POWER_SIZE(size,dimension));
    for(i=1;i<size;++i)
        d->power[i] = d->power[i-1]+dimension;
    for(i=0;i<dimension;++i)
        d->max_power_per_dimension[i] = 0;

    d->product = MAXENTMC_INCREMENT_POINTER(d,MAXENTMC_POWER_SIZE(size,dimension));
    d->product->p1 = (struct maxentmc_power_struct *)p1;
    d->product->p2 = (struct maxentmc_power_struct *)p2;
    d->product->entry = MAXENTMC_INCREMENT_POINTER(d,product_sub_header_length);
    d->product->entry[0].p = MAXENTMC_INCREMENT_POINTER(d,product_header_length);

    /** Copy the data and disassemble the ordered list in the process **/

    for(i=0;i<size;++i){
        d->product->entry[i].np = ordered_power_list->product_data_size;
        if(i)
            d->product->entry[i].p = d->product->entry[i-1].p + d->product->entry[i-1].np;
        size_t j;
        for(j=0;j<dimension;++j)
            d->power[i][j] = ordered_power_list->power[j];
        struct maxentmc_product_list_link_struct * tp = ordered_power_list->product_data;
        j = 0;
        while(tp){
            d->product->entry[i].p[j].i1 = tp->i1;
            d->product->entry[i].p[j].i2 = tp->i2;
            struct maxentmc_product_list_link_struct * tp2 = tp;
            tp = tp2->next;
            free(tp2);
            ++j;
        }
        free(ordered_power_list->power);
        struct maxentmc_product_ordered_list_struct * opl2 = ordered_power_list;
        ordered_power_list = opl2->next;
        free(opl2);
    }

    MAXENTMC_POWER_SET_PROPERTY(MAXENTMC_POWER_ORDERED,d);
    maxentmc_power_update_max_power(d);

    return d;

}

int maxentmc_power_print_product(struct maxentmc_power_struct const * const d, FILE * const out)
{

    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(out);

    if(d->product){
        size_t i;
        for(i=0;i<d->size;++i){
            struct maxentmc_product_link_struct entry = d->product->entry[i];
            fprintf(out,"p[%zu] = p1[%zu]*p2[%zu]",i,entry.p[0].i1,entry.p[0].i2);
            size_t j, np = entry.np;
            for(j=1;j<np;++j)
                fprintf(out," + p1[%zu]*p2[%zu]",entry.p[j].i1,entry.p[j].i2);
            fputs("\n",out);
        }
    }
    else
        fputs("No product attached\n",out);

    return 0;

}


int maxentmc_power_multiply_add(size_t const size,
                                struct maxentmc_power_struct const * const p1, maxentmc_index_t const x1_power_dim, maxentmc_float_t const * const * const x1,
                                struct maxentmc_power_struct const * const p2, maxentmc_index_t const x2_power_dim, maxentmc_float_t const * const * const x2,
                                struct maxentmc_power_struct const * const p, maxentmc_index_t const x_power_dim, maxentmc_float_t * const * const x)
{

    if(size == 0){
        MAXENTMC_MESSAGE(stderr,"error: size is zero");
        return -1;
    }

    MAXENTMC_CHECK_NULL(p1);
    MAXENTMC_CHECK_NULL(x1);
    MAXENTMC_CHECK_NULL(p2);
    MAXENTMC_CHECK_NULL(x2);
    MAXENTMC_CHECK_NULL(p);
    MAXENTMC_CHECK_NULL(x);

    if(p->product == NULL){
        MAXENTMC_MESSAGE(stderr,"error: no product structure present");
        return -1;
    }

    maxentmc_float_t const * const * x1t;
    maxentmc_float_t const * const * x2t;

    if((p->product->p1 == p1) && (p->product->p2 == p2)){
        x1t = x1;
        x2t = x2;
    }
    else{
        if((p->product->p1 == p2) && (p->product->p2 == p1)){
            x1t = x2;
            x2t = x1;
        }
        else{
            MAXENTMC_MESSAGE(stderr,"error: attached powers do not match");
            return -1;
        }
    }

    size_t i;

    if(x_power_dim)
        if(x1_power_dim)
            if(x2_power_dim)
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[i][k] += x1t[i1][k] * x2t[i2][k];
                    }
                }
            else
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[i][k] += x1t[i1][k] * x2t[k][i2];
                    }
                }
        else
            if(x2_power_dim)
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[i][k] += x1t[k][i1] * x2t[i2][k];
                    }
                }
            else
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[i][k] += x1t[k][i1] * x2t[k][i2];
                    }
                }
    else
        if(x1_power_dim)
            if(x2_power_dim)
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[k][i] += x1t[i1][k] * x2t[i2][k];
                    }
                }
            else
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[k][i] += x1t[i1][k] * x2t[k][i2];
                    }
                }
        else
            if(x2_power_dim)
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[k][i] += x1t[k][i1] * x2t[i2][k];
                    }
                }
            else
                for(i=0;i<p->size;++i){
                    size_t j;
                    struct maxentmc_product_link_struct const entry = p->product->entry[i];
                    for(j=0;j<entry.np;++j){
                        size_t const i1 = entry.p[j].i1;
                        size_t const i2 = entry.p[j].i2;
                        size_t k;
                        for(k=0;k<size;++k)
                            x[k][i] += x1t[k][i1] * x2t[k][i2];
                    }
                }


    return 0;

}
