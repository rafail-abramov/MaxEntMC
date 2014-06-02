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
#include "maxentmc_list.h"

struct maxentmc_list_link_struct{

    maxentmc_index_t * powers;
    maxentmc_float_t * data;
    struct maxentmc_list_link_struct * next;

};

/** Basic functions **/

struct maxentmc_list_struct * maxentmc_list_alloc(maxentmc_index_t const dim, size_t const data_size,
                                  enum MAXENTMC_LIST_POWER_ORDER order,
                                  enum MAXENTMC_LIST_POWER_INSERTION_ORDER insertion_order)
{
    if(dim == 0){
        MAXENTMC_MESSAGE(stderr,"error: zero dimension");
        return NULL;
    }
    if(data_size == 0){
        MAXENTMC_MESSAGE(stderr,"error: zero data_size");
        return NULL;
    }
    struct maxentmc_list_struct * d = malloc(sizeof(struct maxentmc_list_struct));
    if(d==NULL){
        MAXENTMC_MESSAGE(stderr,"error: malloc returned NULL");
        return NULL;
    }
    d->dimension = dim;
    d->order = order;
    d->insertion_order = insertion_order;
    d->data_size = data_size;
    d->list_size = 0;
    d->link = NULL;
    return d;
}

int maxentmc_list_clear(struct maxentmc_list_struct * const d)
{
    MAXENTMC_CHECK_NULL(d);

    struct maxentmc_list_link_struct * link = d->link;
    while(link){
        struct maxentmc_list_link_struct * link2 = link;
        link = link2->next;
        free(link2);
    }

    d->list_size = 0;
    d->link = NULL;

    return 0;
}

void maxentmc_list_free(struct maxentmc_list_struct * const d)
{
    if(d){
        maxentmc_list_clear(d);
        free(d);
   }
}

static int maxentmc_list_find_link(struct maxentmc_list_struct const * const d, maxentmc_index_t const * const powers, struct maxentmc_list_link_struct ** const link)
{
    maxentmc_index_t const dim = d->dimension;
    *link = d->link;

    if(*link==NULL)
        return -1; /** The list is empty, attach to the beginning of the list **/

    if(d->order == MAXENTMC_LIST_ORDERED){

        int i = maxentmc_power_comparison(dim,(*link)->powers,powers);

        if(((d->insertion_order == MAXENTMC_LIST_DESCEND) && (i>0)) || ((d->insertion_order == MAXENTMC_LIST_ASCEND) && (i<0)))
            return -1; /** The list is nonempty, but still attach to the beginning of the list **/

        if(i==0)
            return 2; /** Found it, this is the first link, *link points to it **/

        /** If we are here, need to iterate over the list **/

        if(d->insertion_order == MAXENTMC_LIST_DESCEND)
            while((*link)->next && ((i=maxentmc_power_comparison(dim,(*link)->next->powers,powers))<0))
                *link = (*link)->next;
        else
            while((*link)->next && ((i=maxentmc_power_comparison(dim,(*link)->next->powers,powers))>0))
                *link = (*link)->next;

        /** Here we have three possibilities:
            1) The list ended ((*link)->next == NULL) -- then return 0;
            2) The power comparison returned > 0 -- also return 0;
            3) The power comparison returned == 0 -- return 1. **/

        if(i==0)
            return 1;
        else
            return 0;

    }
    else{

        /** The list is nonempty and not ordered, traverse it and see if we already have that power. **/

        /** Check the first link **/

        if(memcmp((*link)->powers,powers,sizeof(maxentmc_index_t)*dim) == 0)
            return 2; /** Found it, this is the first link, *link points to it **/

        /** If we are here, need to iterate over the list **/

        while((*link)->next && memcmp((*link)->next->powers,powers,sizeof(maxentmc_index_t)*dim))
                *link = (*link)->next;

        /** Here we have three possibilities:
            1) The list ended ((*link)->next == NULL) -- then return 0;
            2) The power comparison returned == 0 -- return 1. **/

        if((*link)->next)
            return 1;
        else
            return 0;

    }

}

#define MAXENTMC_LIST_LINK_HEADER_SIZE MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_index_t),sizeof(struct maxentmc_list_link_struct))

static struct maxentmc_list_link_struct * maxentmc_list_alloc_link(struct maxentmc_list_struct const * const d,
                                                                   maxentmc_index_t const * const powers,
                                                                   struct maxentmc_list_link_struct * const link_next)
{
    size_t const link_power_size = MAXENTMC_ALIGNED_SIZE(sizeof(maxentmc_float_t),MAXENTMC_LIST_LINK_HEADER_SIZE+sizeof(maxentmc_index_t)*d->dimension);

    struct maxentmc_list_link_struct * link = malloc(link_power_size+sizeof(maxentmc_float_t)*d->data_size);

    link->powers = MAXENTMC_INCREMENT_POINTER(link,MAXENTMC_LIST_LINK_HEADER_SIZE);

    if(d->data_size)
        link->data = MAXENTMC_INCREMENT_POINTER(link,link_power_size);
    else
        link->data = NULL;

    link->next = link_next;

    memcpy(link->powers,powers,sizeof(maxentmc_index_t)*d->dimension);

    return link;
}

int maxentmc_list_insert_ca(struct maxentmc_list_struct * const d, maxentmc_index_t const * const powers, maxentmc_float_t const * const data)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(powers);
    if(d->data_size){
        MAXENTMC_CHECK_NULL(data);
    }

    struct maxentmc_list_link_struct * link;

    switch(maxentmc_list_find_link(d, powers, &link)){
        case -1:
            d->link = maxentmc_list_alloc_link(d, powers, d->link);
            if(d->link==NULL){
                MAXENTMC_MESSAGE(stderr,"error: cannot allocate link")
                return -1;
            }
            link = d->link;
            ++(d->list_size);
            break;
        case 0:
            link->next = maxentmc_list_alloc_link(d, powers, link->next);
            if(link->next==NULL){
                MAXENTMC_MESSAGE(stderr,"error: cannot allocate link")
                return -1;
            }
            link = link->next;
            ++(d->list_size);
            break;
        case 1:
            link = link->next;
            break;
    }
    if(d->data_size)
        memcpy(link->data,data,sizeof(maxentmc_float_t)*d->data_size);

    return 0;
}

int maxentmc_list_insert(struct maxentmc_list_struct * const d, ...)
{
    MAXENTMC_CHECK_NULL(d);
    maxentmc_index_t powers[d->dimension];
    maxentmc_float_t data[d->data_size];
    size_t i;
    va_list ap;
    va_start(ap,d);
    for(i=0;i<d->dimension;++i)
        powers[i] = va_arg(ap,int);
    for(i=0;i<d->data_size;++i)
        data[i] = va_arg(ap,maxentmc_float_t);
    va_end(ap);
    return maxentmc_list_insert_ca(d,powers,data);
}

int maxentmc_list_delete_ca(struct maxentmc_list_struct * const d, maxentmc_index_t const * const powers)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(powers);
    struct maxentmc_list_link_struct * link, * link2;

    switch(maxentmc_list_find_link(d, powers, &link)){
        case 1:
            link2 = link->next->next;
            free(link->next);
            link->next = link2;
            --(d->list_size);
            break;
        case 2:
            d->link = link->next;
            free(link);
            --(d->list_size);
            break;
        default:
            return -1;
    }

    return 0;
}

int maxentmc_list_delete(struct maxentmc_list_struct * const d, ...)
{
    MAXENTMC_CHECK_NULL(d);
    maxentmc_index_t powers[d->dimension];
    maxentmc_index_t i;
    va_list ap;
    va_start(ap,d);
    for(i=0;i<d->dimension;++i)
        powers[i] = va_arg(ap,int);
    va_end(ap);
    return maxentmc_list_delete_ca(d,powers);
}

int maxentmc_list_print(struct maxentmc_list_struct const * const d, FILE * out)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(out);

    fprintf(out,"dimension = %u, data_size = %zu, list_size = %zu\n",d->dimension,d->data_size,d->list_size);
    struct maxentmc_list_link_struct * link = d->link;

    while(link){
        size_t i;
        for(i=0;i<d->dimension;++i)
            fprintf(out,"%u ",link->powers[i]);
        fputs("|",out);
        for(i=0;i<d->data_size;++i)
            fprintf(out," %g",link->data[i]);
        fputs("\n",out);
        link = link->next;
    }

    return 0;

}

/** Power, vector and matrix creation functions **/

struct maxentmc_power_struct * maxentmc_list_create_power(struct maxentmc_list_struct const * const d)
{
    MAXENTMC_CHECK_NULL_PT(d);

    if(d->list_size==0){
        MAXENTMC_MESSAGE(stderr,"error: list is empty");
        return NULL;
    }

    struct maxentmc_power_struct * p = maxentmc_power_alloc(d->dimension,d->list_size);
    if(p == NULL){
        MAXENTMC_MESSAGE(stderr,"error: could not allocate powers");
        return NULL;
    }

    struct maxentmc_list_link_struct * link = d->link;
    size_t i;
    if((d->order == MAXENTMC_LIST_FORWARD) || ((d->order == MAXENTMC_LIST_ORDERED) && (d->insertion_order == MAXENTMC_LIST_DESCEND)))
        for(i=0;i<d->list_size;++i){
            memcpy(p->power[i],link->powers,sizeof(maxentmc_index_t)*d->dimension);
            link=link->next;
        }
    else
        for(i=0;i<d->list_size;++i){
            memcpy(p->power[d->list_size-i-1],link->powers,sizeof(maxentmc_index_t)*d->dimension);
            link=link->next;
        }

    if(d->order == MAXENTMC_LIST_ORDERED)
        MAXENTMC_POWER_SET_PROPERTY(MAXENTMC_POWER_ORDERED,p);
    else
        MAXENTMC_POWER_CLEAR_PROPERTY(MAXENTMC_POWER_ORDERED,p);

    maxentmc_power_update_max_power(p);

    return p;
}

struct maxentmc_power_vector_struct ** maxentmc_list_create_power_vector_array(struct maxentmc_list_struct const * const d)
{
    MAXENTMC_CHECK_NULL_PT(d);

    if(d->data_size == 0){
        MAXENTMC_MESSAGE(stderr,"error: list has no data");
        return NULL;
    }

/*
    struct maxentmc_power_vector_struct ** array = malloc(sizeof(struct maxentmc_power_vector_struct *)*d->data_size);
    if(array == NULL){
        MAXENTMC_MESSAGE(stderr,"error: malloc returned NULL");
        return NULL;
    }
*/

    struct maxentmc_power_vector_struct * array[d->data_size];

    struct maxentmc_power_struct * p = maxentmc_list_create_power(d);
    if(p==NULL){
        MAXENTMC_MESSAGE(stderr,"error: could not create powers");
        free(array);
        return NULL;
    }

    size_t i = 0;
    do
        array[i++] = maxentmc_power_vector_alloc_from_power(p);
    while(array[i-1] && (i<d->data_size));

    if(array[--i]==NULL){
        MAXENTMC_MESSAGE(stderr,"error: could not allocate array");
        while(i)
            maxentmc_power_vector_free(array[--i]);
    /*
        free(array);
    */
        return NULL;
    }

    struct maxentmc_list_link_struct * link = d->link;
    if((d->order == MAXENTMC_LIST_FORWARD) || ((d->order == MAXENTMC_LIST_ORDERED) && (d->insertion_order == MAXENTMC_LIST_DESCEND)))
        for(i=0;i<d->list_size;++i){
            size_t j;
            for(j=0;j<d->data_size;++j)
                array[j]->gsl_vec.data[i] = link->data[j];
            link=link->next;
        }
    else
        for(i=0;i<d->list_size;++i){
            size_t j;
            for(j=0;j<d->data_size;++j)
                array[j]->gsl_vec.data[d->list_size-i-1] = link->data[j];
            link=link->next;
        }

/*
    return array;
*/

    struct maxentmc_power_vector_struct ** array2;
    int status;

    MAXENTMC_ALLOC(array2,sizeof(struct maxentmc_power_vector_struct *)*d->data_size,status)
    if(status){
        i = d->data_size;
        while(i)
            maxentmc_power_vector_free(array[--i]);
        return NULL;
    }
    memcpy(array2,array,sizeof(struct maxentmc_power_vector_struct *)*d->data_size);

    return array2;
}

int maxentmc_list_create_power_vectors(struct maxentmc_list_struct const * const d, ...)
{
    MAXENTMC_CHECK_NULL(d);
    if(d->data_size==0){
        MAXENTMC_MESSAGE(stderr,"error: list has no data");
        return -1;
    }
    struct maxentmc_power_vector_struct ** array = maxentmc_list_create_power_vector_array(d);
    if(array == NULL){
        MAXENTMC_MESSAGE(stderr,"error: vec_array could not be allocated");
        return -1;
    }
    va_list ap;
    va_start(ap,d);
    size_t i;
    for(i=0;i<d->data_size;++i)
        *va_arg(ap,struct maxentmc_power_vector_struct **)=array[i];
    va_end(ap);
    free(array);
    return 0;
}
