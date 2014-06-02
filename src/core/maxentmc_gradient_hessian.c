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

#include <errno.h>
#include <string.h>
#include "maxentmc_gradient_hessian.h"

struct maxentmc_LGH_power_list_struct {

    struct maxentmc_LGH_power_list_struct * next;
    struct maxentmc_power_struct * p;
    maxentmc_index_t have_L, have_G, have_H;
    size_t L_element, * G_elements, ** H_elements;

};

static maxentmc_index_t maxentmc_LGH_find_zero_power(struct maxentmc_power_struct const * const p, size_t * const element)
{

    if((p==NULL) || (element==NULL)){
        MAXENTMC_MESSAGE(stdout,"error: NULL pointer provided");
        return 0;
    }

    size_t i=0;
    size_t total_power;
    do{
        total_power = 0;
        maxentmc_index_t j;
        for(j=0;j<p->dimension;++j)
            total_power += p->power[i][j];
    }while(total_power && ((++i)<p->size));

    if(total_power)
        return 0;
    else{
        *element = i;
        return 1;
    }

}

/** Allocation/deallocation routines **/

struct maxentmc_LGH_struct * maxentmc_LGH_alloc(struct maxentmc_power_vector_struct const * const constraints)
{
    MAXENTMC_CHECK_NULL_PT(constraints);

    struct maxentmc_power_struct const * const powers = constraints->powers;

    struct maxentmc_LGH_struct * d;

    int status;

    MAXENTMC_ALLOC(d,MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,sizeof(struct maxentmc_LGH_struct)),status)
    if(status)
       return NULL;

    d->powers = (struct maxentmc_power_struct *)powers;
    d->L_element = 0;
    d->have_L = maxentmc_LGH_find_zero_power(powers,&(d->L_element));
    d->LGH_powers = NULL;

    return d;
}

void maxentmc_LGH_free(struct maxentmc_LGH_struct * const d)
{

    if(d){
        /** Destroy list here **/
        struct maxentmc_LGH_power_list_struct * LGH_temp = d->LGH_powers;
        while(LGH_temp){
            struct maxentmc_LGH_power_list_struct * LGH_temp_2 = LGH_temp;
            LGH_temp = LGH_temp_2->next;
            free(LGH_temp_2);
        }
        /** End destroying list **/

        free(d);
    }

}


/** List addition routines **/

int maxentmc_LGH_add_power_vector(struct maxentmc_LGH_struct * const d, struct maxentmc_power_vector_struct const * const v)
{

    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(v);

    struct maxentmc_power_struct const * const p = v->powers;

    /** First, check if p is already here **/

    if(p == d->powers){
        MAXENTMC_MESSAGE(stderr,"error: power is already present");
        return -1;
    }

    struct maxentmc_LGH_power_list_struct * LGH_temp = d->LGH_powers;

    if(LGH_temp){
        if(LGH_temp->p == p){
            MAXENTMC_MESSAGE(stderr,"error: power is already present");
            return -1;
        }
        while(LGH_temp->next){
            LGH_temp = LGH_temp->next;
            if(LGH_temp->p == p){
                MAXENTMC_MESSAGE(stderr,"error: power is already present");
                return -1;
            }
        }
    }

    /** If we are here, identical power was not found **/

    maxentmc_index_t const dim = d->powers->dimension;
    size_t const size = d->powers->size;

    /** Now, check the dimension **/

    if(dim != p->dimension){
        MAXENTMC_MESSAGE(stderr,"error: dimensions do not match");
        return -1;
    }

    /** Dimensions match **/

    /** Check if we can use it for Lagrangian computation **/

    size_t L_element = 0;
    maxentmc_index_t have_L = maxentmc_LGH_find_zero_power(p,&L_element);

    /** Check if we can use it for gradient computation **/
    /** For now, find matching powers by brute force **/

    maxentmc_index_t have_G = 0;
    size_t G_elements[size], i=0;
    do{
        size_t j=0;
        while(memcmp(d->powers->power[i],p->power[j],sizeof(maxentmc_index_t)*dim)
               && ((++j)<p->size));
        if(j<p->size)
            G_elements[i++] = j;
        else
            i = size+1;
    }while(i<size);

    if(i==size)
        have_G = 1;

    /** Check if we can use it for Hessian computation **/
    /** Currently, the check is used via product **/
    maxentmc_index_t have_H = 0;
    if(p->product && (p->product->p1 == d->powers) && (p->product->p2 == d->powers))
        have_H = 1;

    if(!have_L && !have_G && !have_H){
        MAXENTMC_MESSAGE(stderr,"error: provided power cannot be used for anything");
        return -1;
    }

    /** If we are here, there is some use to the power, so need to add it **/

    /** First, we need to figure out the allocation size **/

    size_t link_size = 0, G_increment = 0, H_increment = 0, H0_increment = 0;
    if(!have_G && !have_H){
        link_size = MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,sizeof(struct maxentmc_LGH_power_list_struct));
        MAXENTMC_MESSAGE(stdout,"DEBUG: !have_G && !have_H");
    }
    else{
        /** Either have_G is present, or have_H, or both **/
        if(!have_H){
            /** only have_G **/
            link_size = MAXENTMC_ALIGNED_SIZE(sizeof(size_t),sizeof(struct maxentmc_LGH_power_list_struct));
            G_increment = link_size;
            link_size = MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,link_size+sizeof(size_t)*size);
            /*MAXENTMC_MESSAGE(stdout,"DEBUG: !have_H");*/
        }
        if(!have_G){
            /** only have_H **/
            link_size = MAXENTMC_ALIGNED_SIZE(sizeof(size_t *),sizeof(struct maxentmc_LGH_power_list_struct));
            H_increment = link_size;
            link_size = MAXENTMC_ALIGNED_SIZE(sizeof(size_t),link_size+sizeof(size_t *)*size);
            H0_increment = link_size;
            link_size = MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,link_size+sizeof(size_t)*size*size);
            /*MAXENTMC_MESSAGE(stdout,"DEBUG: !have_G");*/
        }
        if(have_G && have_H){
            /** both have_G and have_H **/
            link_size = MAXENTMC_ALIGNED_SIZE(sizeof(size_t *),sizeof(struct maxentmc_LGH_power_list_struct));
            H_increment = link_size;
            link_size = MAXENTMC_ALIGNED_SIZE(sizeof(size_t),link_size+sizeof(size_t *)*size);
            G_increment = link_size;
            H0_increment = G_increment+sizeof(size_t)*size;
            link_size = MAXENTMC_ALIGNED_SIZE(MAXENTMC_CACHE_LINE_SIZE,link_size+sizeof(size_t)*size*(size+1));
            /*MAXENTMC_MESSAGE(stdout,"DEBUG: have_G && have_H");*/
        }
    }

    if(LGH_temp){
        int status;
        MAXENTMC_ALLOC(LGH_temp->next,link_size,status);
        if(status)
            return -1;
        LGH_temp = LGH_temp->next;
    }
    else{
        int status;
        MAXENTMC_ALLOC(d->LGH_powers,link_size,status);
        if(status)
            return -1;
        LGH_temp = d->LGH_powers;
    }
    LGH_temp->next = NULL;
    LGH_temp->p = (struct maxentmc_power_struct *)p;
    LGH_temp->have_L = have_L;
    LGH_temp->have_G = have_G;
    LGH_temp->have_H = have_H;
    LGH_temp->L_element = L_element;
    if(have_G){
        LGH_temp->G_elements = MAXENTMC_INCREMENT_POINTER(LGH_temp,G_increment);
        memcpy(LGH_temp->G_elements,G_elements,sizeof(size_t)*size);
    }
    else
        LGH_temp->G_elements = NULL;
    if(have_H){
        LGH_temp->H_elements = MAXENTMC_INCREMENT_POINTER(LGH_temp,H_increment);
        LGH_temp->H_elements[0] = MAXENTMC_INCREMENT_POINTER(LGH_temp,H0_increment);
        for(i=1;i<size;++i)
            LGH_temp->H_elements[i] = LGH_temp->H_elements[i-1] + size;
        for(i=0;i<p->size;++i){
            size_t j;
            struct maxentmc_product_link_struct const entry = p->product->entry[i];
            for(j=0;j<entry.np;++j){
                size_t const i1 = entry.p[j].i1;
                size_t const i2 = entry.p[j].i2;
                LGH_temp->H_elements[i1][i2] = i;
            }
        }
    }
    else
        LGH_temp->H_elements = NULL;

    return 0;
}


/** Computation routines **/

int maxentmc_LGH_compute_lagrangian(struct maxentmc_LGH_struct const * const d, struct maxentmc_power_vector_struct const * const moments,
                                  struct maxentmc_power_vector_struct const * const constraints,
                                  struct maxentmc_power_vector_struct const * const multipliers, maxentmc_float_t * const L)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(moments);
    MAXENTMC_CHECK_NULL(constraints);
    MAXENTMC_CHECK_NULL(multipliers);
    MAXENTMC_CHECK_NULL(L);

    /** Check if powers match **/

    if((d->powers != constraints->powers) || (d->powers != multipliers->powers)){
        MAXENTMC_MESSAGE(stderr,"error: powers do not match");
        return -1;
    }

    /** Look for an appropriate power for moments **/

    maxentmc_index_t found = 0;

    /** First, check the constraint/multiplier power **/

    if((d->powers == moments->powers) && d->have_L){
        found = 1;
        *L = moments->gsl_vec.data[d->L_element*moments->gsl_vec.stride];
    }
    else{
        /** Iterate over the list **/

        struct maxentmc_LGH_power_list_struct * LGH_temp = d->LGH_powers;

        while(!found && LGH_temp){
            if((LGH_temp->p == moments->powers) && LGH_temp->have_L){
                found = 1;
                *L = moments->gsl_vec.data[LGH_temp->L_element*moments->gsl_vec.stride];
            }
            else
                LGH_temp = LGH_temp->next;
        }

    }

    if(!found){
        MAXENTMC_MESSAGE(stderr,"error: could not find appropriate power");
        return -1;
    }

    size_t i;

    for(i=0;i<d->powers->size;++i)
        *L -= multipliers->gsl_vec.data[i*multipliers->gsl_vec.stride]*constraints->gsl_vec.data[i*constraints->gsl_vec.stride];

    return 0;
}

int maxentmc_LGH_compute_gradient(struct maxentmc_LGH_struct const * const d, struct maxentmc_power_vector_struct const * const moments,
                                struct maxentmc_power_vector_struct const * const constraints, maxentmc_gsl_vector_t * const G)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(moments);
    MAXENTMC_CHECK_NULL(constraints);
    MAXENTMC_CHECK_NULL(G);

    /** Check if powers match **/

    if(d->powers != constraints->powers){
        MAXENTMC_MESSAGE(stderr,"error: powers do not match");
        return -1;
    }

    if(constraints->gsl_vec.size != G->size){
        MAXENTMC_MESSAGE(stderr,"error: sizes do not match");
        return -1;
    }

    size_t const G_stride = G->stride;

    /** Look for an appropriate power for moments **/

    maxentmc_index_t found = 0;

    /** First, check the constraint/multiplier power **/

    size_t i;

    if(d->powers == moments->powers){
        found = 1;
        for(i=0;i<d->powers->size;++i)
            G->data[G_stride*i] = moments->gsl_vec.data[i*moments->gsl_vec.stride];
    }
    else{
        /** Iterate over the list **/

        struct maxentmc_LGH_power_list_struct * LGH_temp = d->LGH_powers;

        while(!found && LGH_temp){
            if((LGH_temp->p == moments->powers) && LGH_temp->have_G){
                found = 1;
                for(i=0;i<d->powers->size;++i)
                    G->data[G_stride*i] = moments->gsl_vec.data[LGH_temp->G_elements[i]*moments->gsl_vec.stride];
            }
            else
                LGH_temp = LGH_temp->next;
        }
    }

    if(!found){
        MAXENTMC_MESSAGE(stderr,"error: could not find appropriate power");
        return -1;
    }

    for(i=0;i<d->powers->size;++i)
        G->data[G_stride*i] -= constraints->gsl_vec.data[i*constraints->gsl_vec.stride];

    return 0;

}

int maxentmc_LGH_compute_hessian(struct maxentmc_LGH_struct const * const d, struct maxentmc_power_vector_struct const * const moments,
                               maxentmc_gsl_matrix_t * const H)
{
    MAXENTMC_CHECK_NULL(d);
    MAXENTMC_CHECK_NULL(moments);
    MAXENTMC_CHECK_NULL(H);

    /** Check if H is square **/

    if(H->size1 != H->size2){
        MAXENTMC_MESSAGE(stderr,"error: Hessian is not a square matrix");
        return -1;
    }

    /** Check if sizes match **/

    if(d->powers->size != H->size1){
        MAXENTMC_MESSAGE(stderr,"error: sizes do not match");
        return -1;
    }

    size_t const H_tda = H->tda;

    /** Look for an appropriate power for moments **/

    maxentmc_index_t found = 0;

    /** Iterate over the list **/

    struct maxentmc_LGH_power_list_struct * LGH_temp = d->LGH_powers;

    while(!found && LGH_temp){
        if((LGH_temp->p == moments->powers) && LGH_temp->have_H){
            found = 1;
            size_t i;
            for(i=0;i<d->powers->size;++i){
                size_t j;
                for(j=0;j<d->powers->size;++j)
                    H->data[i*H_tda+j] = moments->gsl_vec.data[LGH_temp->H_elements[i][j]*moments->gsl_vec.stride];
            }
        }
        else
            LGH_temp = LGH_temp->next;
    }

    if(!found){
        MAXENTMC_MESSAGE(stderr,"error: could not find appropriate power");
        return -1;
    }

    return 0;
}
