#include <string.h>
#include <netinet/in.h>
#include <netcdf.h>
#include "../core/maxentmc_list.h"
#include "test_read_old_data.h"

double poly_ntoh (double const x)
#if __BYTE_ORDER == __LITTLE_ENDIAN
{

    union {double x; char c[sizeof(double)];} in, out;

    in.x = x;

    size_t j;
    for (j = 0; j < sizeof(double)/2; ++j){
        out.c[j] = in.c[sizeof(double)-1-j];
        out.c[sizeof(double)-1-j] = in.c[j];
    }

    return out.x;

}
#else
{
    return x;
}
#endif


maxentmc_power_vector_t read_old_data(int const file_id, int const var_id, size_t const * const location)
{

    int ndims;
    if(nc_inq_varndims(file_id,var_id,&ndims) != NC_NOERR){
        MAXENTMC_MESSAGE(stderr,"error: nc_inq_varndims failed");
        return NULL;
    }

    if(ndims==0){
        MAXENTMC_MESSAGE(stderr,"error: insufficient number of dimensions")
        return NULL;
    }

    if((ndims>1) && (location == NULL)){
        MAXENTMC_MESSAGE(stderr,"error: location is needed, but NULL provided");
        return NULL;
    }

    int dim_ids[ndims];
    size_t dim_len[ndims], i;
    nc_inq_vardimid(file_id,var_id,dim_ids);
    for(i=0;i<ndims;++i)
        nc_inq_dimlen(file_id,dim_ids[i],dim_len+i);

    for(i=0;i<(ndims-1);++i)
        if(location[i]>=dim_len[i]){
            MAXENTMC_MESSAGE_VARARG(stderr,"error: location[%zu]=%zu exceeds dimension length %zu",i,location[i],dim_len[i]);
            return NULL;
        }

    size_t const raw_data_size = dim_len[ndims-1];
    char raw_data[raw_data_size+sizeof(int)] __attribute__ ((aligned(16)));
    size_t nc_start[ndims], nc_count[ndims];
    for(i=0;i<(ndims-1);++i){
        nc_start[i] = location[i];
        nc_count[i] = 1;
    }
    nc_start[ndims-1] = 0;
    nc_count[ndims-1] = raw_data_size;

    nc_get_vara_text(file_id,var_id,nc_start,nc_count,raw_data+sizeof(int));

    maxentmc_index_t dim, pow;
    /*size_t size;*/

    int * pt_i = (int *)raw_data;

    dim = ntohl(*++pt_i);
    pow = ntohl(*++pt_i);
    /*size = ntohl(*++pt_i);*/
    ++pt_i;

    maxentmc_list_t list = maxentmc_list_alloc(dim,1,MAXENTMC_LIST_ORDERED,MAXENTMC_LIST_ASCEND);

    /** Cycle through the powers now **/

    maxentmc_index_t power[dim], flag = 1;
    memset(power,0,dim*sizeof(maxentmc_index_t));

    size_t c=0, cur_pow = 0;
    double * pt_d = (double *)(++pt_i);

    do{

        maxentmc_float_t x = poly_ntoh(*pt_d++);
        maxentmc_list_insert_ca(list,power,&x);

        ++c;

        if(power[dim-1]<pow){
            ++(power[0]); ++cur_pow;

            for(i=0;cur_pow>pow;++i){
                cur_pow -= power[i]-1;
                power[i]=0;
                ++(power[i+1]);
            }
        }
        else
            flag = 0;

    }while(flag);

    maxentmc_power_vector_t constraints;

    maxentmc_list_create_power_vectors(list,&constraints);

    return constraints;

}
