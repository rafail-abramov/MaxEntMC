#include "test_convert_constraints.h"

int test_convert_constraints(char const * const filename, char const * const varname, char const * const my_varname)
{
    int file_id, var_id;

    if(nc_open(filename,NC_NOWRITE,&file_id) != NC_NOERR){
        fputs("Wrong filename\n",stderr);
        return -1;
    }

    if(nc_inq_varid(file_id,varname,&var_id) != NC_NOERR){
        fputs("Wrong varname\n",stderr);
        return -1;
    }

    int ndims=0;
    nc_inq_varndims(file_id,var_id,&ndims);

    if(ndims==0 || ndims>2){
        fputs("Wrong number of dimensions\n",stderr);
        return -1;
    }

    char my_filename[1024];
    FILE * out;
    int dimids[ndims];
    size_t i, dimlen;
    maxentmc_power_vector_t v;
    maxentmc_index_t dim, pow;
    switch(ndims){
        case 1:
            v = read_old_data(file_id,var_id,NULL);
            dim = maxentmc_power_vector_get_dimension(v);
            maxentmc_power_vector_get_max_power(v,&pow);
            if(my_varname)
                sprintf(my_filename,"constraints_%s_dim%u_pow%u.dat",my_varname,dim,pow);
            else
                sprintf(my_filename,"constraints_dim%u_pow%u.dat",dim,pow);
            out = fopen(my_filename,"wb");
            if(out){
                maxentmc_power_vector_fwrite_power(v,out);
                maxentmc_power_vector_fwrite_values(v,out);
            }
            fclose(out);
            maxentmc_power_vector_free(v);
            break;
        case 2:
            nc_inq_vardimid(file_id,var_id,dimids);
            nc_inq_dimlen(file_id,dimids[0],&dimlen);
            for(i=0;i<dimlen;++i){
                v = read_old_data(file_id,var_id,&i);
                dim = maxentmc_power_vector_get_dimension(v);
                maxentmc_power_vector_get_max_power(v,&pow);
                if(my_varname)
                    sprintf(my_filename,"./constraints_%s_dim%u_pow%u_%zu.dat",my_varname,dim,pow,i+1);
                else
                    sprintf(my_filename,"./constraints_dim%u_pow%u_%zu.dat",dim,pow,i+1);
                out = fopen(my_filename,"wb");
                if(out){
                    maxentmc_power_vector_fwrite_power(v,out);
                    maxentmc_power_vector_fwrite_values(v,out);
                }
                fclose(out);
                maxentmc_power_vector_free(v);
            }
            break;
    }

    return 0;

}

