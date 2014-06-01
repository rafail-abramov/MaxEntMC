#ifndef MAXENTMC_READ_OLD_DATA_H_INCLUDED
#define MAXENTMC_READ_OLD_DATA_H_INCLUDED

#include "../core/maxentmc_vector.h"

struct maxentmc_power_vector_struct * read_old_data(int const file_id, int const var_id, size_t const * const location);

#endif // MAXENTMC_READ_OLD_DATA_H_INCLUDED
