/**************************************************************************\
** 
 ** supomega.h
 **
 **  This file contains prototypes for functions which interface
 **  between the SD solution code and the random variable 
 **  realization generation code.
 **  
 ** History:
 **   ?? ??? ?? - <Tim Burton> - created.
 **
 \**************************************************************************/

/**************************************************************************\
**  The function get_omega_idx() receives as parameters an
 ** array to be loaded with the indices corresponding to the 
 ** values of the stochastic elements to be used in the next
 ** run of the sub_problem, an array of members of a marriage
 ** for determining the areas over each stochastic element to
 ** be sampled and a value indicating the number of members 
 ** listed in the array.
 ** get_omega_idx() then determines which values will be returned
 ** and loads the array "observ" with the corresponding indices.
 **
 **  The function also return the sum of the indices as its 
 ** return value.
 \**************************************************************************/
/* modified by Yifan 2012.07.02 */

#ifndef SUPOMEGA_H_
#define SUPOMEGA_H_
#include "sdglobal.h"

sd_small get_omega_idx(sdglobal_type* sd_global, sd_small *observ, sd_small *members,
		sd_small num_members, sd_long *RUN_SEED);

/**************************************************************************\
**  The function get_omega_vals() receives an array of indices
 ** to values in omegas and loads the array RT with the
 ** corresponding values of omega.
 \**************************************************************************/

double get_omega_vals(sdglobal_type* sd_global, sd_small *observ, double *RT);
double get_omega_vals_from_file(sdglobal_type *sd_global, int obs_idx, double *RT, sd_long *fidx);

/**************************************************************************\
**  The function get_omega_col() loads a list of col numbers
 ** corresponding to the stochastic elements contained contained
 ** in the struct omegas.
 \**************************************************************************/

sd_small get_omega_col(sdglobal_type* sd_global, sd_small *cols);

/**************************************************************************\
**  The function get_omega_col() loads a list of col numbers
 ** corresponding to the stochastic elements contained contained
 ** in the struct omegas.
 \**************************************************************************/

sd_small get_omega_row(sdglobal_type* sd_global, sd_small *rows);

void sort_omegas(sdglobal_type* sd_global);
void cipher_omegas(sdglobal_type* sd_global);

#endif
