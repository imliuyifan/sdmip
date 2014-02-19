/*
 ** This file is a mere shadow of the actual supomega.c file.  
 ** It assumes that, in each call, there is only one member
 ** of the given cell, and that that member encompasses all
 ** of the probability for every random variable in omega.
 */

#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "sdglobal.h"
#include "supomega.h"
#include "utility.h"
#include "rvgen.h"
#include "log.h"

/****************************************************************\
**  The function get_omega_idx() receives as parameters an
 **  array to be loaded with the indices corresponding to the 
 **  values of the stochastic elements to be used in the next
 **  run of the sub_problem, an array of members of a marriage
 **  for determining the areas over each stochastic element to
 **  be sampled and a value indicating the number of members 
 **  listed in the array.
 **  get_omega_idx() then determines which values will be returned
 **  and loads the array "observ" with the corresponding indices.
 **
 **  The function also return the sum of the indices as its 
 **  return value.
 **
 ** Currently, the members parameters are ignored, for 1 cell 
 ** implementation.
 \****************************************************************/
/* modified by Yifan 2012.07.02 */
sd_small get_omega_idx(sdglobal_type* sd_global, sd_small *observ, sd_small *members,
		sd_small num_members, sd_long *RUN_SEED)
{
	sd_small i, j;
	sd_small sum = 0;
	double val;
	/*  locate values of stochastic elements in distribution array */
	for (i = 0; i < sd_global->omegas.num_omega; i++)
	{
		val = scalit(0, 1, RUN_SEED);
		for (j = 0; val > sd_global->omegas.omega_probs[i][j]; j++)
			/* loop until value falls below the cdf at j */;
		sd_global->omegas.indices[i] = j;
		sum += j;
	}

	/* Reduce the full array of indices down to an encrypted array */
	encode(sd_global->omegas.indices, sd_global->omegas.key, observ,
			sd_global->omegas.num_omega);

#ifdef ENCODE
	{	int cnt;

		printf("\nOmega indices: ");
		for(cnt = 0; cnt < 12 && cnt < sd_global->omegas.num_omega; cnt++)
		printf("%d ", sd_global->omegas.indices[cnt]);
		printf("\nOmega values: ");
		for(cnt = 0; cnt < 12 && cnt < sd_global->omegas.num_omega; cnt++)
		printf("%lf ", sd_global->omegas.omega_vals[cnt][sd_global->omegas.indices[cnt]]);
		printf("\nOmega encoded: ");
		for(cnt = 0; cnt < 12 && cnt < sd_global->omegas.num_cipher; cnt++)
		printf("%d ", observ[cnt]);
		printf("\n");
	}
#endif

	return sum; /* This is the 1-norm of the uncoded indices */
}

/****************************************************************\
**  The function get_omega_vals() receives an array of indices
 **  to values in omegas and loads the array RT with the 
 **  corresponding values of omega.
 \****************************************************************/
double get_omega_vals(sdglobal_type* sd_global, sd_small *observ, double *RT)
{
	int i;

	decode(observ, sd_global->omegas.key, sd_global->omegas.indices, sd_global->omegas.num_omega);
	for (i = 0; i < sd_global->omegas.num_omega; i++)
		RT[i] = sd_global->omegas.omega_vals[i][sd_global->omegas.indices[i]];

#ifdef ENCODE
	{	int cnt;

		printf("\nRT Omega Vals: ");
		for(cnt = 0; cnt < 12 && cnt < sd_global->omegas.num_omega; cnt++)
		printf("%lf ", RT[cnt]);
		printf("\nOmega observ: ");
		for(cnt = 0; cnt < 12 && cnt < sd_global->omegas.num_cipher; cnt++)
		printf("%d ", observ[cnt]);
		printf("\n");
	}
#endif

	return 0.0; /* In place of 1-norm */
}

/* Yifan 2012.05.21 */

/****************************************************************\
 **  The function get_omega_from_file() get omegas from an external
 **  file. The position of the omega is stored in the _omegas.fidx_
 \****************************************************************/

double get_omega_vals_from_file(sdglobal_type *sd_global, int obs_idx, double *RT, sd_long *fidx)
{
  int i,status;
  
  fseek(sd_global->fptrOMEGA, fidx[obs_idx] , SEEK_SET);
  
  for(i=0; i < sd_global->omegas.num_omega; i++){
    status = fscanf(sd_global->fptrOMEGA, "%lf",&RT[i]);
    if (status<=0) {
      printf("Data error in Omega file.");
      exit(1);
    }
  }
  
  /* Deviation from mean will be used to calculate delta */
  for(i=0; i < sd_global->omegas.num_omega; i++){
    RT[i] -= sd_global->omegas.mean[i];
  }
  
  /* Get ready for reading the next omega */
  fidx[obs_idx+1] = ftell(sd_global->fptrOMEGA);
  
  return 0.0; /* In place of 1-norm */
}

/* Yifan 2012.05.21 */

/****************************************************************\
**  The function get_omega_col() loads a list of col numbers
 **  corresponding to the stochastic elements contained contained
 **  in the struct omegas.
 \****************************************************************/

sd_small get_omega_col(sdglobal_type* sd_global, sd_small *cols)
{
	int i;

	for (i = 0; i < sd_global->omegas.num_omega; i++)
		cols[i] = sd_global->omegas.col[i];

	return 0; /* In place of 1-norm */
}

/****************************************************************\
**  The function get_omega_col() loads a list of row numbers
 **  corresponding to the stochastic elements contained contained
 **  in the struct omegas.
 \****************************************************************/

sd_small get_omega_row(sdglobal_type* sd_global, sd_small *rows)
{
	int i;

	for (i = 0; i < sd_global->omegas.num_omega; i++)
		rows[i] = sd_global->omegas.row[i];

	return 0; /* In place of 1-norm */
}

/*
 ** This function will sort the data inside the _omegas_ structure
 ** according to the _col_ array.  In this way, the elements which
 ** correspond to the right hand side, R(omega), will appear at the
 ** beginning of all arrays, while elements of T(omega) will appear
 ** afterwards.
 */
void sort_omegas(sdglobal_type* sd_global)
{
	int i, j; /* loop counters */
	int min; /* smallest value in rest of array */
	int min_idx; /* location of minimum of the array */

	int itemp; /* used to swap integer elements */
	double *dptemp; /* used to swap double pointer elements */

	for (i = 0; i < sd_global->omegas.num_omega - 1; i++)
	{
		min = sd_global->omegas.col[i];
		min_idx = i;

		for (j = i + 1; j < sd_global->omegas.num_omega; j++)
			if (sd_global->omegas.col[j] < min)
			{
				min = sd_global->omegas.col[j];
				min_idx = j;
			}

		if (min_idx != i)
		{
			swap(sd_global->omegas.col[i], sd_global->omegas.col[min_idx],
					itemp);
			swap(sd_global->omegas.row[i], sd_global->omegas.row[min_idx],
					itemp);
			swap(sd_global->omegas.num_vals[i],
					sd_global->omegas.num_vals[min_idx], itemp);
			swap(sd_global->omegas.omega_vals[i],
					sd_global->omegas.omega_vals[min_idx], dptemp);
			swap(sd_global->omegas.omega_probs[i],
					sd_global->omegas.omega_probs[min_idx], dptemp);
		}
	}

}

/*
 ** This function will prepare the global _omegas_ structure for 
 ** encoding / decoding observations.  It should be called as soon as
 ** the rest of the structure is initialized.  It will create a key
 ** to be used in the coding process, and will intialize the _num_ciphers_
 ** field as the number of ints in the coded observation array.
 */
void cipher_omegas(sdglobal_type* sd_global)
{
	/* 
	 ** Need one key for every random variable, so you know where and
	 ** how each one has been encoded in the cipher array.
	 */
	sd_global->omegas.key = arr_alloc(sd_global->omegas.num_omega, one_key);

	/* Temporary array to use when encoding / decoding omega indices */
	sd_global->omegas.indices = arr_alloc(sd_global->omegas.num_omega, sd_small);

	/* Now form a key based on the range of each random variable in omegas */
	sd_global->omegas.num_cipher = form_key(sd_global->omegas.key,
			sd_global->omegas.num_vals, sd_global->omegas.num_omega);
}
