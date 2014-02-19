/***********************************************************************\
**
 ** delta.c
 **
 **
 **
 ** It is assumed that the zeroth element of every vector contains
 ** the 1-norm of that vector.  If the 1-norm is never necessary,
 ** the zeroth element may be un-initialized; however, the vector
 ** must start at 1. 
 **
 ** calc_delta_row()
 ** calc_delta_col()
 ** print_delta()
 ** new_delta()
 ** free_delta()
 **
 **
 ** History:
 **   04 Nov 91 - <Jason Mai> - created.
 **   03 Dec 91 - <Jason Mai> - added new_delta and printing functions
 **   14 Dec 91 - <Jason Mai> - prepared for testing / compilation. 
 **   18 Feb 92 - <Jason Mai> - updated with new structures (esp. num_type) 
 **   Should not be used without the consent of either Suvrajeet Sen
 **   or Jason Mai
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "utility.h"
#include "delta.h"
#include "omega.h"
#include "log.h"
#include "sdglobal.h"

/***********************************************************************\
** This function calculates a new row in the delta structure, based
 ** on a new dual vector, lambda_pi, by calculating lambda_pi X R
 ** and lambda_pi X T for all previous realizations of R(omega) and 
 ** T(omega).  It is assumed that the lambda vector is distinct
 ** from all previous ones, and thus a new row is warranted.
 \***********************************************************************/
void calc_delta_row(sdglobal_type* sd_global, delta_type *delta,
		lambda_type *lambda, omega_type *omega, num_type *num, int pi_idx)
{
	int obs;
	sparse_vect Romega;
	sparse_matrix Tomega;
	vector lamb_pi;
	vector pi_cross_T;

#ifdef TRACE
	printf("Inside calc_delta_row\n");
#endif

	if (sd_global->MALLOC)
	{
		printf("before init_R_T_omega\n");
		malloc_verify();
	}

	/* Initialize all vectors for calculations */
	init_R_T_omega(&Romega, &Tomega, omega, num);

	if (!(delta->val[pi_idx] = arr_alloc(num->iter, pi_R_T_type)))
		err_msg("Allocation", "calc_delta_row", "delta->val");

	lamb_pi = expand_vect(lambda->val[pi_idx], lambda->row, num->rv_rows,
			num->sub_rows);

	/* For all observations, calculate pi X R and pi X T */
	for (obs = 0; obs < omega->most; obs++)
		if (valid_omega_idx(omega, obs))
		{
			get_R_T_omega(sd_global, omega, obs);

			/* Multiply the new dual vector by previous observations of omega */
			/* Reduce the vector resulting from Pi x T to its sparse form */
			delta->val[pi_idx][obs].R = PIxR(lamb_pi, &Romega);
			pi_cross_T = PIxT(lamb_pi, &Tomega, num->mast_cols);
			delta->val[pi_idx][obs].T = reduce_vect(pi_cross_T, delta->col,
					num->rv_cols);
			mem_free(pi_cross_T);
		}

	mem_free(lamb_pi);
}

/***********************************************************************\
** This function calculates a new column in the delta structure,
 ** based on a new observation of omega.  Thus, lambda_pi X R and 
 ** lambda_pi X T are calculated for all values of lambda_pi, for
 ** the new R(omega) and T(omega).  Room in the array has already
 ** been allocated, so the function only fills it, in the column 
 ** specified by _obs_.  It is assumed that this observation is distinct 
 ** from all previous ones, and thus a new column must be calculated.  
 \***********************************************************************/
void calc_delta_col(sdglobal_type* sd_global, delta_type *delta,
		lambda_type *lambda, omega_type *omega, num_type *num, int obs)
{
	int pi_idx;
	sparse_vect Romega;
	sparse_matrix Tomega;
	vector lamb_pi;
	vector pi_cross_T;

#ifdef TRACE
	printf("Inside calc_delta_col\n");
#endif

	/* Initialize vectors for calculations */
	init_R_T_omega(&Romega, &Tomega, omega, num);
	get_R_T_omega(sd_global, omega, obs);

	/* For all dual vectors, lambda(pi), calculate pi X Romega and pi X Tomega */
	for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++)
	{
		/* Retrieve a new (sparse) dual vector, and expand it into a full vector */
		lamb_pi = expand_vect(lambda->val[pi_idx], lambda->row, num->rv_rows,
				num->sub_rows);

		/* Multiply the dual vector by the observation of Romega and Tomega */
		/* Reduce PIxT from its full vector form into a sparse vector */
		delta->val[pi_idx][obs].R = PIxR(lamb_pi, &Romega);
		pi_cross_T = PIxT(lamb_pi, &Tomega, num->mast_cols);
		delta->val[pi_idx][obs].T = reduce_vect(pi_cross_T, delta->col,
				num->rv_cols);
		mem_free(pi_cross_T);
		mem_free(lamb_pi);
		/*
		 ** Sloppy and slow to alloc & free pi_cross_T and lamb_pi every iteration!
		 */

	}
}

/***********************************************************************\
** This function frees a row of the delta structure, and all the 
 ** dynamically allocated memory associated with it.  Once the row has
 ** been deleted, the last row in the delta matrix is copied into its
 ** place.
 \***********************************************************************/
void drop_delta_row(delta_type *delta, lambda_type *lambda, omega_type *omega,
		int row)
{
	int col;

#ifdef TRACE
	printf("Inside drop_delta_row()\n");
#endif

	/* Free each calculation of Pi x T */
	for (col = 0; col < omega->most; col++)
		if (valid_omega_idx(omega, col))
			mem_free(delta->val[row][col].T);

	/* Free the specified row */
	mem_free(delta->val[row]);

	/* Copy the last row into the position of the vacated row */
	delta->val[row] = delta->val[lambda->cnt];
}

/***********************************************************************\
** This function removes a column from the delta structure.  Columns
 ** are harder to get rid of than rows -- they can't be deallocated.
 ** So, this function frees each Pi x T stored in the column, and fills
 ** it with NULL pointers, but it does NOT free the column itself.
 \***********************************************************************/
void drop_delta_col(delta_type *delta, lambda_type *lambda, int col)
{
	int row;

#ifdef TRACE
	printf("Inside drop_delta_col()\n");
#endif

	for (row = 0; row < lambda->cnt; row++)
	{
		mem_free(delta->val[row][col].T);
		delta->val[row][col].T = NULL;
		delta->val[row][col].R = 0.0;
	}
}

/***********************************************************************\
**
 \***********************************************************************/

/***********************************************************************\
**
 \***********************************************************************/

/***********************************************************************\
** This function is intended for debugging purposes.  It prints
 ** one delta_R and delta_T pair along with the columns of delta_T.
 \***********************************************************************/
void print_delta(delta_type *delta, num_type *num, int idx, int obs)
{
	int cnt;

	printf("\nDelta (%d,%d) :: R: %f \nDelta T: ", idx, obs,
			delta->val[idx][obs].R);
	for (cnt = 0; cnt <= num->rv_cols; cnt++)
		printf("%f ", delta->val[idx][obs].T[cnt]);
	printf("\nDelta cols:");
	for (cnt = 0; cnt <= num->rv_cols; cnt++)
		printf("%d ", delta->col[cnt]);
	printf("\n");
}

/***********************************************************************\
** This function creates a new delta structure with arrays of the specified
 ** size and returns a pointer to it.  Note that the pi X T vectors
 ** themselves are not allocated, since they will not all be filled with
 ** values.  (they are only filled as they are produced).
 ** Not even the arrays of pi_R_T_types are allocated, as this also
 ** occurs in calc_delta_row().  However, the column coordinates of the
 ** eventual multiplications are initialized, since they are known.
 \***********************************************************************/
delta_type *new_delta(int num_iter, coord_type *coord)
{
	delta_type *d;

#ifdef TRACE
	printf("Inside new_delta\n");
#endif

	if (!(d = (delta_type *) mem_malloc (sizeof(delta_type))))
		err_msg("Allocation", "new_delta", "d");

	if (!(d->val = arr_alloc(num_iter, pi_R_T_type*)))
		err_msg("Allocation", "new_delta", "d->val");

	d->col = coord->delta_col;

	return d;
}

/***********************************************************************\
** This function frees all the data associated with the delta
 ** three-dimensional matrix.  Since the size of the matrix is
 ** determined by the number of distinct realizations and dual
 ** vectors, these counts must be passed in.  Note that delta->col
 ** is not freed, since it belongs to the coord structure, and may
 ** be used by the next cell / soln.
 \***********************************************************************/
void free_delta(delta_type *delta, omega_type *omega, int num_lamb)
{
	int cnt, idx;

#ifdef TRACE
	printf("Inside free_delta\n");
#endif

	for (cnt = 0; cnt < num_lamb; cnt++)
	{
		for (idx = 0; idx < omega->most; idx++)
			if (valid_omega_idx(omega, idx))
				mem_free(delta->val[cnt][idx].T);
		mem_free(delta->val[cnt]);
	}mem_free(delta->val);
	mem_free(delta);
}
