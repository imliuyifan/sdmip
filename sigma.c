/***********************************************************************\
**
 ** sigma.c
 **
 **
 **
 ** It is assumed that the zeroth element of every vector contains
 ** the 1-norm of that vector.  If the 1-norm is never necessary,
 ** the zeroth element may be un-initialized; however, the vector
 ** must start at 1. 
 **
 ** calc_sigma()
 ** print_sigma()
 ** new_sigma()
 ** free_sigma()
 **
 ** History:
 **   ?? Nov 1991 - <Jason Mai> - created.
 **   17 Dec 1991 - <Jason Mai> - fixed compilation errors.
 **   19 Feb 1992 - <Jason Mai> - updated with new structures.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "utility.h"
#include "sigma.h"
#include "log.h"
#include "sdglobal.h"

/***********************************************************************\
** This function calculates Pi X R and Pi X T for a new Pi.  Then it
 ** compares these results to all previous calculations of Pi X R and 
 ** Pi X T.  If it finds a duplicate, it returns the index of the duplicate.
 ** Otherwise, it adds the two to the end of the sigma structure, and 
 ** returns the index of the last element.  Note that the random 
 ** components of R and T (which depend on omega) are not used.  It is
 ** assumed that the means of the random variables are stored at their  
 ** corresponding location in Rbar and Tbar.  Also note that each entry
 ** in sigma has only one entry in lambda associated with it, even though
 ** it *is* possible to have a duplicate sigma and distinct lambda.  So,
 ** if lambda is new, a new entry in sigma is automatically calculated.
 \***********************************************************************/
int calc_sigma(sdglobal_type* sd_global, cell_type *c, sigma_type *sigma,
		num_type *num, vector pi_k, sparse_vect *Rbar, sparse_matrix *Tbar,
		int lamb_idx, BOOL new_lamb, BOOL *new_sigma)
{
	double pi_R; /* scalar value of Pi x Rbar */
	double Mu_R; /* added by Yifan to store Mu x R */
	vector pi_T; /* reduced vector product of Pi x Tbar */
	vector temp; /* temporary full product of Pi x Tbar */
	int cnt; /* steps through sigma structure */

#ifdef TRACE
	printf("Inside calc_sigma\n");
#endif

	/* Add these new values of pi_R and pi_T and store index to lambda Yifan*/
	Mu_R = compute_Mu(c->subprob, num->sub_cols);
	//printf("Mu_R = %f \n", Mu_R);

	pi_R = PIxR(pi_k, Rbar) + Mu_R;

	temp = PIxT(pi_k, Tbar, num->mast_cols);
	pi_T = reduce_vect(temp, sigma->col, num->nz_cols);
	mem_free(temp);

	if (!new_lamb)
	{
		/* Compare pi_R and pi_T with all previous values */
		/* The TOLERANCE in the following lines are substituted by
		 TOLERANCE.   zl */

		for (cnt = 0; cnt < sigma->cnt; cnt++)
			/* Add <= and DBL_ABS in case pi_R is zero or negative 04/25/2013 Yifan */
			if (DBL_ABS(pi_R - sigma->val[cnt].R)
					<= sd_global->config.TOLERANCE * DBL_ABS(pi_R))
				if (equal_arr(pi_T, sigma->val[cnt].T, num->nz_cols,
						sd_global->config.TOLERANCE))
				{
					if (sigma->lamb[cnt] == lamb_idx)
					{
						mem_free(pi_T);
						*new_sigma = FALSE;
						return cnt;
					}

				}
	}

	if (sd_global->MALLOC)
	{
		printf("Before assignment of new sigma\n");
		printf("Sigma->cnt is: %d", sigma->cnt);
		/*malloc_verify();*/
	}

	*new_sigma = TRUE;
	sigma->val[sigma->cnt].R = pi_R;
	sigma->val[sigma->cnt].T = pi_T;
	sigma->lamb[sigma->cnt] = lamb_idx;
	sigma->ck[sigma->cnt] = c->k;

	if (sd_global->MALLOC)
	{
		printf("Before returning\n");
		printf("Sigma->cnt is: %d", sigma->cnt);
		/* malloc_verify(); */
		printf("After malloc_verify() before returning\n");
	}
	return sigma->cnt++;
}

/***********************************************************************\
** This function prints a given entry in the sigma matrix.
 ** It is meant for debugging purposes.
 \***********************************************************************/
void print_sigma(sigma_type *sigma, num_type *num, int idx)
{
	int cnt;

	printf("\nSigma %d:: lamb:%d R:%f T:", idx, sigma->lamb[idx],
			sigma->val[idx].R);
	for (cnt = 0; cnt <= num->nz_cols; cnt++)
		printf("%f ", sigma->val[idx].T[cnt]);
	printf("\nSigma cols:: ");
	for (cnt = 0; cnt <= num->nz_cols; cnt++)
		printf("%d ", sigma->col[cnt]);
	printf("\n");
}

/***********************************************************************\
** This function creates a new sigma structure, and allocates memory
 ** for the arrays associated with it.  It returns a pointer to this
 ** structure.  Some pi X T vectors are also allocated, according to
 ** the num_vals parameter  (num_vals is expected to be less than
 ** num_sigmas, so that there is room for further work).  Note that
 ** memory for sigma->col is not allocated, but is taken from prob.
 \***********************************************************************/
sigma_type *new_sigma(int num_iter, int num_nz_cols, int num_pi,
		coord_type *coord)
{
	sigma_type *sigma;
	int cnt;

#ifdef TRACE
	printf("Inside new_sigma\n");
#endif

	if (!(sigma = (sigma_type *) mem_malloc (sizeof(sigma_type))))
		err_msg("Allocation", "new_sigma", "sigma");

	if (!(sigma->lamb = arr_alloc(num_iter, int)))
		err_msg("Allocation", "new_sigma", "sigma->lamb");

	if (!(sigma->ck = arr_alloc(num_iter, int)))
		err_msg("Allocation", "new_sigma", "sigma->ck");

	if (!(sigma->val = arr_alloc(num_iter, pi_R_T_type)))
		err_msg("Allocation", "new_sigma", "sigma->val");

	for (cnt = 0; cnt < num_pi && cnt < num_iter; cnt++)
		if (!(sigma->val[cnt].T = arr_alloc(num_nz_cols+1, double)))
			err_msg("Allocation", "new_sigma", "sigma->val[cnt]");

	sigma->col = coord->sigma_col;
	sigma->cnt = num_pi;

	return sigma;
}

/***********************************************************************\
** This function frees all the data associated with the sigma
 ** structure, including every Pi x Tbar vector.  It then frees
 ** the structure itself.  Note it does NOT free sigma->col, since
 ** this is contained permanently in the prob structure.
 \***********************************************************************/
void free_sigma(sigma_type *sigma)
{
	int cnt;

#ifdef TRACE
	printf("Inside free_sigma\n");
#endif

	mem_free(sigma->lamb);
	mem_free(sigma->ck);
	//added by Yifan to clean iteration number(c->k) records
	for (cnt = 0; cnt < sigma->cnt; cnt++)
		mem_free(sigma->val[cnt].T);
	mem_free(sigma->val);
	mem_free(sigma);
}

/***********************************************************************\
** This function writes the entire sigma structure to a file.
 ** It may be used as the basis for sending cell messages in a
 ** parallel version of the algorithm.
 \***********************************************************************/
void write_sigma(FILE *fptr, sigma_type *sigma, num_type *num)
{
	int idx, cnt;

	fprintf(fptr, "SIGMA\n");
	fprintf(fptr, "CNT\n");
	fprintf(fptr, "%d \n", sigma->cnt);

	/* Print the columns of each entry in Pi x T */
	fprintf(fptr, "COLS");
	for (idx = 0; idx < num->nz_cols; idx++)
	{
		if (!(idx % NUM_INTS))
			fprintf(fptr, "\n");
		fprintf(fptr, "%d ", sigma->col[idx]);
	}

	/* Print the lambda associated with each sigma entry */
	fprintf(fptr, "\nLAMB");
	for (cnt = 0; cnt < sigma->cnt; cnt++)
	{
		if (!(cnt % NUM_INTS))
			fprintf(fptr, "\n");
		fprintf(fptr, "%d ", sigma->lamb[cnt]);
	}

	/* Print the contents of the sigma matrix */
	fprintf(fptr, "\nVAL");
	for (cnt = 0; cnt < sigma->cnt; cnt++)
	{
		fprintf(fptr, "\n%lf", sigma->val[cnt].R);
		for (idx = 0; idx < num->nz_cols; idx++)
		{
			if (!(idx % NUM_DBLS))
				fprintf(fptr, "\n");
			fprintf(fptr, "%lf ", sigma->val[cnt].T[idx]);
		}
		fprintf(fptr, "\n");
	}
}

