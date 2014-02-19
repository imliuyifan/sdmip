/***********************************************************************\
**
 ** lambda.c
 **
 **
 **
 ** It is assumed that the zeroth element of every vector contains
 ** the 1-norm of that vector.  If the 1-norm is never necessary,
 ** the zeroth element may be un-initialized; however, the vector
 ** must still start at 1. 
 **
 **
 ** calc_lambda()
 ** print_lambda()
 ** new_lambda()
 ** free_lambda()
 **
 ** History:
 **   ?? Nov 1991 - <Jason Mai> - created.
 **   16 Dec 1991 - <Jason Mai> - fixed compilation errors.
 **   18 Feb 1992 - <Jason Mai> - updated with new structures (esp. num_type). 
 **   Should NOT be usewd without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "utility.h"
#include "lambda.h"
#include "log.h"
#include "sdglobal.h"

/***********************************************************************\
** This function stores a new lambda_pi vector in the lambda
 ** structure.  Each lambda_pi represents only those dual variables 
 ** whose rows in the constraint matrix have random elements.  Thus 
 ** the (full) dual vector, Pi,  passed to the function is converted 
 ** into the sparse vector lambda_pi.  This vector is then compared with
 ** all previous lambda_pi vectors, searching for a duplication.
 ** If a duplicate is found, the vector is not added to the 
 ** structure, and the function returns the index of the duplicate 
 ** vector.  Otherwise, it adds the vector to the end of the structure, 
 ** and returns an index to the last element in lambda.
 \***********************************************************************/
int calc_lambda(sdglobal_type* sd_global, lambda_type *lambda, num_type *num,
		vector Pi, BOOL *new_lamb)
{
	int pi_idx; /* steps through lambda structure */
	int length; /* length of each lambda_pi vector */
	double *lambda_pi; /* sparse vector of dual variables */

#ifdef TRACE
	printf("Inside calc_lambda\n");
#endif

	length = num->rv_rows;

	/* Pull out only those elements in dual vector which have rv's */
	lambda_pi = reduce_vect(Pi, lambda->row, length);

	/* Compare resulting lambda_pi with all previous vectors */
	for (pi_idx = 0; pi_idx < lambda->cnt; pi_idx++)
		if (equal_arr(lambda_pi, lambda->val[pi_idx], length,
				sd_global->config.TOLERANCE))
		{
			mem_free(lambda_pi);
			*new_lamb = FALSE;
			return pi_idx;
		}

	/* Add the vector to lambda struct */
	lambda->val[lambda->cnt] = lambda_pi;

	*new_lamb = TRUE;
	return lambda->cnt++;
}

/***********************************************************************\
**
 \***********************************************************************/
void print_lambda(lambda_type *lambda, num_type *num, int idx)
{
	int cnt;

	printf("\nLambda %d:: ", idx);
	for (cnt = 0; cnt <= num->rv_rows; cnt++)
		printf("%f  ", lambda->val[idx][cnt]);
	printf("\nLambda Rows: ");
	for (cnt = 0; cnt <= num->rv_rows; cnt++)
		printf("%d ", lambda->row[cnt]);
	printf("\n");
}

/***********************************************************************\
** This function allocates a new lambda structure, with room
 ** for num_lambdas lambda vectors of size vect_size.  It returns a 
 ** pointer to the structure.  Only some of the individual lambda vectors 
 ** are expected to be allocated (according to the num_vect parameter)
 ** so that there is room for new lambdas to be created. 
 \***********************************************************************/
lambda_type *new_lambda(int num_iter, int num_lambda, int num_rv_rows,
		coord_type *coord)
{
	lambda_type *lambda;
	int cnt;

#ifdef TRACE
	printf("Inside new_lambda\n");
#endif

	if (!(lambda = (lambda_type *) mem_malloc (sizeof(lambda_type))))
		err_msg("Allocation", "new_lambda", "lambda");

	if (!(lambda->val = arr_alloc(num_iter, vector)))
		err_msg("Allocation", "new_lambda", "lambda->val");

	for (cnt = 0; cnt < num_lambda; cnt++)
		if (!(lambda->val[cnt] = arr_alloc(num_rv_rows+1, double)))
			err_msg("Allocation", "new_lambda", "lambda->val[cnt]");

	lambda->cnt = num_lambda;
	lambda->row = coord->lambda_row;

	return lambda;
}

/***********************************************************************\
** Note that lambda->row is not freed, since it belongs to someone else.
 \***********************************************************************/
void free_lambda(lambda_type *lambda)
{
	int cnt;

#ifdef TRACE
	printf("Inside free_lambda\n");
#endif

	for (cnt = 0; cnt < lambda->cnt; cnt++)
		mem_free(lambda->val[cnt]);

	mem_free(lambda->val);
	mem_free(lambda);
}

/***********************************************************************\
** This function writes the entire lambda structure to a file.
 ** It may be used as the basis for sending cell messages in a
 ** parallel version of the algorithm.
 \***********************************************************************/
void write_lambda(FILE *fptr, lambda_type *lambda, num_type *num)
{
	int idx, cnt;

	fprintf(fptr, "LAMBDA\n");
	fprintf(fptr, "CNT\n");
	fprintf(fptr, "%d \n", lambda->cnt);

	/* Print the rows in which lambda(pi) elements occur */
	fprintf(fptr, "ROWS");
	for (idx = 0; idx < num->rv_rows; idx++)
	{
		if (!(idx % NUM_INTS))
			fprintf(fptr, "\n");
		fprintf(fptr, "%d ", lambda->row[idx]);
	}

	/* Print the contents of the lambda collection */
	fprintf(fptr, "\nVAL");
	for (cnt = 0; cnt < lambda->cnt; cnt++)
	{
		for (idx = 0; idx < num->rv_rows; idx++)
		{
			if (!(idx % NUM_DBLS))
				fprintf(fptr, "\n");
			fprintf(fptr, "%lf ", lambda->val[cnt][cnt]);
		}
		fprintf(fptr, "\n");
	}
}

