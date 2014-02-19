/***************************************************************************\
**
 ** memory.c
 **
 **  This file contains routines for keeping track of and reducing
 **  the amount of memory in use by the SD algorithm. 
 **
 ** History:
 **   28 Oct 92 - <Jason Mai> - created for the second time.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 \***************************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "utility.h"
#include "delta.h"
#include "memory.h"
#include "omega.h"
#include "log.h"
#include "cuts.h"
#include "sdglobal.h"
#include <time.h>
#include <limits.h>

/***************************************************************************\
** This function prints statistics about the SD solution process to the
 ** screen and to various data files.
 \***************************************************************************/
void write_statistics(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
#ifdef WRITE
	FILE *f_cnt;
    FILE *f_iter;
    FILE *f_i;
	int cnt;
#endif

//	double time_diff;
	/* int		cnt, idx;*/
	clock_t cur_clock;
	double conf_int[2];
	/*added by Yifan to print out 1st stage solutions xlm files*/
	/*
	 char		mname[20] = "m   .xlm";
	 char    sname[20] = "s   .xlm";
	 static int    mnum = 0;
	 static int    snum = 0;
	 */

#ifdef TRACE
	printf("Inside write_statistics()\n");
#endif

	/* Find the time taken for this iteration */
	cur_clock = clock();
//	time_diff = cur_clock - sd_global->LAST_CLOCK;
	sd_global->LAST_CLOCK = cur_clock;

#ifdef WRITE
	/* Write out the size of each structure, plus memory allocated */
	f_cnt = fopen(CNT_DAT, "a");
	fprintf(f_cnt, "%d, %lld, %d, %d, %d, %d, %lf;\n", c->k, sd_global->MEM_USED,
			c->sigma->cnt, c->lambda->cnt, s->omega->cnt, c->cuts->cnt,
			time_diff / 1000000);
	fclose(f_cnt);
#endif

#ifdef RUN
	if (!(c->k % sd_global->config.PRINT_CYCLE))
	{
		/* Write information about the solution to the screen */
		/* write_iter(stdout, p, c, s); */
		printf("\nIteration %-4d:", c->k);
	}
#endif


	if (!(c->k % sd_global->config.PRINT_CYCLE)
			&& sd_global->config.EVAL_RUN_FLAG == 1)
	{
		evaluate_inc(sd_global, c, p, s, s->incumb_x, NULL, conf_int, 0);
	}

#ifdef WRITE
	if(!(c->k % sd_global->config.PRINT_CYCLE))
	{
		/* Write the same stuff to a data file */
		f_iter = fopen(ITER_DAT, "a");
		write_iter(f_iter, p, c, s);
		fclose(f_iter);
	}
    f_i = fopen("ITERD", "a");
    fprintf(f_i, "%f,%f\n",s->incumb_est, s->pi_ratio[c->k%sd_global->config.SCAN_LEN]);
    fclose(f_i);
#endif

#ifdef TRACE
	printf("Exiting write_statistics()\n");
#endif
}

/***************************************************************************\
** This function reduces the amount of data being stored by various 
 ** structures of the SD algorithm.  Each one (cuts, omega, lambda, 
 ** sigma, and delta) is analyzed for redundant, irrelevant, or worthless
 ** data which is thinned out according to various tolerances.
 \***************************************************************************/
void thin_data(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{

#ifdef TRACE
	printf("Inside thin_data()\n");
#endif

	/* Reset the amount of memory allocated during this iteration */
	sd_global->MEM_USED = 0;

	/* If its possible to thin out any cuts, then try it */
	if (c->k >= sd_global->config.DROP_TIME)
		thin_cuts(sd_global, p, c, s);

	/* If its time to reduce the data structures, do so */
	if (c->k >= sd_global->config.START_THIN
			&& !((c->k - sd_global->config.START_THIN)
					% sd_global->config.THIN_CYCLE))
	{
		thin_pi(p, c, s);

		thin_omega(sd_global, p, c, s);

	}

#ifdef TRACE
	printf("Exiting thin_data()\n");
#endif
}

/***************************************************************************\
** This function reduces the number of dual vectors stored in the
 ** lambda, sigma, and delta structures.  It drops any dual vector which 
 ** has not been used in the formation of a single one of the current cuts 
 ** (that is, it hasn't been selected as an istar for any cut).
 \***************************************************************************/
void thin_pi(prob_type *p, cell_type *c, soln_type *s)
{
	int *sig_histo, *lamb_histo;
	int cnt, idx, pi_idx;
	int sig_cnt = 0;
	int lamb_cnt = 0;
//	int scnt;

#ifdef TRACE
	printf("Inside thin_pi()\n");
#endif

	/* Allocate arrays of counters for every Pi in sigma and lambda */
	if (!(sig_histo = arr_alloc(c->sigma->cnt, int)))
		err_msg("Allocation", "thin_pi()", "sig_histo");
	if (!(lamb_histo = arr_alloc(c->lambda->cnt, int)))
		err_msg("Allocation", "thin_pi()", "lamb_histo");

	/* Loop through all the cuts to count the number of times each Pi is used */
	for (cnt = 0; cnt < c->cuts->cnt; cnt++)
	{
		for (idx = 0; idx < c->cuts->val[cnt]->omega_cnt; idx++)
		{
			pi_idx = c->cuts->val[cnt]->istar[idx];

#ifdef PIS
			if (pi_idx < 0 || pi_idx >= c->sigma->cnt)
			printf("!!! pi_idx = %d !!!\n", pi_idx);

			if (c->sigma->lamb[pi_idx] < 0 || c->sigma->lamb[pi_idx] >= c->lambda->cnt)
			printf("!!! lamb[pi_idx] = %d !!!\n", c->sigma->lamb[pi_idx]);
#endif

			/* Don't count -1 istars */
			if (pi_idx >= 0)
			{
				++sig_histo[pi_idx];
				++lamb_histo[c->sigma->lamb[pi_idx]];
			}
		}
	}

//	scnt = c->sigma->cnt;
//	lcnt = c->lambda->cnt;

#ifdef PIS
	/*
	 printf("Sigma histo: ");
	 for (cnt = 0; cnt < c->sigma->cnt; cnt++)
	 printf("%d (%d): ", sig_histo[cnt], c->sigma->lamb[cnt]);
	 printf("\n");

	 printf("Lambda histo: ");
	 for (cnt = 0; cnt < c->lambda->cnt; cnt++)
	 printf("%d, ", lamb_histo[cnt]);
	 printf("\n");
	 */
#endif

	/* Remove any Pi in sigma whose count is zero */
	for (pi_idx = c->sigma->cnt - 1; pi_idx >= 0; pi_idx--)
	{
		if (sig_histo[pi_idx] == 0)
		{
			++sig_cnt;
			drop_sigma(c->sigma, c->cuts, pi_idx);
		}
	}

	/* Remove any Pi in lambda whose count was zero */
	for (pi_idx = c->lambda->cnt - 1; pi_idx >= 0; pi_idx--)
	{
		if (lamb_histo[pi_idx] == 0)
		{
			++lamb_cnt;
			drop_lambda(c->lambda, s->delta, s->omega, c->sigma, pi_idx);
		}
	}

#ifdef RUN
	printf("thin_pi() dropped %d sigmas and %d lambdas.\n", sig_cnt, lamb_cnt);
#endif

#ifdef WRITE
	/* Print the histogram arrays to output files */
	write_histo_file(sig_histo, scnt, c->k, S_HIST_DAT);
	write_histo_file(lamb_histo, lcnt, c->k, L_HIST_DAT);
#endif

	mem_free(sig_histo);
	mem_free(lamb_histo);
}

/***************************************************************************\
** This function steps through the observations in omega and eliminates
 ** any which do not seem necessary.  This involves freeing the entry
 ** in omega and eliminating the correpsonding columns of delta.
 \***************************************************************************/
void thin_omega(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	int obs, dupl;
	int omeg_cnt = 0;

#ifdef TRACE
	printf("Inside thin_omega()\n");

#endif

	/* Loop through the latest observations of omega, starting at the end */
	for (obs = s->omega->most - 1; obs >= s->omega->last; obs--)
	{

#ifdef OMEG
		printf("Checking omega %d.\n");
#endif

		/* First make sure the observation hasn't been dropped already */
		if (valid_omega_idx(s->omega, obs))
		{
			/* Check whether this observation has been duplicated by another */
			if (duplicate_omega(sd_global, p, c, s, obs, &dupl))
			{
				/* Drop this omega in favor of its duplicate */
				drop_omega(s->omega, s->delta, c->lambda, c->cuts, obs, dupl);
				omeg_cnt++;
			}
		}
	}
	/* Next time we'll stop where we started this time */
	s->omega->last = s->omega->most;

#ifdef RUN
	printf("thin_omega() dropped %d omegas.\n", omeg_cnt);
#endif
}

/***************************************************************************\
** This function drops a designated row from the sigma structure.  The
 ** row is freed, the count is decremented, and the last entry in sigma
 ** is swapped into the vacated position.  It is possible that some cuts
 ** are using this sigma as an istar, but we ignore that for now, since
 ** we know that the rule in memory.c only drops Pi's which have NO istars.
 \***************************************************************************/
void drop_sigma(sigma_type *sigma, cut_type *cuts, int idx)
{
	int cnt, obs;

#ifdef TRACE
	printf("Inside drop_sigma()\n");
#endif

	/* One less sigma exists */
	--sigma->cnt;

	/* Free the data associated with the dropped one */
	mem_free(sigma->val[idx].T);

	/* Swap the last entry into the open position */
	sigma->lamb[idx] = sigma->lamb[sigma->cnt];
	sigma->val[idx] = sigma->val[sigma->cnt];

	/* Update the istars of any cut which referenced the swapped sigma */
	for (cnt = 0; cnt < cuts->cnt; cnt++)
		for (obs = 0; obs < cuts->val[cnt]->omega_cnt; obs++)
			if (cuts->val[cnt]->istar[obs] == sigma->cnt)
				cuts->val[cnt]->istar[obs] = idx;

	/* Here you should do something to fix up the istars of dropped sigma */
}

/***************************************************************************\
** This function drops a designated row from the lambda structure.  The
 ** specified row is freed, the lambda count is decremented, and the last
 ** entry in lambda is swapped into the vacated position.  It is likely
 ** that some entry in sigma was associated with the lambda that got 
 ** swapped, so this function searches the sigma structure and updates
 ** any sigma->lamb references to it.  Finally, the dropped row in lambda
 ** corresponds (1-to-1) to a row in delta, so this function calls on
 ** drop_delta_row() to free it.
 \***************************************************************************/
void drop_lambda(lambda_type *lambda, delta_type *delta, omega_type *omega,
		sigma_type *sigma, int idx)
{
	int cnt;

#ifdef TRACE
	printf("Inside drop_lambda()\n");
#endif

	/* One less lambda exists */
	--lambda->cnt;

	/* Free the designated lambda vector */
	mem_free(lambda->val[idx]);

	/* Swap the last entry in lambda into the emptied position */
	lambda->val[idx] = lambda->val[lambda->cnt];

	/* Update the rest of the world, to make it look like lambda never existed */
	drop_delta_row(delta, lambda, omega, idx);

	/* Change any references in sigma to the swapped entry in lambda */
	for (cnt = 0; cnt < sigma->cnt; cnt++)
		if (sigma->lamb[cnt] == lambda->cnt)
			sigma->lamb[cnt] = idx;

	/* Change any references in sigma to the dropped entry in lambda */
	/* Sigma always gets dropped if lambda does, so ignore this */
}

/***************************************************************************\
** This function removes an entry from the omega structure, in favor of
 ** another entry in the structure.  It deletes the omega row, increments
 ** the weight on the favored omega, changes the omega->filter array,
 ** and removes the corresponding column from delta.  Any cuts which had
 ** istar entries for the dropped omega are filled with a -1 in that position.
 \***************************************************************************/
void drop_omega(omega_type *omega, delta_type *delta, lambda_type *lambda,
		cut_type *cuts, int drop, int keep)
{
	int cnt;

#ifdef TRACE
	printf("Inside drop_omega()\n");
#endif

#ifdef OMEG
	printf("Dropping omega %d in favor of %d.\n", drop, keep);
#endif

	/* Free the entry in omega */
	mem_free(omega->idx[drop]);

	/* Free the column in delta */
	drop_delta_col(delta, lambda, drop);

	/* Fix the cuts' istars, so the dropped omega is erased */
	for (cnt = 0; cnt < cuts->cnt; cnt++)
		if (drop < cuts->val[cnt]->omega_cnt)
			cuts->val[cnt]->istar[drop] = DROPPED;

	/* Update the counters of the omega structure */
	omega->weight[keep] += omega->weight[drop];
	omega->filter[drop] = UNUSED;
	omega->weight[drop] = 0;
	if (omega->next > drop)
		omega->next = drop;
	omega->cnt--;
}

/***************************************************************************\
** This function determines whether the specified observation of omega
 ** is essentially a duplicate of some previous observation.  It returns
 ** TRUE if (according to some criteria) it's a duplicate; FALSE otherwise.
 ** It initializes the _dupl_ parameter with the index of the duplicate.
 \***************************************************************************/
BOOL duplicate_omega(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int obs, int *dupl)
{
#ifdef OMEG
	int cnt;
#endif
	int omeg_idx;
	int *marked1, *marked2;

	/* Allocate arrays large enough to flag any pi in sigma */
	if (!(marked1 = arr_alloc(c->sigma->cnt, int)))
		err_msg("Allocation", "duplicate_omega", "marked1");
	if (!(marked2 = arr_alloc(c->sigma->cnt, int)))
		err_msg("Allocation", "duplicate_omega", "marked2");

	/* Mark all of the pi's used as istars for omega _obs_ */
	mark_used_pis(marked1, c->cuts, obs);

#ifdef OMEG
	printf("Marked1 array: ");
	for (cnt = 0; cnt < c->sigma->cnt; cnt++)
	printf("%d ", marked1[cnt]);
	printf(" MMM\n");
#endif

	/* Loop through all the previous omegas */
	for (omeg_idx = obs - 1; omeg_idx >= 0; omeg_idx--)
	{
		if (valid_omega_idx(s->omega, omeg_idx))
		{
			/* Copy the original marked array of the _obs_ omega */
			copy_arr(marked2, marked1, c->sigma->cnt - 1);

			/* Add marks for all pi's used as istars of the current observation */
			mark_used_pis(marked2, c->cuts, omeg_idx);

			/* Compare the two observations in the marked rows of delta */
			if (duplic_delta_col(sd_global, s->delta, c->sigma, p->num, marked2,
					omeg_idx, obs))
			{
				*dupl = omeg_idx;
				return TRUE;
			}
		}
	}

	mem_free(marked1);
	mem_free(marked2);
	*dupl = -1;

	return FALSE;
}

/***************************************************************************\
** This function determines whether two columns in delta are essentially
 ** the same by checking the rows specified by in _marked_ array.
 ** It returns TRUE if the two columns were within the tolerance for
 ** every row in the marked array; FALSE otherwise.
 \***************************************************************************/
BOOL duplic_delta_col(sdglobal_type* sd_global, delta_type *delta,
		sigma_type *sigma, num_type *num, int *marked, int omeg_idx, int obs)
{
	int cnt, pi;

#ifdef LOOP
	printf("Inside duplic_delta_col()\n");
#endif

	/* Loop through all the rows in the marked array */
	for (cnt = 0; cnt < sigma->cnt; cnt++)
		if (marked[cnt] == USED)
		{
			/* Locate row in delta */
			pi = sigma->lamb[cnt];

			/* Check the Pi x Romega field for significant difference */
			if (DBL_ABS(delta->val[pi][omeg_idx].R - delta->val[pi][obs].R)
					> DBL_ABS(sd_global->config.THIN_TOLER * delta->val[pi][obs].R))
				return FALSE;

			/* Check the Pi x Tomega vector for significant difference */
			if (!equal_arr(delta->val[pi][obs].T, delta->val[pi][omeg_idx].T,
					num->rv_cols, sd_global->config.THIN_TOLER))
				return FALSE;
		}

	return TRUE;
}

/***************************************************************************\
** This function tallies those pi's which are referenced by the istar 
 ** field of the current cuts.  For each cut, it evaluates the _obs_ entry
 ** in the istar array and marks that field in the _marked_ array as USED.
 ** Entries in the marked array which aren't used as istars are left alone.
 \***************************************************************************/
void mark_used_pis(int *marked, cut_type *cuts, int obs)
{
	int cnt, pi;

#ifdef LOOP
	printf("Inside mark_used_pis()\n");
#endif

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		if (obs < cuts->val[cnt]->omega_cnt)
		{
			pi = cuts->val[cnt]->istar[obs];
			if (pi >= 0)
				marked[pi] = USED;
		}
}

/***************************************************************************\
** This function prints out a histogram array to a specified file.
 ** It is printed in matlab plot-ready format.
 \***************************************************************************/
void write_histo_file(int *histo, int arr_cnt, int iter, char *fname)
{
	int cnt;
	FILE *f_out;

#ifdef TRACE
	printf("Inside write_histo_file()\n");
#endif

	f_out = fopen(fname, "w");
	if (f_out == NULL)
	{
		printf("!!! Error printing histogram data !!!\n");
		return;
	}

	/*
	 fprintf(f_out, "%% Iteration #%d\n", iter);
	 fprintf(f_out, "array = [\n");
	 */
	for (cnt = 0; cnt < arr_cnt; cnt++)
		fprintf(f_out, "%d;\n", histo[cnt]);
	/*
	 fprintf(f_out, "]\n");
	 */
	fclose(f_out);
}

/***************************************************************************\
** This function prints out iteration data about the SD solution process.
 ** This includes the objective function estimates, gamma, and the
 ** candidate & incumbent solution vectors.
 \***************************************************************************/
void write_iter(sdglobal_type* sd_global, FILE *fout, prob_type *p,
		cell_type *c, soln_type *s)
{

#ifdef TRACE
	printf("Inside write_iter()\n");
#endif

	if (fout == NULL)
	{
		printf("!!! Error printing iteration data !!!\n");
		return;
	}

	fprintf(fout, "\nIteration %d.  Gamma = %lf\n", c->k, s->gamma);

	/* Added to assure the validity of zero as subproblem LB.
	 * Lei 09/13/05. */
	if (sd_global->config.SUB_LB_CHECK)
		fprintf(fout, "Subproblem LB = %lf\n", s->sub_lb_checker);

	fprintf(fout, "f_k(candid_x) = %lf, f_k(incumb_x) = %lf\n", s->candid_est,
			s->incumb_est);
	/* Write first stage decisions */
	/*
	 fprint_vect(fout, s->candid_x, p->num->mast_cols, "Candidate X");
	 fprint_vect(fout, s->incumb_x, p->num->mast_cols, "Incumbent X");
	 fprintf(fout, "\n");
	 */

	fprintf(fout, "COUNTS: Lambda=%d.  Sigma=%d.  Omega=%d  Cuts=%d.\n",
			c->lambda->cnt, c->sigma->cnt, s->omega->cnt, c->cuts->cnt);
}
