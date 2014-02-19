/***********************************************************************\
**
 ** soln.c
 **
 ** This file contains the functions used in preparing for a new solution.
 ** The omega and delta structures, which depend upon realizations of the
 ** random variables, are needed temporarily.  While a cell is being solved,
 ** they are stored in _soln_, but may be discarded afterwards.
 **
 **
 ** History:
 **   16 Mar 1992 - <Jason Mai> - created.
 **   20 Mar 1992 - <Jason Mai> - debugged.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "utility.h"
#include "delta.h"
#include "omega.h"
#include "log.h"
#include "quad.h"
#include "cuts.h"
#include "testout.h"
#include "improve.h"
#include "master.h"
#include "subprob.h"
#include "solver.h"
#include "batch.h"
#include "sdglobal.h"

/***********************************************************************\
** This function finds all rows from the subproblem (R and T) which
 ** contain random elements.  To accomplish this, it makes a set of 
 ** all elements in the omega.row array, such that repeated elements
 ** only occur once.  This resulting array (set) is returned. 
 ** These rows will become the coordinates of lambda(pi).
 **
 ** This assumes that get_omega_row will fill in the 1-norm. 
 ** Remember the 1-norm when trying to figure this code out!
 \***********************************************************************/
int *find_rows(int num_elem, int *num_rows, int *omega_row)
{
	int *rv_row;
	int i, j;
	int len;

#ifdef TRACE
	printf("Inside find_rows\n");
#endif

	len = 0;
	if (!(rv_row = arr_alloc(num_elem+1, int)))
		err_msg("Allocation", "find_rows", "rv_row");

	/* Copy over distinct elements of omega.row */
	for (i = 0; i < num_elem; i++)
	{
		for (j = 0; j < len && omega_row[i + 1] != rv_row[j]; j++)
			/* no loop body */;
		if (j == len)
			rv_row[++len] = omega_row[i + 1];
	}

	/* Shrink the array down to the number of distinct elements found */
	rv_row = (int *) mem_realloc(rv_row, (len+1)*sizeof(int));

	rv_row[0] = 0;
	*num_rows = len;

	return rv_row;
}

/***********************************************************************\
** This function is very similar to find_rows, in that it may be used to 
 ** find columns in the Tomega matrix which contain random elements.  
 ** However, it is also used to find non-zero columns in the Tbar matrix.
 ** It excludes column 0, since this is in the R vector, not the T matrix.
 ** Also, it takes advantage of the known fact that the elements of 
 ** T have been sorted in increasing order by column.  Thus, to 
 ** determine the distinct-ness of a given element, we need only to 
 ** compare it with the most recent element.  It returns the resulting
 ** array.  These columns will become the coordinates of Pi x Tomega.
 \***********************************************************************/
int *find_cols(int num_elem, int *num_cols, int *omega_col)
{
	int *rv_col;
	int i;
	int len;
	int last;

#ifdef TRACE
	printf("Inside find_cols\n");
#endif

	len = 0;
	if (!(rv_col = arr_alloc(num_elem+1, int)))
		err_msg("Allocation", "find_cols", "rv_col");

	last = 0; /* If the last observation was 0, we won't store new ones */

	/* Copy over all the distinct elements of omega_col */
	for (i = 0; i < num_elem; i++)
		if (omega_col[i + 1] != last)
			last = rv_col[++len] = omega_col[i + 1];

	/* Shrink the array down to the number of distinct elements found */
	rv_col = (int *) mem_realloc(rv_col, (len+1)*sizeof(int));

	rv_col[0] = 0;
	*num_cols = len;

	return rv_col;
}

/***********************************************************************\
** This function allocates memory for a new soln data structure.
 ** It's fields are also allocated, and initialized.
 \***********************************************************************/
soln_type * new_soln(sdglobal_type* sd_global, prob_type *p, vector x_k)
{
	soln_type *s;
	int length;
	int i;
	double Sm;

#ifdef TRACE
	printf("Inside new_soln\n");
#endif

	if (!(s = (soln_type *) mem_malloc (sizeof(soln_type))))
		err_msg("Allocation", "new_soln", "s");

	/* Initialize the time_type struct to record the time spent on different
	 procedures of the SD algorithm. added by zl, 06/29/04. */
	if (!(s->run_time = (time_type *) mem_malloc (sizeof(time_type))))
		err_msg("Allocation", "time_type", "s->run_time");
	if (!(s->Batch_pi = (vector *) mem_calloc (BATCH_SIZE, sizeof(vector))))
		err_msg("Allocation", "new_soln", "Batch_pi");
	for (i = 0; i < BATCH_SIZE; i++)
	{
		s->Batch_pi[i] = arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	}
	if (!(s->Batch_dj = (vector *) mem_calloc (BATCH_SIZE, sizeof(vector))))
		err_msg("Allocation", "new_soln", "Batch_pi");
	for (i = 0; i < BATCH_SIZE; i++)
	{
		s->Batch_dj[i] = arr_alloc(p->num->mast_cols+2,double);
	}

	if (!(s->xc_height = arr_alloc(BATCH_SIZE, double)))
		err_msg("Allocation", "new_soln", "xc_height");
	/*beta is used to store cut coefficients*/
	if (!(s->beta = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "new_soln", "beta");

	if (!(s->dual_statble_flag = (BOOL *) mem_malloc (sizeof(BOOL))))
		err_msg("Allocation", "BOOL_type", "s->dual_statble_flag");
	if (!(s->pi_ratio = arr_alloc (sd_global->config.SCAN_LEN, double)))
		err_msg("Allocation", "double_type", "s->pi_ratio");

	/*added by Yifan, 09/29/2011*/
	for (i = 0; i < sd_global->config.SCAN_LEN; i++)
		s->pi_ratio[i] = 0;

	s->run_time->total_time = 0.0;
	s->run_time->iteration_time = 0.0;
	s->run_time->iteration_accum = 0.0;
	s->run_time->soln_master_iter = 0.0;
	s->run_time->soln_subprob_iter = 0.0;
	s->run_time->full_test_iter = 0.0;
	s->run_time->argmax_iter = 0.0;
	s->run_time->soln_master_accum = 0.0;
	s->run_time->soln_subprob_accum = 0.0;
	s->run_time->full_test_accum = 0.0;
	s->run_time->argmax_accum = 0.0;
	s->max_ratio = 0.0;
	s->min_ratio = 1.0;

	/* Initialize the boolean _optimal_flag_. zl, 06/30/04. */
	s->optimality_flag = FALSE;
	s->smpl_test_flag = FALSE; /* 08/17/04. */
	s->smpl_ever_flag = FALSE; /* 08/20/04. */
	*s->dual_statble_flag = FALSE; /*added by Yifan 09/27/2011*/
	s->incumbent_change = FALSE;

	/* Initialize the lower bound checker of the subproblem objective
	 function values. zl, 07/01/04. */
	s->sub_lb_checker = sd_global->config.SUBPROB_LB;

	length = p->num->iter + p->num->iter / p->tau + 1;
	s->omega = new_omega(p->num->iter, p->num->rv, p->coord);
	s->delta = new_delta(length, p->coord);

	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	s->feasible_delta = new_delta(length, p->coord);

	/* Make initial allocation of the x vectors -- not freed until the end */
	s->incumb_x = duplic_arr(x_k, p->num->mast_cols);
	s->incumb_d = duplic_arr(x_k, p->num->mast_cols); /* modified by Yifan 2012.09.23 */
	s->incumb_avg = duplic_arr(x_k, p->num->mast_cols); /* modified by Yifan 2012.09.23 */
	s->candid_x = duplic_arr(x_k, p->num->mast_cols);
    s->mip_x = duplic_arr(x_k, p->num->mast_cols);  /* modified by Yifan 2013.9.15 */

	/* Initial allocation for the subproblem dual vector, including 1-norm */
	s->Pi = arr_alloc(p->num->sub_rows+1, double);
	/* Allocating space for the master problem dual vector & slacks JH 4/8/98 */
	/* Yifan 03/09/2012 Master_pi structure need to be updated at the begining of the optimality full test*/
	s->Master_pi = arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	s->Master_dj = arr_alloc(p->num->mast_cols+2,double);

	/* modified by Yifan 2012.10.05 */
	s->Batch_pi_mean = arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	s->Batch_dj_mean = arr_alloc(p->num->mast_cols+2,double);
	s->Batch_pi_stdev = arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	s->Batch_dj_stdev = arr_alloc(p->num->mast_cols+2,double);
	/* modified by Yifan 2012.10.05 */
	s->R_Master_pi_mean = arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	s->R_Master_dj_mean = arr_alloc(p->num->mast_cols+2,double);
	s->R_Master_pi_stdev =
			arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	s->R_Master_dj_stdev = arr_alloc(p->num->mast_cols+2,double);

	/* Yifan one for zero norm and another one for eta column */

	s->incumb_cut = 0;
	s->last_update = 0;

	if (sd_global->config.LB_TYPE == 0)
	{
		Sm = 0.0;
	}
	else
	{
		Sm = sd_global->Eta0;
	}

	/* This can't be zero... maybe start with C x X_k ? */
	s->gamma = 0.0;
	s->candid_est = Sm + CxX(p->c, s->candid_x, p->num->mast_cols);
	s->incumb_est = s->candid_est;
	s->incumb_stdev = 0.0;

	return s;
}

/***********************************************************************\
** This function frees the structures contained in the soln
 ** data structure, and then frees the soln itself.
 \***********************************************************************/
void free_soln(prob_type *p, cell_type *c, soln_type *s)
{
	int i;
#ifdef TRACE
	printf("Inside free_soln\n");
#endif
	free_delta(s->delta, s->omega, c->lambda->cnt);
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	free_delta(s->feasible_delta, s->omega, c->feasible_lambda->cnt);
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	free_omega(s->omega);
	mem_free(s->run_time);
	/* added by zl, 06/29/04. */
	mem_free(s->dual_statble_flag);
	/* added by Yifan */
	mem_free(s->incumb_x);
	mem_free(s->incumb_d);
	/* modified by Yifan 2012.09.23 */
	mem_free(s->incumb_avg);
	/* modified by Yifan 2012.09.23 */
	mem_free(s->beta);
	/* modified by Yifan 2012.09.28 */
	mem_free(s->xc_height);
	/* modified by Yifan 2012.09.28 */

	for (i = 0; i < BATCH_SIZE; i++)
	{
		mem_free(s->Batch_pi[i]);
		mem_free(s->Batch_dj[i]);
	}mem_free(s->Batch_pi);
	/* modified by Yifan 2012.09.28 */
	mem_free(s->Batch_dj);
	/* modified by Yifan 2012.09.28 */
	mem_free(s->candid_x);
    mem_free(s->mip_x);
	mem_free(s->pi_ratio);
	/* added by Yifan, 09/27/2011*/
	mem_free(s->Pi);
	mem_free(s->Master_pi);
	mem_free(s->Master_dj);

	/* modified by Yifan 2012.10.05 */
	mem_free(s->Batch_pi_mean);
	mem_free(s->Batch_dj_mean);
	mem_free(s->Batch_pi_stdev);
	mem_free(s->Batch_dj_stdev);
	/* modified by Yifan 2012.10.05 */
	mem_free(s->R_Master_pi_mean);
	mem_free(s->R_Master_dj_mean);
	mem_free(s->R_Master_pi_stdev);
	mem_free(s->R_Master_dj_stdev);

	mem_free(s);
}

/***********************************************************************\
** This function initializes the time_type struct, which is to record
 ** the time spent on differnt procedures of the algorithm. zl, 06/30/04
 \***********************************************************************/

/***********************************************************************\
** This function is called when the SLP has been solved.  It prints 
 ** the solution to the screen, and prints out other information to 
 ** an output file.  It has available all data (prob, cell, soln)
 ** associated with the problem and its solution. 
 \***********************************************************************/
void print_soln(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, char *fname)
{
	int i;
	char filename[NAME_SIZE * 2];
	char soln_file[NAME_SIZE * 2]; /* added by zl. */
	FILE *f_out;
	FILE *f_sol; /* added by zl. */

#ifdef TRACE
	printf("Inside print_soln\n");
#endif

	/* Print an answer to the screen */
	printf("\nSolution for %s:", fname);
	printf("\nX_opt = ");
	for (i = 0; i < p->num->mast_cols; i++)
		printf("%f ", s->incumb_x[i + 1]);
	printf("\n");

	/* Print the first stage solution to an output file */
	strcpy(soln_file, fname);
	strcat(soln_file, ".soln.out");

	if (!(f_sol = fopen(soln_file, "w")))
	{
		printf("Error writing output file %s.", soln_file);
		return;
	}

	fprintf(f_sol, "First stage solution :: \n");
	for (i = 1; i <= p->num->mast_cols; i++)
		fprintf(f_sol, " %s = %f\n", p->master->cname[i - 1], s->incumb_x[i]);

	fprintf(f_sol, "\nEstimated first stage objective function value: %f.\n",
			s->incumb_est);

	/* Print the statistic data to an output file */
	strcpy(filename, fname);
	strcat(filename, ".out"); /* modified by zl. */

	if (!(f_out = fopen(filename, "w")))
	{
		printf("Error writing output file %s.", filename);
		return;
	}
	/* added by zl, 08/02/04. */
	fprintf(f_out, "Solver \tSD\n\n");

	/* Be careful of the valid LB on the subproblem objective function
	 values. zl, 07/01/04. */
	if (s->sub_lb_checker < sd_global->config.SUBPROB_LB)
		fprintf(f_out, "Warning: the lowest subprob obj function value = %f\n",
				s->sub_lb_checker);
	else
		fprintf(f_out, "%f is a valid LB on subprob obj function values\n",
				sd_global->config.SUBPROB_LB);

	/* added by zl, 06/30/04. */
	fprintf(f_out, "Full test tolerance = %f\n", sd_global->config.EPSILON);
	if (s->smpl_ever_flag)
	{
		if (s->optimality_flag)
		{
			if (sd_global->config.TEST_TYPE == 0)
				fprintf(f_out, "Optimality test: Smpl_passed\n");
			else
				fprintf(
						f_out,
						"Optimality test: Full_passed  #passed: %d, #failed: %d.\n",
						s->passed, sd_global->config.M - s->passed);

			fprintf(f_out, "Inc_UB = %f, LB = %f, error = %f\n\n",
					s->incumb_est, s->candid_est,
					(s->incumb_est - s->candid_est) / s->incumb_est);
		}
		else
		{
			if (!(s->smpl_test_flag))
			{
				fprintf(f_out, "Optimality test: Smpl_ever\n");
				fprintf(f_out, "Inc_UB = %f, LB = %f, error = %f\n\n",
						s->incumb_est, s->candid_est,
						(s->incumb_est - s->candid_est) / s->incumb_est);
			}
			else
			{
				fprintf(f_out, "Optimality test: Full_failed\n");
				fprintf(f_out, "Inc_UB = %f, LB = %f, error = %f\n\n",
						s->incumb_est, s->candid_est, s->full_test_error);
			}
		}
	}
	else
	{
		fprintf(f_out, "Optimality test: Smpl_never\n");
		fprintf(f_out, "Inc_UB = %f, LB = %f, error = %f\n\n", s->incumb_est,
				s->candid_est, (s->incumb_est - s->candid_est) / s->incumb_est);
	}

	/* modified by zl. 
	 fprint_vect(f_out, s->incumb_x, p->num->mast_cols, "Final incumbent");
	 */
	fprintf(f_out, "First stage Solution :: ");
	for (i = 1; i <= p->num->mast_cols; i++)
		fprintf(f_out, "%f ", s->incumb_x[i]);
	fprintf(f_out, "\n");

	fprintf(f_out, "Estimated first stage objective: %lf\n\n", s->incumb_est);

	/*  fprintf(f_out, "Actual objective at incumbent: unknown\n"); */

	fprintf(f_out, "Number of iterations:   c->k        = %d\n", c->k);
	fprintf(f_out, "Number of cuts formed:  cuts->cnt   = %d\n", c->cuts->cnt);
	/* LP_cnt and LP_test added by zl. 06/30/02 */
	fprintf(f_out, "Number of LPs solved:   c->LP_cnt   = %d\n", c->LP_cnt);
	fprintf(f_out, "Number of LPs solved\n");
	fprintf(f_out, "  for testing:          c->LP_test  = %d\n", c->LP_test);
	fprintf(f_out, "Number of dual vectors: lambda->cnt = %d\n",
			c->lambda->cnt);
	fprintf(f_out, "Number of observations: omega->cnt  = %d\n", s->omega->cnt);
	fprintf(f_out, "Number of sigma pairs:  sigma->cnt  = %d\n", c->sigma->cnt);
	/* added by zl. 06/18/02 */
	fprintf(f_out,
			"Number of iterations since the last time the incumbent cut\n");
	fprintf(f_out, "  was updated:  c->k - s->incumb_k  = %d\n",
			c->k - s->incumb_k);

	/* Added by zl, 06/30/04. */
	fprintf(f_out, "\nTotal time \t\t\t= %.2f\n", s->run_time->total_time);
	fprintf(f_out, "Accumulated iteration time \t= %.2f\n",
			s->run_time->iteration_accum);
	fprintf(f_out, "Accumulated master soln time \t= %.2f\n",
			s->run_time->soln_master_accum);
	fprintf(f_out, "Accumulated subprob soln time \t= %.2f\n",
			s->run_time->soln_subprob_accum);
	fprintf(f_out, "Accumulated full test time \t= %.2f\n",
			s->run_time->full_test_accum);
	fprintf(f_out, "Accumulated argmax time \t= %.2f\n",
			s->run_time->argmax_accum);
	fclose(f_out);
	fclose(f_sol); /* added by zl. */

	/* Print the contents of the cell to a file */
	/*
	 write_cell(p, c, fname);
	 */

}

/*  status = 0 for each incumbent solution during replication
    status = 1 for compromise solution
    status = 2 for mean solution
*/
int print_detailed_soln(sdglobal_type* sd_global, soln_type *s, prob_type *p,
		char *fname, int status)
{
	FILE *fp;
	char fp_fname[NAME_SIZE * 2];
	int i, soln_status = 0, report_average = 0;
	double mean, ll, ul;
	vector A_x;
	vector x;
	vector pi_mean;
	vector dj_mean;
	vector pi_stdev;
	vector dj_stdev;
	strcpy(fp_fname, fname);

	if (!(A_x = arr_alloc(p->num->mast_rows + 1, double)))
		err_msg("Allocation", "print_soln", "A_x");

	if (status == 0)
	{
		strcat(fp_fname, ".detailed_rep_soln.out");
		fp = fopen(fp_fname, "a");
	}
	else
	{
		strcat(fp_fname, ".detailed_soln.out");
		fp = fopen(fp_fname, "a");
	}

	/* Print out problem name */
	fprintf(fp, "%-40s%s\n", "Problem:", fname == NULL ? "" : fname);
	/* Print out number of rows (including objective)*/
	fprintf(fp, "%-40s%d\n", "First Stage Rows:", p->num->mast_rows + 1);
	/* Print out number of columns (excluding rhs)*/
	fprintf(fp, "%-40s%d\n", "First Stage Columns:", p->num->mast_cols);
	/* Print out number of non-zeros */
	fprintf(fp, "%-40s%d\n", "First Stage Non-zeros:", p->A->cnt);

	if (status != 0)
	{
      if (sd_global->average_flag == 1) {
        /* report_average = test_average_soln(s, p); */
        report_average = 1;
      }
      if (report_average){
			fprintf(fp, "Mean solution is recommended for this instance.\n");
      }
	}

	if (status == 0)
	{
		soln_status = 3;
		mean = s->rep_mean;
		ul = s->rep_mean + 1.96 * s->rep_stdev;
		ll = s->rep_mean - 1.96 * s->rep_stdev;
		x = s->incumb_x;
		pi_mean = s->Master_pi;
		dj_mean = s->Master_dj;
	}

	else
	{
		if (BATCH_SIZE >= 30)
		{
			s->Obj_lb_U = s->Obj_lb_mean + 1.96 * s->Obj_lb_stdev;
			s->Obj_lb_L = s->Obj_lb_mean - 1.96 * s->Obj_lb_stdev;
		}
		else
		{
			/*Needs update here, 2.093 is for BATCH_SIZE = 20 */
			s->Obj_lb_U = s->Obj_lb_mean + 2.093 * s->Obj_lb_stdev;
			s->Obj_lb_L = s->Obj_lb_mean - 2.093 * s->Obj_lb_stdev;
		}
        if (status == 3) {
          /* In this case, 3 replications are finished and thus the multiplier is for t dist with 2 degrees of freedom */
          s->Obj_lb_U = s->Obj_lb_mean + 4.271 * s->Obj_lb_stdev;
          s->Obj_lb_L = s->Obj_lb_mean - 4.271 * s->Obj_lb_stdev;
        }

		//if (s->Obj_comp_mean >= s->Obj_lb_L && s->Obj_comp_mean <= s->Obj_lb_U) {
		if (status == 1)
		{
			soln_status = 1;
			mean = s->Obj_comp_mean;
			ul = s->Obj_comp_mean + 1.96 * s->Obj_comp_stdev;
			ll = s->Obj_comp_mean - 1.96 * s->Obj_comp_stdev;
			x = s->incumb_d;
			pi_mean = s->Batch_pi_mean;
			dj_mean = s->Batch_dj_mean;
			pi_stdev = s->Batch_pi_stdev;
			dj_stdev = s->Batch_dj_stdev;
		}
		// else if (s->Obj_mean_mean >= s->Obj_lb_L && s->Obj_mean_mean <= s->Obj_lb_U){
		else if (status == 2 || status == 3)
		{
			soln_status = 2;
			mean = s->Obj_mean_mean;
			ul = s->Obj_mean_mean + 1.96 * s->Obj_mean_stdev;
			ll = s->Obj_mean_mean - 1.96 * s->Obj_mean_stdev;
			x = s->incumb_avg;
			pi_mean = s->R_Master_pi_mean;
			dj_mean = s->R_Master_dj_mean;
			pi_stdev = s->R_Master_pi_stdev;
			dj_stdev = s->R_Master_dj_stdev;
		}
		else
		{
			soln_status = 0;
			return 1;
		}
	}

	TxX_plus(p->A, x, A_x);

	/* Print out solution status Compromise/Mean */
	if (status == 0)
	{
		fprintf(fp, "Replication No. %d\n", p->current_batch_id);
	}
	if (status == 1 || status == 2)
	{
		fprintf(fp, "%-40s%d\n", "Number of replications:", BATCH_SIZE);
	}
    if (status == 3) {
      fprintf(fp, "%-40s%d\n", "Number of replications:", 3);
    }
	fprintf(
			fp,
			"%-40s%s\n",
			"Status:",
			soln_status == 1 ? "COMPROMISE SOLUTION" :
			soln_status == 2 ? "MEAN SOLUTION" :
			soln_status == 3 ? "REPLICATION SOLUTION" :
			soln_status == 0 ? "NO  SOLUTION" : "???");
	/* Print out Objective Upper bound and its 95% CI's */
	if (sd_global->config.EVAL_FLAG == 0 && soln_status == 3)
	{
		fprintf(fp, "%-40s\n", "Total Objective Function Upper Bound:");
	}
	else
	{
		fprintf(fp, "%-40s%.2f,[%f,%f],%.2f\n",
				"Total Objective Function Upper Bound:", mean, ll, ul,
				(ul - mean) / 1.96);
	}

	/* Print out Objective Lower bound and its 95% CI's */
	if (status == 0)
	{
		fprintf(fp, "%-40s%f\n", "Total Objective Function Lower Bound:",
				s->incumb_est);
	}
	else
	{
		fprintf(fp, "%-40s%.2f,[%f,%f],%.2f\n",
				"Total Objective Function Lower Bound:", s->Obj_lb_mean,
				s->Obj_lb_L, s->Obj_lb_U, s->Obj_lb_stdev);
	}

	/* Print our row infomation */
	fprintf(fp, "\nFirst Stage Solutions:\n");
	fprintf(fp, "   No.   Row name   Activity      Lower bound  "
			" Upper bound   Dual          Dual STDEV\n");
	fprintf(fp, "------ ------------ ------------- ------------- "
			"------------- ------------- -------------\n");

	/* This is the Objective Row */
	fprintf(fp, "%6d ", 0);
	fprintf(fp, "%-12s ", p->master->objname == NULL ? "" : p->master->objname);
	fprintf(fp, "%13.6g ", fabs(mean) <= 1e-5 ? 0.0 : mean);
	fprintf(fp, "\n");

	for (i = 0; i < p->num->mast_rows; i++)
	{
		fprintf(fp, "%6d ", i + 1);
		if (p->master->rname[i] == NULL || strlen(p->master->rname[i]) <= 12)
			fprintf(fp, "%-12s ",
					p->master->rname[i] == NULL ? "" : p->master->rname[i]);
		else
			fprintf(fp, "%s\n%20s", p->master->rname[i], "");

		fprintf(fp, "%13.6g ", fabs(A_x[i + 1]) <= 1e-5 ? 0.0 : A_x[i + 1]);
		if (p->master->senx[i] == 'G' || p->master->senx[i] == 'E')
			fprintf(fp, "%13.6g ", p->master->rhsx[i]);
		else
			fprintf(fp, "%13s ", "");

		if (p->master->senx[i] == 'L' || p->master->senx[i] == 'E')
			fprintf(fp, "%13.6g ", p->master->rhsx[i]);
		else
			fprintf(fp, "%13s ", "");
		fprintf(fp, "%13.6g ",
				fabs(pi_mean[i + 1]) <= 1e-5 ? 0.0 : pi_mean[i + 1]);
		if (soln_status != 3)
		{
			fprintf(fp, "%13.6g ",
					pi_stdev[i + 1] <= 1e-5 ? 0.0 : pi_stdev[i + 1]);
		}
		fprintf(fp, "\n");
	}

	/* Print our column infomation */
	fprintf(fp, "\n");
	fprintf(fp, "   No. Column name  Activity      Lower bound  "
			" Upper bound   Reduced Cost  RC STDEV\n");
	fprintf(fp, "------ ------------ ------------- ------------- "
			"------------- ------------- -------------\n");

	for (i = 0; i < p->num->mast_cols; i++)
	{
		fprintf(fp, "%6d ", i + 1);
		if (p->master->cname[i] == NULL || strlen(p->master->cname[i]) <= 12)
			fprintf(fp, "%-12s ",
					p->master->cname[i] == NULL ? "" : p->master->cname[i]);
		else
			fprintf(fp, "%s\n%20s", p->master->cname[i], "");

		fprintf(fp, "%13.6g ", fabs(x[i + 1]) <= 1e-5 ? 0.0 : x[i + 1]);

		if (p->master->bdl[i] < -1e+20)
		{
			fprintf(fp, "%13s", "");
		}
		else
		{
			fprintf(fp, "%13.6g ", p->master->bdl[i]);
		}

		if (p->master->bdu[i] >= 1e+20)
		{
			fprintf(fp, "%13s", "");
		}
		else
		{
			fprintf(fp, "%13.6g ", p->master->bdu[i]);
		}

		fprintf(fp, "%13.6g ",
				fabs(dj_mean[i + 1]) <= 1e-5 ? 0.0 : dj_mean[i + 1]);
		if (soln_status != 3)
		{
			fprintf(fp, "%13.6g ",
					fabs(dj_stdev[i + 1]) <= 1e-5 ? 0.0 : dj_stdev[i + 1]);
		}

		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n");

	fclose(fp);
	mem_free(A_x);
	return 0;
}

int test_average_soln(sdglobal_type* sd_global, soln_type *s, prob_type *p)
{
	int i, j;
	double a, b;

	for (j = 1; j <= p->master->mac; j++)
	{
		a = 0.0;
		for (i = 0; i < BATCH_SIZE; i++)
		{
			a +=
					DBL_ABS((sd_global->batch_incumb->incumb_x[i][j] - s->incumb_avg[j]));
		}
		a /= 20;
		b = DBL_ABS(sd_global->config.MEAN_DEV * s->incumb_avg[j]);
		if (a > 0.001 && b > 0.001 && a > b)
		{
			return 0;
		}
	}
	return 1;
}

void record_master_dual_and_obj(sdglobal_type* sd_global, prob_type *prob,
		soln_type *soln, FILE *master_dual, FILE *master_obj)
{
	int i;

	master_dual = fopen("Master_Dual.out", "a");
	fprintf(master_dual, "Here are the duals for batch%d\n",
			prob->current_batch_id);
	sd_global->batch_incumb->R_Master_pi[prob->current_batch_id][0] =
			soln->Master_pi[0];
	for (i = 1; i <= prob->num->mast_rows + prob->num->max_cuts; i++)
	{
		sd_global->batch_incumb->R_Master_pi[prob->current_batch_id][i] =
				soln->Master_pi[i];
		fprintf(master_dual, "%f\n", soln->Master_pi[i]);
	}
	fprintf(master_dual, "Here are the dual slacks for batch%d\n",
			prob->current_batch_id);
	sd_global->batch_incumb->R_Master_dj[prob->current_batch_id][0] =
			soln->Master_dj[0];
	for (i = 1; i <= prob->num->mast_cols; i++)
	{
		sd_global->batch_incumb->R_Master_dj[prob->current_batch_id][i] =
				soln->Master_dj[i];
		fprintf(master_dual, "%f\n", soln->Master_dj[i]);
	}
	fclose(master_dual);

	master_obj = fopen("Master_Obj.out", "a");
	fprintf(master_obj, "Here is the optimal value\n%f\n", soln->opt_value);
	fclose(master_obj);
}

void process_batch_prob(sdglobal_type* sd_global, prob_type *prob,
		soln_type *soln, double *obj, FILE *batch_dual, FILE *batch_d,
		FILE *batch_obj)
{
	int i, j, k;
	write_prob(sd_global->batch_problem, "final-batch-prob.lp");
	solve_problem(sd_global, sd_global->batch_problem);
	*obj = get_objective(sd_global->batch_problem);
	for (prob->num->batch_id = 0; prob->num->batch_id < BATCH_SIZE;
			prob->num->batch_id++)
	{
		get_dual(soln->Batch_pi[prob->num->batch_id], sd_global->batch_problem,
				prob->num, prob->num->mast_rows);
		get_dual_slacks(soln->Batch_dj[prob->num->batch_id],
				sd_global->batch_problem, prob->num, prob->num->mast_cols);

		for (k = 0; k <= prob->num->mast_rows + prob->num->max_cuts; k++)
		{
			soln->Batch_pi[prob->num->batch_id][k] = BATCH_SIZE
					* soln->Batch_pi[prob->num->batch_id][k];
		}
		for (k = 0; k <= prob->num->mast_cols; k++)
		{
			soln->Batch_dj[prob->num->batch_id][k] = BATCH_SIZE
					* soln->Batch_dj[prob->num->batch_id][k];
		}

		batch_dual = fopen("Batch_dual.out", "a");
		fprintf(batch_dual, "Here are the duals for batch%d\n",
				prob->num->batch_id);
		for (i = 1; i <= prob->num->mast_rows + prob->num->max_cuts; i++)
		{
			fprintf(batch_dual, "%f\n", soln->Batch_pi[prob->num->batch_id][i]);
		}
		fprintf(batch_dual, "Here are the dual slacks for batch%d\n",
				prob->num->batch_id);
		for (i = 1; i <= prob->num->mast_cols; i++)
		{
			fprintf(batch_dual, "%f\n", soln->Batch_dj[prob->num->batch_id][i]);
		}
		fclose(batch_dual);
	}

	/* modified by Yifan 2012.10.05 */
	/* Calculate mean and stdard deviation of duals and reduced costs for compromise solution*/
	calc_mean_stdev(soln->Batch_pi, soln->Batch_pi_mean, soln->Batch_pi_stdev,
			prob->num->mast_rows, BATCH_SIZE);
	calc_mean_stdev(soln->Batch_dj, soln->Batch_dj_mean, soln->Batch_dj_stdev,
			prob->num->mast_cols, BATCH_SIZE);

	/* modified by Yifan 2012.10.05 */
	/* Calculate mean and stdard deviation of duals and reduced costs for mean solution*/
	calc_mean_stdev(sd_global->batch_incumb->R_Master_pi,
			soln->R_Master_pi_mean, soln->R_Master_pi_stdev,
			prob->num->mast_rows, BATCH_SIZE);
	calc_mean_stdev(sd_global->batch_incumb->R_Master_dj,
			soln->R_Master_dj_mean, soln->R_Master_dj_stdev,
			prob->num->mast_cols, BATCH_SIZE);

	calc_var(sd_global, sd_global->Obj_lb, &soln->Obj_lb_mean,
			&soln->Obj_lb_stdev, BATCH_SIZE);

	get_primal(soln->incumb_d, sd_global->batch_problem, prob->num->mast_cols);
	/* Get the unified X by adding the first incumbent solution x1 and d1 */
	batch_d = fopen("Batch_x.out", "a");
	fprintf(batch_d, "Here are the batch X:\n");
	for (i = 1; i <= prob->num->mast_cols; i++)
	{
		soln->incumb_d[i] += sd_global->batch_incumb->incumb_x[0][i];
		fprintf(batch_d, "%f\n", soln->incumb_d[i]);
	}
	/* c_xc is used to estimate the lower bound of the problem */
	soln->c_xc = CxX(prob->c, soln->incumb_d, prob->num->mast_cols);
	get_beta_x(sd_global, soln, soln->beta, sd_global->batch_problem, prob->num,
			prob->num->mast_cols);
	fprintf(batch_d, "Here are the Average X of replications:\n");
	for (i = 1; i <= prob->num->mast_cols; i++)
	{
		soln->incumb_avg[i] = 0.0;
		for (j = 0; j < BATCH_SIZE; j++)
		{
			soln->incumb_avg[i] += sd_global->batch_incumb->incumb_x[j][i];
		}
		soln->incumb_avg[i] /= (1.0 * BATCH_SIZE);
	}
	for (i = 1; i <= prob->num->mast_cols; i++)
	{
		fprintf(batch_d, "%f\n", soln->incumb_avg[i]);
	}
	fclose(batch_d);

	batch_obj = fopen("Batch_Obj.out", "a");
	fprintf(batch_obj, "Here is the optimal value\n%f\n", soln->opt_value);
	fclose(batch_obj);

}

/* modified by Yifan 2013.02.15 */
void calc_lowerbound(sdglobal_type* sd_global, prob_type *prob, soln_type *soln)
{
	double b_lambda[BATCH_SIZE];
	double alpha_theta[BATCH_SIZE];
//	double L_Mu[BATCH_SIZE];
	double U_Mu[BATCH_SIZE];
	double l[BATCH_SIZE];
	double lrep[BATCH_SIZE][BATCH_SIZE];
	double lavg[BATCH_SIZE];
	//double *rhs1;
	int idx;
	int cnt;
	int rep;
	int row_num;
	FILE *clb;

	/*
	 if (!(rhs1 = arr_alloc(prob->num->mast_rows+1, double)))
	 err_msg ("Allocation", "calc_lowerbound", "rhs1");

	 for (cnt=0; cnt<prob->num->mast_rows; cnt++) {
	 rhs1[cnt+1] = prob->master->rhsx[cnt];
	 }

	 TxX(prob->A, sd_global->batch_incumb->incumb_x[0], rhs1);*/
	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		b_lambda[cnt] = 0.0;
		alpha_theta[cnt] = 0.0;
//		L_Mu[cnt] = 0.0;
		U_Mu[cnt] = 0.0;
		l[cnt] = 0.0;
	}

	/* 1. Lower bound from compromise problem */
	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		for (idx = 0; idx < prob->num->mast_rows; idx++)
		{
			/*Be careful here, rhsx starts from array location 0,
			 while R_Master_pi starts from array location 1  2013/01/22 Yifan*/
			// b_lambda += rhs1[idx] * sd_global->batch_incumb->R_Master_pi[0][idx+1];
			b_lambda[cnt] += prob->master->rhsx[idx]
					* soln->Batch_pi[cnt][idx + 1];
		}
		for (idx = 0; idx < sd_global->bcuts->batch[cnt]->cnt; idx++)
		{
			row_num = sd_global->bcuts->batch[cnt]->val[idx]->row_num;
			alpha_theta[cnt] +=
					(sd_global->bcuts->batch[cnt]->val[idx]->alpha_incumb
							+ CxX(sd_global->bcuts->batch[cnt]->val[idx]->beta,
									sd_global->batch_incumb->incumb_x[cnt],
									prob->num->mast_cols))
							* soln->Batch_pi[cnt][row_num + 1];
		}

		for (idx = 0; idx < prob->num->mast_cols; idx++)
		{
			U_Mu[cnt] += sd_global->batch_incumb->incumb_x[cnt][idx + 1]
					* soln->Batch_dj[cnt][idx + 1];
		}
		l[cnt] = b_lambda[cnt] + alpha_theta[cnt] + U_Mu[cnt];
	}

	/* 2. Lower bound from average dual from replications */
	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		lavg[cnt] = 0.0;
		b_lambda[cnt] = 0.0;
		alpha_theta[cnt] = 0.0;
//		L_Mu[cnt] = 0.0;
		U_Mu[cnt] = 0.0;
	}
	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{

		for (rep = 0; rep < BATCH_SIZE; rep++)
		{
			for (idx = 0; idx < prob->num->mast_rows; idx++)
			{
				b_lambda[cnt] += prob->master->rhsx[idx]
						* sd_global->batch_incumb->R_Master_pi[rep][idx + 1];
			}
			for (idx = 0; idx < sd_global->bcuts->batch[cnt]->cnt; idx++)
			{
				row_num = sd_global->bcuts->batch[cnt]->val[idx]->row_num;
				alpha_theta[cnt] +=
						(sd_global->bcuts->batch[cnt]->val[idx]->alpha_incumb
								+ CxX(
										sd_global->bcuts->batch[cnt]->val[idx]->beta,
										sd_global->batch_incumb->incumb_x[cnt],
										prob->num->mast_cols))
								* sd_global->batch_incumb->R_Master_pi[rep][row_num
										+ 1];
			}

			for (idx = 0; idx < prob->num->mast_cols; idx++)
			{
				U_Mu[cnt] += sd_global->batch_incumb->incumb_x[cnt][idx + 1]
						* sd_global->batch_incumb->R_Master_dj[rep][idx + 1];
			}
		}
		lavg[cnt] = (b_lambda[cnt] + alpha_theta[cnt] + U_Mu[cnt]) / BATCH_SIZE;
	}

	/* 3. Lower bound from each replication */
	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		for (rep = 0; rep < BATCH_SIZE; rep++)
		{
			lrep[rep][cnt] = 0.0;
		}
	}

	for (rep = 0; rep < BATCH_SIZE; rep++)
	{
		/* At the begining of each replication, clean up the array  */
		for (cnt = 0; cnt < BATCH_SIZE; cnt++)
		{
			b_lambda[cnt] = 0.0;
			alpha_theta[cnt] = 0.0;
//			L_Mu[cnt] = 0.0;
			U_Mu[cnt] = 0.0;
		}
		for (cnt = 0; cnt < BATCH_SIZE; cnt++)
		{
			for (idx = 0; idx < prob->num->mast_rows; idx++)
			{
				b_lambda[cnt] += prob->master->rhsx[idx]
						* sd_global->batch_incumb->R_Master_pi[rep][idx + 1];
			}
			for (idx = 0; idx < sd_global->bcuts->batch[cnt]->cnt; idx++)
			{
				row_num = sd_global->bcuts->batch[cnt]->val[idx]->row_num;
				alpha_theta[cnt] +=
						(sd_global->bcuts->batch[cnt]->val[idx]->alpha_incumb
								+ CxX(
										sd_global->bcuts->batch[cnt]->val[idx]->beta,
										sd_global->batch_incumb->incumb_x[cnt],
										prob->num->mast_cols))
								* sd_global->batch_incumb->R_Master_pi[rep][row_num
										+ 1];
			}

			for (idx = 0; idx < prob->num->mast_cols; idx++)
			{
				U_Mu[cnt] += sd_global->batch_incumb->incumb_x[cnt][idx + 1]
						* sd_global->batch_incumb->R_Master_dj[rep][idx + 1];
			}
			lrep[rep][cnt] = b_lambda[cnt] + alpha_theta[cnt] + U_Mu[cnt];
		}
	}

	if (!(clb = fopen("complb.txt", "a")))
	{
		printf("WARNING :: eval.dat file not opened\n");
		return;
	}

	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		fprintf(clb, "%f\t", l[cnt]);
	}
	fprintf(clb, "\n\n");

	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		fprintf(clb, "%f\t", lavg[cnt]);
	}

	for (rep = 0; rep < BATCH_SIZE; rep++)
	{
		fprintf(clb, "\n\n");
		for (cnt = 0; cnt < BATCH_SIZE; cnt++)
		{
			fprintf(clb, "%f\t", lrep[rep][cnt]);
		}
	}

	fclose(clb);
}

/************************************************************************\
 ** This function takes the SD code into Feasibility mode (solve_cell() take
 the SD into Optimality mode). The SD will not return to optimality mode
 until the candidate and incumbent solution are both feasibile.
 \************************************************************************/
void resolve_infeasibility(sdglobal_type* sd_global, cell_type *c, prob_type *p,
		soln_type *s, int omeg_idx, BOOL new_omega, double *phi)
{
	/* QP master will be solved in Feasibility mode instead of Optimality mode */
	c->opt_mode = FALSE;

	if (new_omega == TRUE)
	{
		/* Update optimality cut data structure if new omega is observered. */
		//?????calc_delta_col(s->delta, c->lambda, s->omega, p->num, omeg_idx);
		form_fea_cut(sd_global, p, c, s, omeg_idx, s->candid_x, new_omega);
	}

	if (c->incumb_infea == TRUE)
	{
		/*why this line of code will cause problem??????*/
		replace_incumbent(sd_global, p, c, s, phi);
		printf("Note: Incumbent replaced by candidate! ri-1\n");
	}

	while (1)
	{
		/* Solve QP master */
		/* Yifan 03/26/2012 Test for New Ideas on quad_scalar */
		c->quad_scalar = sd_global->config.MIN_QUAD_SCALAR;
		construct_QP(p, c, c->quad_scalar);
		solve_QP_master(sd_global, p, c, s); /* QP solves in feasibility mode */

		c->fea_count++;
		/* Solve Sub problem with candidate solution */
		solve_subprob(sd_global, p, c, s, s->candid_x, omeg_idx);

		if (p->subprob->feaflag == TRUE)
		{
			break;
		}

		form_fea_cut(sd_global, p, c, s, omeg_idx, s->candid_x, FALSE);

		if (c->incumb_infea == TRUE)
		{
			replace_incumbent(sd_global, p, c, s, phi);
			printf("Note: Incumbent replaced by candidate! ri-2\n");
		}
		else
		{
			printf("Incumbent feasibile but candidate infeasible!\n");
		}

	}

	form_fea_cut(sd_global, p, c, s, omeg_idx, s->candid_x, FALSE);

	if (c->incumb_infea == TRUE)
	{
		replace_incumbent(sd_global, p, c, s, phi);
		/*added by Yifan 04/22/2013 Immediately solve master after incumbent change */

		printf("Note: Incumbent replaced by candidate! ri-3\n");
	}

	s->incumbent_change = TRUE;

	c->opt_mode = TRUE;

	printf("fea_count is : %d\n", c->fea_count);

}

void replace_incumbent(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, double *phi)
{
	/*
	 double candid_est;

	 if (c->cuts->cnt > 0) {
	 candid_est = max_cut_height(c->cuts, s->candid_x, c, p->num);
	 candid_est += CxX(p->c, s->candid_x, p->num->mast_cols);

	 s->incumb_est = max_cut_height(c->cuts, s->incumb_x, c, p->num);
	 s->incumb_est += CxX(p->c, s->incumb_x, p->num->mast_cols);
	 }
	 else {
	 candid_est = s->candid_est;
	 }*/

	/* Yifan 03/25/2012 Test for New Ideas on quad_scalar */
	/*calc_quad_scalar(p, c, s, candid_est);*/

	new_incumbent(p, c, s, s->candid_est);

	/* Yifan 03/25/2012 Test for New Ideas on quad_scalar */

	if (c->k > 1 && s->norm_d_k > sd_global->config.TOLERANCE)
		if (s->norm_d_k >= sd_global->config.R3 * s->norm_d_k_1)
		{

			c->quad_scalar *= sd_global->config.R2 * sd_global->config.R3
					* s->norm_d_k_1 / s->norm_d_k;

			c->quad_scalar =
					min(sd_global->config.MAX_QUAD_SCALAR, c->quad_scalar );

			c->quad_scalar =
					max(sd_global->config.MIN_QUAD_SCALAR, c->quad_scalar);

			if (sd_global->config.MASTER_TYPE == SDQP)
				construct_QP(p, c, c->quad_scalar);
		}

	if (sd_global->config.MASTER_TYPE == SDQP)
	{
		change_rhs(p, c, s);
		change_bounds(p, c, s);
	}

	phi[0] = s->incumb_est; /*  SS messing around  */
	phi[1] = sd_global->config.SMOOTH_PARM * s->incumb_est
			+ (1 - sd_global->config.SMOOTH_PARM) * phi[1]; /*  SS messing around  */
	phi[2] = 0.0;
	phi[3] = (double) c->k;
	/*printf("*****   NEW INCUMBENT FROM CANDIDATE  *****\n");*/
	printf("+");
	fflush(stdout);

	/* Added to assure the validity of zero as subproblem LB.
	 * Lei 09/13/05 */
	if (sd_global->config.SUB_LB_CHECK)
		printf("Subproblem LB = %lf\n", s->sub_lb_checker);

	s->norm_d_k_1 = s->norm_d_k; /* zl. 06/10/02 */

	/* Since it is now replaced by a candidate, we assume it is feasibile now*/
	c->incumb_infea = FALSE;
	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	s->gamma = 0.0;
}

void calc_quad_scalar(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, double candid_est)
{
	double *q_vector, *d; /* vector: c + Bk_theta - A_Trans_lambda. */
//	double q_term;
	double dir_derivative = 0.0;
	double gamma;
	int cnt;

	if (!(q_vector = arr_alloc(p->num->mast_cols+2, double)))
		err_msg("Allocation", "cal_temp_lb", "q_vector");
	if (!(d = arr_alloc(p->num->mast_cols+2, double)))
		err_msg("Allocation", "cal_temp_lb", "d");

//	q_term = calc_q_vec(p, c, s, q_vector);

	for (cnt = 1; cnt <= p->num->mast_cols; cnt++)
	{
		d[cnt] = s->incumb_x[cnt] - s->candid_x[cnt];
		dir_derivative += q_vector[cnt] * d[cnt];
	}

	/* Both gamma and dir_derivative are negative */

	if (dir_derivative < 0)
	{
		gamma = candid_est - s->incumb_est;
		c->quad_scalar =
				max(sd_global->config.MIN_QUAD_SCALAR, dir_derivative/gamma);

		if (sd_global->config.MASTER_TYPE == SDQP)
			construct_QP(p, c, c->quad_scalar);
	}

}

double calc_q_vec(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, double *q_vec)
{
	double *bk; /* vector: b - A*incumb_x. */
	double *lambda; /* vector: the dual of the primal constraints. */
	double bk_lambda; /* scalar: bk*lambda. */

	sparse_matrix *A_Trans; /* sparse_matrix: the transpose of A. */
	double *A_Trans_lambda; /* vector: - A_Trans * lambda. */

	double *theta; /* the dual of the reformed cut constraints. */
	double *Vk; /* the vector of scalars of cut constraints.
	 Vk = alpha - beta * incumb_x  */
//	double Vk_theta; /* scalar: Vk*theta. */

	double *Bk_theta; /* vector: Bk_Transpose * theta, where Bk_Transpose is
	 the matrix of cut coefficients. */
	double *Bk_col; /* vector: Store one column of Bk while calculating
	 Bk_theta. */
	double q_term; /* scalar: q_vec * q_vec. */

	double eta_zero;
	int cnt;
	int i;

#ifdef TRACE
	printf ("Inside calc_q_vec.\n");
#endif

	if (!(bk =
			arr_alloc(p->num->mast_rows + c->feasible_cuts_added->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "bk");
	if (!(lambda =
			arr_alloc(p->num->mast_rows + c->feasible_cuts_added->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "lambda");
	/* Yifan 03/19/2012 Test update A and B matrix*/
	if (!(A_Trans_lambda = arr_alloc(p->num->mast_cols+2, double)))
		err_msg("Allocation", "cal_temp_lb", "A_lambda");
	if (!(A_Trans = (sparse_matrix *) mem_malloc(sizeof(sparse_matrix))))
		err_msg("Allocation", "cal_temp_lb", "A_Trans");
	if (!(A_Trans->val = arr_alloc(p->A->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "A_Trans->val");
	if (!(A_Trans->row = arr_alloc(p->A->cnt+1, int)))
		err_msg("Allocation", "cal_temp_lb", "A_Trans->row");
	if (!(A_Trans->col = arr_alloc(p->A->cnt+1, int)))
		err_msg("Allocation", "cal_temp_lb", "A_Trans->col");
	if (!(theta = arr_alloc(c->cuts->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "theta");
	if (!(Vk = arr_alloc(c->cuts->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "Vk");
	/* Yifan 03/19/2012 Test update A and B matrix*/
	if (!(Bk_theta = arr_alloc(p->num->mast_cols+2, double)))
		err_msg("Allocation", "cal_temp_lb", "Bk_theta");
	if (!(Bk_col = arr_alloc(c->cuts->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "Bk_col");

	eta_zero = (sd_global->Abar
			- CxX(sd_global->Bbar, s->incumb_x, p->num->mast_cols));

	/* 1. bk */
	for (cnt = 0; cnt < p->num->mast_rows - 1; cnt++)
	{
		bk[cnt + 1] = p->master->rhsx[cnt];
	}

	TxX(p->A, s->incumb_x, bk);

	/* Yifan 03/18/2012 update rhs of Lower Bound constraint */
	if (sd_global->config.LB_TYPE == 0)
	{
		bk[p->num->mast_rows] = 0.0;
	}
	else
		bk[p->num->mast_rows] = sd_global->Eta0;

	/* Yifan 03/11/2012 Updated for Feasibility Cuts*/
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		bk[p->num->mast_rows + cnt + 1] =
				c->feasible_cuts_added->val[cnt]->alpha
						- CxX(c->feasible_cuts_added->val[cnt]->beta,
								s->incumb_x, p->num->mast_cols);
	}

	/* 2. lambda*/

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	/* Obtain lambda from s->Master_pi of original master problem constraints. */
	/* Dual values' sign need to be flipped here before assigning to lambda*/
	for (cnt = 0; cnt < p->num->mast_rows - 1; cnt++)
	{
		lambda[cnt + 1] = -s->Master_pi[cnt + 1];
	}

	/* CPLEX return a positive dual multipliers for this eta0 equality constraint*/
	/* Yifan 03/19/2012 Test removing eta0*/
	lambda[p->num->mast_rows] = s->Master_pi[p->num->mast_rows];

	/* Yifan 03/11/2012 Updated for Feasibility Cuts NOTICEHERE*/
	/* Obtain lambda from s->Master_pi of  master problem feasibility cuts. */
	/* Dual values' sign need to be flipped here before assigning to lambda*/
	/* Yifan 03/19/2012 Test removing eta0*/
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		lambda[p->num->mast_rows + cnt + 1] =
				-s->Master_pi[c->feasible_cuts_added->val[cnt]->row_num + 1];
	}

	/* 3. bk_lambda = bk * lambda.*/
	bk_lambda = CxX(bk, lambda,
			p->num->mast_rows - 1 + c->feasible_cuts_added->cnt);

	/* Yifan 03/14/2012 Make sure addition or subtraction*/
	/* Take care of the bound constraints */
	for (i = 1; i <= p->num->mast_cols + 1; i++)
		bk_lambda += s->Master_dj[i] * s->incumb_x[i];

	/* 4. A_Trans */
	A_Trans->cnt = p->A->cnt;

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	for (cnt = 1; cnt <= A_Trans->cnt; cnt++)
	{
		A_Trans->val[cnt] = p->A->val[cnt];
		A_Trans->row[cnt] = p->A->col[cnt];
		A_Trans->col[cnt] = p->A->row[cnt];
	}

	/* 5. A_Trans_lambda */

	/* Calculate - A_Trans * lambda due to _TxX()_ function */
	TxX(A_Trans, lambda, A_Trans_lambda);

	/* Yifan 03/11/2012 Updated for Lower Bound constraint*/
	/* Yifan 03/19/2012 Test update A_trans_lambda vector, eta0 equality constraint*/
	for (i = 0; i < p->num->mast_cols; i++)
	{
		A_Trans_lambda[i + 1] -= sd_global->Bbar[i + 1]
				* lambda[p->num->mast_rows];
	}

	/* Yifan 03/11/2012 Updated for Feasibility Cuts*/
	for (i = 0; i < p->num->mast_cols; i++)
	{
		for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
		{
			A_Trans_lambda[i + 1] -= c->feasible_cuts_added->val[cnt]->beta[i
					+ 1] * lambda[p->num->mast_rows + cnt + 1];
		}
	}

	/* Yifan 03/19/2012 A_trans_lambda vector, eta0 equality constraint*/
	A_Trans_lambda[p->num->mast_cols + 1] = -lambda[p->num->mast_rows];

	/* Yifan 03/14/2012 Make sure addition or subtraction*/
	/* These are for bound constraints */
	for (i = 0; i < p->num->mast_cols + 1; i++)
	{
		A_Trans_lambda[i + 1] -= s->Master_dj[i + 1];
	}

	/* 6. theta and 7. Vk */
	for (cnt = 0; cnt < c->cuts->cnt; cnt++)
	{
		/* 6. theta*/
		/* Yifan 03/14/2012 Duals for optimality cuts*/
		theta[cnt + 1] = ((double) c->k / (double) c->cuts->val[cnt]->cut_obs)
				* s->Master_pi[c->cuts->val[cnt]->row_num + 1];

		/* 7. Vk */
		/* Yifan 03/14/2012 Updated for optimality cut height*/
		Vk[cnt + 1] = c->cuts->val[cnt]->alpha
				- CxX(c->cuts->val[cnt]->beta, s->incumb_x, p->num->mast_cols);
		Vk[cnt + 1] *= (double) c->cuts->val[cnt]->cut_obs / (double) c->k;
		Vk[cnt + 1] += (1 - (double) c->cuts->val[cnt]->cut_obs / (double) c->k)
				* eta_zero;
	}

	/* 8. Vk_theta. */
//	Vk_theta = CxX(Vk, theta, c->cuts->cnt);
	/* 9. Bk and 10. Bk_theta. */
	/*
	 ** Calculate Bk_theta. Bk is the transpose of the matrix of cut
	 ** coefficients {T->beta}, where T are the reformed cuts.
	 */
	for (i = 1; i <= p->num->mast_cols; i++)
	{

		/* 9. Bk. */
		for (cnt = 0; cnt < c->cuts->cnt; cnt++)
		{
			/* Yifan 03/14/2012 might need updates here!!!*/
			Bk_col[cnt + 1] = c->cuts->val[cnt]->beta[i];
			Bk_col[cnt + 1] *= (double) c->cuts->val[cnt]->cut_obs
					/ (double) c->k;
		}

		/* 10. Bk_theta. */
		Bk_theta[i] = CxX(Bk_col, theta, c->cuts->cnt);
	}

	/* Yifan 03/19/2012 Update Bk_theta vector's eta0 column*/
	i = p->num->mast_cols + 1;
	/* Since the coef of eta0 column in opt cuts are zeros*/
	Bk_theta[i] = 0.0;

	/* Calculate the quadratic vector, q_vec.
	 Note: A_Trans_lambda = - A_Trans * lambda. */
	/* 11. q_vec */
	for (i = 1; i <= p->num->mast_cols; i++)
	{

		q_vec[i] = p->c[i] - Bk_theta[i] - A_Trans_lambda[i];

	}

	/* Yifan 03/19/2012 Update q_vec vector*/
	i = p->num->mast_cols + 1;
	q_vec[i] = -Bk_theta[i] - A_Trans_lambda[i];
	q_term = CxX(q_vec, q_vec, p->num->mast_cols + 1);

	return q_term;

}

/************************************************************************\
** This function represents the SD algorithm, as solved for a
 ** single cell.  It creates temporary data structures required for
 ** the solution, then proceeds with the iterative algorithm.
 ** It terminates when an optimal solution has been discovered,
 ** or when the limit on the number of iterations has been reached.
 ** It assumes that the cell passed to it has already been allocated
 ** and initialized, and likewise that the problem has been initialized.
 ** status = 0 for replications
 ** status = 1 for compromise evaluation
 ** status = 2 for mean solution evaluation
 ** status = 3 for early termintated mean solution evaluation
 \************************************************************************/
void evaluate_inc(sdglobal_type* sd_global, cell_type *cell, prob_type *prob,
		soln_type *s, vector x_k, char *fname, double *conf_int, int status)
{
	soln_type *soln;
	int omeg_idx;
	/* int  i; */
    BOOL new_omega;
	FILE *fout /* , *fin, *fx */;
    FILE *time_sample;  /* modified by Yifan 2013.06.30 */
	double ans;
	double mean, stdev;
	double vari, temp;
	/* int		iteration; */
	int count;
	double cx;
	clock_t eval_start_time;
	clock_t eval_end_time;
	double total_time;


#ifdef DEBUG
	int idx;
#endif

#ifdef TRACE
	printf("Inside evaluate_inc\n");
#endif

	 /*x_k[1] = 2;
	 x_k[2] = 6;
	 x_k[3] = 5;
	 x_k[4] = 4;
	 x_k[5] = 193.47;
	 x_k[6] = 180.06;
	 x_k[7] = 185.91;
	 x_k[0] = x_k[1] + x_k[2] + x_k[3] + x_k[4] + x_k[5] + x_k[6] + x_k[7];*/
	eval_start_time = clock(); /*Yifan 2012-09-17*/

	soln = new_soln(sd_global, prob, x_k);
	if (status == 1)
	{
		printf("\n\n Begin evaluation of compromise solution \n");
	}
	else if (status == 2 || status == 3)
	{
		printf("\n\n Begin evaluation of average solution \n");
	}
	else
	{
		printf("\n\n Begin evaluation of incumbent solution \n");
	}
	count = 0;
	mean = 0.0;
	vari = 0.0;
	stdev = 10000000.0;
	temp = 0.0;

	cx = CxX(prob->c, soln->candid_x, prob->num->mast_cols);
    
	change_solver_primal(cell->subprob);
#ifdef OMEGA_FILE
  soln->omega->fidx[0] = s->omega->fidx[cell->k];
#endif
	while (3.92 * stdev > sd_global->config.EVAL_ERROR * DBL_ABS(mean)
			|| count < 100)
	{
		/* comment out by Yifan
		 if (!(count % 100))
		 printf("..%d.", count);
		 */
      
#ifdef OMEGA_FILE
      get_observ(sd_global, soln->omega, prob->num, &new_omega);        /* Yifan 2012.05.21 */
#else
      omeg_idx = generate_observ(sd_global, soln->omega, prob->num, &new_omega, &(sd_global->config.EVAL_SEED1));
#endif
      
		
		if (!solve_subprob(sd_global, prob, cell, soln, soln->candid_x, 0))
		{
			cplex_err_msg(sd_global, "Subproblem", prob, cell, soln);
			break;
		}
      
#ifdef REC_OMEGA
      int i;
      /* Yifan 05/10/2012 Scenario Reader */
      FILE *omg;
      omg = fopen("sampleOmegas_eval", "a");
      for(i=1; i <= sd_global->omegas.num_omega; i++){
        /* Note that RT start from 1 but mean start from 0*/
        fprintf(omg, "%f\t", soln->omega->RT[i]+sd_global->omegas.mean[i-1]);
      }
      fprintf(omg, "\n");
      fclose(omg);
#endif
      
      /* Yifan 2012.05.21 */
      /* Push the pointer to the next line of the omega */
#ifdef OMEGA_FILE
      if (soln->omega->fidx[0] == soln->omega->fidx[1]) {
        err_msg("Need more samples", "evaluate.c", "evaluate_inc()");   /* Generate more samples ... */
      }
      soln->omega->fidx[0] = soln->omega->fidx[1];
#endif
      /* Yifan 2012.05.21 */

		ans = get_objective(cell->subprob);

		if (count == 0)
		{
			mean = ans;
		}
		else
		{
			temp = mean;
			mean = mean + (ans - mean) / (double) (count + 1);
			vari = (1 - 1 / (double) count) * vari
					+ (count + 1) * (mean - temp) * (mean - temp);
			stdev = sqrt(vari / (double) count);
		}

		count++;

		/* Print the results every once in a while for long runs */
		if (!(count % 250))
		{
			printf(".");
			fflush(stdout);
		}
		if (!(count % 10000))
			printf("\n\nobs:%d mean:%lf   error: %lf \n 0.95 CI: [%lf , %lf]\n",
					count, cx + mean, 3.92 * stdev / mean,
					cx + mean - 1.96 * stdev, cx + mean + 1.96 * stdev);

		/* Squash the omega structure back down to nothing */
		soln->omega->cnt = 0;
		soln->omega->last = 0;
		soln->omega->most = 0;
		soln->omega->next = 0;
		soln->omega->filter[soln->omega->next] = UNUSED;
		soln->omega->weight[soln->omega->next] = 0;
		mem_free(soln->omega->idx[0]);
	}

	printf("\n");
	mean += cx;
	conf_int[0] = mean - 1.96 * stdev;
	conf_int[1] = mean + 1.96 * stdev;

	if (!(fout = fopen("eval.dat", "a")))
	{
		printf("WARNING :: eval.dat file not opened\n");
		return;
	}
    
	/*if (!(fx = fopen("solutions.dat", "a")))
	 { printf("WARNING :: solutions.dat file not opened\n"); return; }*/

	eval_end_time = clock();
	total_time = (double) (eval_end_time - eval_start_time);
	/* Print the value of the solution to a file and the screen */
	printf("\nFinal Estimate");
	printf("\n obs:%d mean:%lf \n 0.95 C.I.: [%lf , %lf] \n\n", count, mean,
			mean - 1.96 * stdev, mean + 1.96 * stdev);
	//fprintf(fout, "%d, %lf, %lf;\n", count, mean, stdev);
	fprintf(fout, "%d\t%lf\t%lf\t[%lf, %lf]\t%lf\ttime:%lf\t%d\n", cell->k,
			mean, stdev, mean - 1.96 * stdev, mean + 1.96 * stdev,
			s->incumb_est, total_time / CLOCKS_PER_SEC, count);
    
    /* modified by Yifan 2013.06.30 */
    if (!(time_sample = fopen("time_sample.out", "a"))) {
        printf("WARNING :: time_sample.out file not opened\n");
        return;
    }
    fprintf(time_sample, "%lf\t%d\n",total_time / CLOCKS_PER_SEC, count);
    fclose(time_sample);

	if (status == 1)
	{
		s->Obj_comp_mean = mean;
		s->Obj_comp_stdev = stdev;
	}
	else if (status == 2 || status == 3)
	{
		s->Obj_mean_mean = mean;
		s->Obj_mean_stdev = stdev;
	}
	else
	{
		s->rep_mean = mean;
		s->rep_stdev = stdev;
	}

	fclose(fout);
	free_soln(prob, cell, soln);
#ifdef TRACE
	printf("Exiting evaluate_inc\n");
#endif
	//fclose(fx);

	//fclose(fin);

}

