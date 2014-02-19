/****************************************************************************\
   This file includes all the routines necessary for solving the master
 problem using the regularized QP method. 
 \****************************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "log.h"
#include "sdglobal.h"

/***************************************************************************\
  Change the problem type to quadratic for regularized method.
 \***************************************************************************/
void change_prob_type(one_problem *m)
{
	int status = 0;

#ifdef TRACE
	printf("Inside change_prob_type.\n");
#endif

	status = change_probtype(m, PROB_QP); /* 2011.10.30 */

	if (status)
	{
		fprintf(stderr, "Failed to change master problem type to QP.\n");
		exit(1);
	}
}

/***************************************************************************\
  Consturct the Q diagonal matrix and copy it for quadratic problem. 
 \***************************************************************************/
void construct_QP(prob_type *p, cell_type *c, double sigma)
{
	int status = 0;
	int idx;
	double *qsepvec;

#ifdef TRACE
	printf("Inside construct_QP.\n");
#endif
	/* Yifan 03/21/2012 modified for lower bound*/
	if (!(qsepvec = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "construct_QP", "qsepvec");

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (idx = 0; idx < p->num->mast_cols; idx++)
	{
		qsepvec[idx] = 0.5 * sigma;
		/*  qsepvec[idx] = sigma; */
	}

	/* This is for eta column */
	qsepvec[p->num->mast_cols] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	/* 2011.10.30 - call copy_qp_separable() */
	status = copy_qp_separable(c->master, qsepvec); //modified by Yifan 09/14/2011
	mem_free(qsepvec);
	//added by Yifan to free memory 09/02/11

	if (status)
	{
		fprintf(stderr, "Failed to copy Q matrix.\n");
		exit(1);
	}
#ifdef CAL_CHECK
	write_prob(c->master, "quad.lp"); //modified by Yifan 09/14/2011
#endif
}

/***************************************************************************\
 Desturct the Q diagonal matrix and copy it for quadratic problem.
 \***************************************************************************/
void destruct_QP(prob_type *p, cell_type *c)
{
	int status = 0;
	int idx;
	double *qsepvec;
    
#ifdef TRACE
	printf("Inside destruct_QP.\n");
#endif
	/* Yifan 03/21/2012 modified for lower bound*/
	if (!(qsepvec = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "construct_QP", "qsepvec");
    
	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (idx = 0; idx < p->num->mast_cols; idx++)
	{
		qsepvec[idx] = 0;
	}
    
	/* This is for eta column */
	qsepvec[p->num->mast_cols] = 0.0;
    
	/* Now copy the empty Q matrix for QP problem. */
    
	status = copy_qp_separable(c->master, qsepvec); //modified by Yifan 09/14/2011
	mem_free(qsepvec);
	//added by Yifan to free memory 09/02/11
    
	if (status)
	{
		fprintf(stderr, "Failed to copy Q matrix.\n");
		exit(1);
	}
#ifdef CAL_CHECK
	write_prob(c->master, "quad_after_destruction.lp"); //modified by Yifan 09/14/2011
#endif
}

/* added by Yifan to construct batch mean QP problem 2012-09-09*/
void construct_batch_QP(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		double sigma)
{
	int status = 0;
	int idx;
	double *qsepvec;
	double avg_sigma = 0.0;

#ifdef TRACE
	printf("Inside construct_QP.\n");
#endif
	/* Yifan 03/21/2012 modified for lower bound*/
	if (!(qsepvec = arr_alloc(BATCH_SIZE * (p->num->mast_cols+1), double)))
		err_msg("Allocation", "construct_QP", "qsepvec");

	/* added by Yifan to test the impact of different weights of quadratic multiplier */
	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
		avg_sigma += sd_global->quad_v[idx];
	}
	avg_sigma /= BATCH_SIZE;

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (idx = 0; idx < BATCH_SIZE * p->num->mast_cols; idx++)
	{
//		/*qsepvec[idx] = sigma;*/
//		j = idx / BATCH_SIZE;
		/* added by Yifan to test the impact of different weights of qudratic multiplier  */
		qsepvec[idx] = 0.5 * avg_sigma / BATCH_SIZE;
		/* qsepvec[idx] = 0.5 * sd_global->quad_v[j]/ BATCH_SIZE; */
	}

	/* This is for eta column */
	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
		qsepvec[BATCH_SIZE * p->num->mast_cols + idx] = 0.0;
	}

	/* Now copy the Q matrix for QP problem. */
	/* 2011.10.30 - call copy_qp_separable() */
	status = copy_qp_separable(sd_global->batch_problem, qsepvec); //modified by Yifan 09/14/2011
	mem_free(qsepvec);
	//added by Yifan to free memory 09/02/11

	if (status)
	{
		fprintf(stderr, "Failed to copy Q matrix.\n");
		exit(1);
	}
	/* write_prob(sd_global->batch_problem, "quad.lp");*/
}

/****************************************************************************\
  In the regularized QP method, we need to change the rhs of x to d. The 
 A * x = b
 eta + beta * x >= alpha
 Since x = xbar + d, the corresponding changes will therefore be:
 A * d = b - A * xbar
 eta + beta * d >= alpha - beta * xbar
 But as long as the incumbent sulotion does not change, b - A * xbar and
 alpha - beta * xbar (for the existing cuts) won't change. So we only need 
 to change it when the incumbent changes. 

 On the other hand, in each iteration, a new cut will be added (and/or
 some cuts may be dropped) and therefore we need to shift the rhs of the 
 added cut from _alpha_ to _alpha - beta * xbar_, which has taken care of 
 in the routine add_cut() in cutsc. We donot need to worry about the shift
 of rhs for the dropped cuts. 

 This function performs the change of rhs when the incumbent changes, as 
 described above. 
 \****************************************************************************/

void change_rhs(prob_type *p, cell_type *c, soln_type *s)
{
	int status = 0;
	int cnt;
	int idx;
	double *rhs1;
	double *rhs;
	int *indices;

#ifdef TRACE
	printf("Inside change_rhs.\n");
#endif

	if (!(rhs1 = arr_alloc(p->num->mast_rows+1, double)))
		err_msg("Allocation", "change_rhs", "rhs1");
	if (!(rhs =
			arr_alloc(p->num->mast_rows+c->cuts->cnt+c->feasible_cuts_added->cnt, double)))
		err_msg("Allocation", "change_rhs", "rhs");
	if (!(indices =
			arr_alloc(p->num->mast_rows+c->cuts->cnt+c->feasible_cuts_added->cnt, int)))
		err_msg("Allocation", "change_rhs", "indices");

	/*** new rhs = b - A * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** new rhs = b - A * xbar *****\n");
#endif
	/* Be careful with the one_norm!! In the TxX() routine, it assumes the 0th
	 element is reserved for the 1_norm, in the returned vector, the T sparse 
	 vector, and the x vector. */

	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		rhs1[cnt + 1] = p->master->rhsx[cnt];

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "rhs[%d] = %f, ", cnt+1, rhs1[cnt+1]);
#endif
	}

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");

	/* Print out the A matrix. */

	fprintf(g_FilePointer, "++++++ A matrix ++++++");
	for (cnt = 1; cnt<=p->A->cnt; cnt++)
	{
		fprintf(g_FilePointer, "cnt = %d, val = %f, row = %d, col = %d\n",
				cnt, p->A->val[cnt], p->A->row[cnt], p->A->col[cnt]);
	}
#endif

	TxX(p->A, s->incumb_x, rhs1);

#ifdef CAL_CHECK
	/* Print out incumbent X. */
	for (cnt=0; cnt<p->num->mast_cols; cnt++)
	printf ("X[%d] = %f, ", cnt+1, s->incumb_x[cnt+1]);
	printf ("\n");

	/* Print out rhs - AxX. */
	for (cnt=1; cnt<=p->num->mast_rows; cnt++)
	printf ("rhs[%d] = %f, ", cnt, rhs1[cnt]);
	printf("\n");

	for (cnt=0; cnt<p->num->mast_cols; cnt++)
	fprintf (g_FilePointer, "X[%d] = %f, ", cnt+1, s->incumb_x[cnt+1]);
	fprintf (g_FilePointer, "\n");

	/* Print out rhs - AxX. */
	for (cnt=1; cnt<=p->num->mast_rows; cnt++)
	fprintf (g_FilePointer, "rhs[%d] = %f, ", cnt, rhs1[cnt]);
	fprintf(g_FilePointer, "\n");
#endif

	/* Reassign values of rhs1 and print out the new rhs. */
	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		indices[cnt] = cnt;
		rhs[cnt] = rhs1[cnt + 1];

#ifdef CAL_CHECK
		fprintf (g_FilePointer, "rhs[%d] = %f, ", cnt+1, rhs[cnt]);
#endif
	}

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/*** new rhs = alpha - beta * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** rhs = alpha - beta * xbar *****");
#endif
	for (cnt = 0; cnt < c->cuts->cnt; cnt++)
	{
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\nCut #%d:\n", cnt);
#endif
		rhs[p->num->mast_rows + cnt] = c->cuts->val[cnt]->alpha;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "row_num = %d, alpha[%d] = %f\n",
				c->cuts->val[cnt]->row_num, cnt, c->cuts->val[cnt]->alpha);
#endif
		for (idx = 1; idx <= p->num->mast_cols; idx++)
		{
#ifdef CAL_CHECK
			fprintf (g_FilePointer, "beta[%d] = %f, ", idx,
					c->cuts->val[cnt]->beta[idx]);
#endif
			rhs[p->num->mast_rows + cnt] -= c->cuts->val[cnt]->beta[idx]
					* s->incumb_x[idx];
		}
		indices[p->num->mast_rows + cnt] = c->cuts->val[cnt]->row_num;
		/* Yifan 04/04/2012 record "alpha - beta x incumb_x" for later use */
		/* in updating "alpha + (1 - t/k) x eta0" */
		c->cuts->val[cnt]->alpha_incumb = rhs[p->num->mast_rows + cnt];
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "\nindices[%d] = %d, rhs[%d] = %f.",
				p->num->mast_rows+cnt, indices[p->num->mast_rows+cnt],
				p->num->mast_rows+cnt, rhs[p->num->mast_rows+cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif
	/**Do the same thing for feasibility cut**/
	/*** new rhs = alpha - beta * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** rhs = alpha - beta * xbar for feasibility cut*****");
#endif
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\nCut #%d:\n", cnt);
#endif
		rhs[p->num->mast_rows + c->cuts->cnt + cnt] =
				c->feasible_cuts_added->val[cnt]->alpha;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "row_num = %d, alpha[%d] = %f\n",
				c->feasible_cuts_added->val[cnt]->row_num, cnt, c->feasible_cuts_added->val[cnt]->alpha);
#endif
		for (idx = 1; idx <= p->num->mast_cols; idx++)
		{
#ifdef CAL_CHECK
			fprintf (g_FilePointer, "beta[%d] = %f, ", idx,
					c->feasible_cuts_added->val[cnt]->beta[idx]);
#endif
			rhs[p->num->mast_rows + c->cuts->cnt + cnt] -=
					c->feasible_cuts_added->val[cnt]->beta[idx]
							* s->incumb_x[idx];
		}
		indices[p->num->mast_rows + c->cuts->cnt + cnt] =
				c->feasible_cuts_added->val[cnt]->row_num;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "\nindices[%d] = %d, rhs[%d] = %f.",
				p->num->mast_rows+c->cuts->cnt+cnt, indices[p->num->mast_rows+c->cuts->cnt+cnt],
				p->num->mast_rows+c->cuts->cnt+cnt, rhs[p->num->mast_rows+c->cuts->cnt+cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/* Now we change the rhs of the master problem. */
	status = change_rhside(c->master,
			p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt,
			indices, rhs); /* 2011.10.30 */

	if (status)
	{
		fprintf(stderr, "Failed to change the rhs in CPLEX.\n");
		exit(1);
	}

	mem_free(rhs1);
	mem_free(rhs);
	mem_free(indices);
}

void change_rhs_b(prob_type *p, cell_type *c, soln_type *s)
{
	int status = 0;
	int cnt;
	double *rhs1;
	double *rhs;
	int *indices;

#ifdef TRACE
	printf("Inside change_rhs_b.\n");
#endif

	if (!(rhs1 = arr_alloc(p->num->mast_rows+1, double)))
		err_msg("Allocation", "change_rhs_b", "rhs1");
	if (!(rhs = arr_alloc(p->num->mast_rows, double)))
		err_msg("Allocation", "change_rhs_b", "rhs");
	if (!(indices = arr_alloc(p->num->mast_rows, int)))
		err_msg("Allocation", "change_rhs_b", "indices");

	/*** new rhs = b - A * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** new rhs = b - A * xbar *****\n");
#endif
	/* Be careful with the one_norm!! In the TxX() routine, it assumes the 0th
	 element is reserved for the 1_norm, in the returned vector, the T sparse 
	 vector, and the x vector. */

	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		rhs1[cnt + 1] = p->master->rhsx[cnt];
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "rhs[%d] = %f, ", cnt+1, rhs1[cnt+1]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	TxX(p->A, s->incumb_x, rhs1);

#ifdef CAL_CHECK
	/* Print out incumbent X. */

	for (cnt=1; cnt<=p->num->mast_cols; cnt++)
	fprintf (g_FilePointer, "X[%d] = %f, ", cnt, s->incumb_x[cnt]);
	fprintf (g_FilePointer, "\n");

	/* Print out rhs - AxX. */

	for (cnt=1; cnt<=p->num->mast_rows; cnt++)
	fprintf (g_FilePointer, "rhs[%d] = %f, ", cnt, rhs[cnt]);
	fprintf(g_FilePointer, "\n");
#endif

	/* Reassign values of rhs1 and print out the new rhs. */
	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		indices[cnt] = cnt;
		rhs[cnt] = rhs1[cnt + 1];
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "rhs[%d] = %f, ", cnt, rhs[cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif
	/* Now we change the rhs of the master problem. */
	status = change_rhside(c->master, p->num->mast_rows, indices, rhs); /* 2011.10.30 */
	if (status)
	{
		fprintf(stderr, "Failed to change the rhs in CPLEX.\n");
		exit(1);
	}

	mem_free(rhs1);
	mem_free(rhs);
	mem_free(indices);
}

/****************************************************************************\
  In the regularized QP method, we need to change the rhs of x to d. The 
 A * x = b
 eta + beta * x >= alpha
 Since x = xbar + d, the corresponding changes will therefore be:
 A * d = b - A * xbar
 beta * d >= alpha - beta * xbar
 But as long as the incumbent sulotion does not change, b - A * xbar won't
 change. So we only need to change it when the incumbent changes. On the
 other hand, we need to change alpha - beta * xbar in every iteration, since
 the cuts are updated (add/drop) in each iteration. 
 
 This function performs the second changes of rhs, alpha - beta * xbar, as 
 described above. 
 \****************************************************************************/
void change_rhs_cuts(prob_type *p, cell_type *c, soln_type *s, cut_type *cuts)
{
	int status = 0;
	int cnt;
	int idx;
	double *rhs;
	int *indices;

#ifdef TRACE
	printf("Inside change_rhs.\n");
#endif

	if (!(rhs = arr_alloc(cuts->cnt, double)))
		err_msg("Allocation", "change_rhs_cuts", "rhs");
	if (!(indices = arr_alloc(cuts->cnt, int)))
		err_msg("Allocation", "change_rhs_cuts", "indices");

#ifdef CAL_CHECK
	/* Print out incumbent X. */

	fprintf(g_FilePointer, "\n***** incumb_x[] in change_rhs_cuts *****\n");
	for (cnt=1; cnt<=p->num->mast_cols; cnt++)
	fprintf (g_FilePointer, "X[%d] = %f, ", cnt, s->incumb_x[cnt]);
	fprintf (g_FilePointer, "\n");
#endif 

	/*** new rhs = alpha - beta * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** rhs = alpha - beta * xbar *****");
#endif
	for (cnt = 0; cnt < cuts->cnt; cnt++)
	{
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\nCut #%d:\n", cnt);
#endif
		rhs[cnt] = cuts->val[cnt]->alpha;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "row_num = %d, alpha[%d] = %f\n",
				cuts->val[cnt]->row_num, cnt, cuts->val[cnt]->alpha);
#endif
		for (idx = 1; idx <= p->num->mast_cols; idx++)
		{
#ifdef CAL_CHECK
			fprintf (g_FilePointer, "beta[%d] = %f, ", idx,
					cuts->val[cnt]->beta[idx]);
#endif
			rhs[cnt] -= cuts->val[cnt]->beta[idx] * s->incumb_x[idx];
		}
		indices[cnt] = cuts->val[cnt]->row_num;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "\nindices[%d] = %d, rhs[%d] = %f.",
				p->num->mast_rows+cnt, indices[cnt],
				p->num->mast_rows+cnt, rhs[cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/* Now we change the rhs of the master problem. */
	status = change_rhside(c->master, c->cuts->cnt, indices, rhs); /* 2011.10.30 */

	if (status)
	{
		fprintf(stderr, "Failed to change the rhs in CPLEX.\n");
		exit(1);
	}

	mem_free(rhs);
	mem_free(indices);
}

/****************************************************************************\
  This function changes the (lower) bounds of the variables, while changing
 from x to d. The lower bounds of d varibles are -xbar (incumbent solution). 
 \****************************************************************************/
void change_bounds(prob_type *p, cell_type *c, soln_type *s)
{
	int status = 0;
	int cnt;
	double *lbounds;
	double *ubounds;
	int *lindices;
	int *uindices;
	char *llu;
	char *ulu;

#ifdef TRACE
	printf ("Inside change bounds.\n");
#endif

	if (!(lbounds = arr_alloc(p->num->mast_cols, double)))
		err_msg("Allocation", "change_bounds", "lbounds");
	if (!(lindices = arr_alloc(p->num->mast_cols, int)))
		err_msg("Allocation", "change_bounds", "lindices");
	if (!(llu = arr_alloc(p->num->mast_cols, char)))
		err_msg("Allocation", "change_bounds", "llu");

	if (!(ubounds = arr_alloc(p->num->mast_cols, double)))
		err_msg("Allocation", "change_bounds", "ubounds");
	if (!(uindices = arr_alloc(p->num->mast_cols, int)))
		err_msg("Allocation", "change_bounds", "uindices");
	if (!(ulu = arr_alloc(p->num->mast_cols, char)))
		err_msg("Allocation", "change_bounds", "ulu");

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/* Change the Upper Bound GJH 06/24/11 */
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
	{
		ubounds[cnt] = c->master->bdu[cnt] - s->incumb_x[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "indices[%d] = %d, lu[%d] = %c, bounds[%d] = %f,\n",cnt, uindices[cnt], cnt, ulu[cnt], cnt, ubounds[cnt]);
#endif
	}

	status = change_bound(c->master, p->num->mast_cols, uindices, ulu, ubounds); /* 2011.10.30 */
	if (status)
	{
		fprintf(stderr, "Failed to change bounds in CPLEX.\n");
		exit(1);
	}

	/* Change the Lower Bound GJH 06/24/11 */
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
	{
		lbounds[cnt] = -s->incumb_x[cnt + 1]; //Yifan, notice this!!!
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "indices[%d] = %d, lu[%d] = %c, bounds[%d] = %f,\n",cnt, lindices[cnt], cnt, llu[cnt], cnt, lbounds[cnt]);
#endif
	}

	status = change_bound(c->master, p->num->mast_cols, lindices, llu, lbounds); /* 2011.10.30 */
	if (status)
	{
		fprintf(stderr, "Failed to change bounds in CPLEX.\n");
		exit(1);
	}

	mem_free(lbounds);
	mem_free(lindices);
	mem_free(llu);
	mem_free(ubounds);
	mem_free(uindices);
	mem_free(ulu);

#ifdef TRACE
	printf ("Exiting change bounds.\n");
#endif

}

/****************************************************************************\
 When we print the final problem, we need to change the rhs from d back to x. The
 Since x = xbar + d, and
 A * d = b - A * xbar
 eta + beta * d >= alpha - beta * xbar
 
 So the corresponding change sholud be the followings:
 A * x = (b - A * xbar) + A * xbar
 eta + beta * x >= (alpha - beta * xbar) + beta * xbar
 
 This function performs the change of rhs when the incumbent changes, as
 described above.
 \****************************************************************************/

void change_rhs_back(prob_type *p, cell_type *c, soln_type *s)
{
	int status = 0;
	int cnt;
	int idx;
	double *rhs1;
	double *rhs;
	int *indices;

#ifdef TRACE
	printf("Inside change_rhs.\n");
#endif

	if (!(rhs1 = arr_alloc(p->num->mast_rows+1, double)))
		err_msg("Allocation", "change_rhs", "rhs1");
	if (!(rhs =
			arr_alloc(p->num->mast_rows+c->cuts->cnt+c->feasible_cuts_added->cnt, double)))
		err_msg("Allocation", "change_rhs", "rhs");
	if (!(indices =
			arr_alloc(p->num->mast_rows+c->cuts->cnt+c->feasible_cuts_added->cnt, int)))
		err_msg("Allocation", "change_rhs", "indices");

	/*** new rhs = b - A * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** new rhs = b - A * xbar *****\n");
#endif
	/* Be careful with the one_norm!! In the TxX() routine, it assumes the 0th
	 element is reserved for the 1_norm, in the returned vector, the T sparse
	 vector, and the x vector. */

	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		rhs1[cnt + 1] = p->master->rhsx[cnt];

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "rhs[%d] = %f, ", cnt+1, rhs1[cnt+1]);
#endif
	}

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");

	/* Print out the A matrix. */

	fprintf(g_FilePointer, "++++++ A matrix ++++++");
	for (cnt = 1; cnt<=p->A->cnt; cnt++)
	{
		fprintf(g_FilePointer, "cnt = %d, val = %f, row = %d, col = %d\n",
				cnt, p->A->val[cnt], p->A->row[cnt], p->A->col[cnt]);
	}
#endif

	//TxX_plus(p->A, s->incumb_x, rhs1);

#ifdef CAL_CHECK
	/* Print out incumbent X. */
	for (cnt=0; cnt<p->num->mast_cols; cnt++)
	printf ("X[%d] = %f, ", cnt+1, s->incumb_x[cnt+1]);
	printf ("\n");

	/* Print out rhs - AxX. */
	for (cnt=1; cnt<=p->num->mast_rows; cnt++)
	printf ("rhs[%d] = %f, ", cnt, rhs1[cnt]);
	printf("\n");

	for (cnt=0; cnt<p->num->mast_cols; cnt++)
	fprintf (g_FilePointer, "X[%d] = %f, ", cnt+1, s->incumb_x[cnt+1]);
	fprintf (g_FilePointer, "\n");

	/* Print out rhs + AxX. */
	for (cnt=1; cnt<=p->num->mast_rows; cnt++)
	fprintf (g_FilePointer, "rhs[%d] = %f, ", cnt, rhs1[cnt]);
	fprintf(g_FilePointer, "\n");
#endif

	/* Reassign values of rhs1 and print out the new rhs. */
	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		indices[cnt] = cnt;
		rhs[cnt] = rhs1[cnt + 1];

#ifdef CAL_CHECK
		fprintf (g_FilePointer, "rhs[%d] = %f, ", cnt+1, rhs[cnt]);
#endif
	}

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/*** new rhs = alpha + beta * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** rhs = alpha - beta * xbar *****");
#endif
	for (cnt = 0; cnt < c->cuts->cnt; cnt++)
	{
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\nCut #%d:\n", cnt);
#endif
		rhs[p->num->mast_rows + cnt] = c->cuts->val[cnt]->alpha;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "row_num = %d, alpha[%d] = %f\n",
				c->cuts->val[cnt]->row_num, cnt, c->cuts->val[cnt]->alpha);
#endif
		for (idx = 1; idx <= p->num->mast_cols; idx++)
		{
#ifdef CAL_CHECK
			fprintf (g_FilePointer, "beta[%d] = %f, ", idx,
					c->cuts->val[cnt]->beta[idx]);
#endif
			//rhs[p->num->mast_rows+cnt] += c->cuts->val[cnt]->beta[idx] * s->incumb_x[idx];
		}
		indices[p->num->mast_rows + cnt] = c->cuts->val[cnt]->row_num;
		/* Yifan 04/04/2012 record "alpha - beta x incumb_x" for later use */
		/* in updating "alpha + (1 - t/k) x eta0" */
		c->cuts->val[cnt]->alpha_incumb = rhs[p->num->mast_rows + cnt];
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "\nindices[%d] = %d, rhs[%d] = %f.",
				p->num->mast_rows+cnt, indices[p->num->mast_rows+cnt],
				p->num->mast_rows+cnt, rhs[p->num->mast_rows+cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif
	/**Do the same thing for feasibility cut**/
	/*** new rhs = alpha - beta * xbar ***/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***** rhs = alpha - beta * xbar for feasibility cut*****");
#endif
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\nCut #%d:\n", cnt);
#endif
		rhs[p->num->mast_rows + c->cuts->cnt + cnt] =
				c->feasible_cuts_added->val[cnt]->alpha;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "row_num = %d, alpha[%d] = %f\n",
				c->feasible_cuts_added->val[cnt]->row_num, cnt, c->feasible_cuts_added->val[cnt]->alpha);
#endif
		for (idx = 1; idx <= p->num->mast_cols; idx++)
		{
#ifdef CAL_CHECK
			fprintf (g_FilePointer, "beta[%d] = %f, ", idx,
					c->feasible_cuts_added->val[cnt]->beta[idx]);
#endif
			//rhs[p->num->mast_rows+c->cuts->cnt+cnt] += c->feasible_cuts_added->val[cnt]->beta[idx] * s->incumb_x[idx];
		}
		indices[p->num->mast_rows + c->cuts->cnt + cnt] =
				c->feasible_cuts_added->val[cnt]->row_num;
#ifdef CAL_CHECK
		fprintf (g_FilePointer, "\nindices[%d] = %d, rhs[%d] = %f.",
				p->num->mast_rows+c->cuts->cnt+cnt, indices[p->num->mast_rows+c->cuts->cnt+cnt],
				p->num->mast_rows+c->cuts->cnt+cnt, rhs[p->num->mast_rows+c->cuts->cnt+cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/* Now we change the rhs of the master problem. */
	status = change_rhside(c->master,
			p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt,
			indices, rhs); /* 2011.10.30 */

	if (status)
	{
		fprintf(stderr, "Failed to change the rhs in CPLEX.\n");
		exit(1);
	}

	mem_free(rhs1);
	mem_free(rhs);
	mem_free(indices);
}

/****************************************************************************\
 This function changes the (lower and upper) bounds of the variables back
 to its origianl values.
 \****************************************************************************/
void change_bounds_back(prob_type *p, cell_type *c, soln_type *s)
{
	int status = 0;
	int cnt;
	double *lbounds;
	double *ubounds;
	int *lindices;
	int *uindices;
	char *llu;
	char *ulu;
    
#ifdef TRACE
	printf ("Inside change bounds back.\n");
#endif
    
	if (!(lbounds = arr_alloc(p->num->mast_cols, double)))
		err_msg("Allocation", "change_bounds", "lbounds");
    if (!(lindices = arr_alloc(p->num->mast_cols, int)))
        err_msg("Allocation", "change_bounds", "lindices");
    if (!(llu = arr_alloc(p->num->mast_cols, char)))
        err_msg("Allocation", "change_bounds", "llu");
    
    if (!(ubounds = arr_alloc(p->num->mast_cols, double)))
        err_msg("Allocation", "change_bounds", "ubounds");
    if (!(uindices = arr_alloc(p->num->mast_cols, int)))
        err_msg("Allocation", "change_bounds", "uindices");
    if (!(ulu = arr_alloc(p->num->mast_cols, char)))
        err_msg("Allocation", "change_bounds", "ulu");
    
#ifdef CAL_CHECK
    fprintf(g_FilePointer, "\n");
#endif
    
    for (cnt = 0; cnt < p->num->mast_cols; cnt++)
    {
        ubounds[cnt] = c->master->bdu[cnt];
        uindices[cnt] = cnt;
        ulu[cnt] = 'U';
#ifdef CAL_CHECK
        fprintf(g_FilePointer, "indices[%d] = %d, lu[%d] = %c, bounds[%d] = %f,\n",cnt, uindices[cnt], cnt, ulu[cnt], cnt, ubounds[cnt]);
#endif
    }
    
	status = change_bound(c->master, p->num->mast_cols, uindices, ulu, ubounds);
	if (status)
	{
		fprintf(stderr, "Failed to change bounds in CPLEX.\n");
		exit(1);
	}
    
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
	{
		lbounds[cnt] = c->master->bdl[cnt]; //Yifan, notice this!!!
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "indices[%d] = %d, lu[%d] = %c, bounds[%d] = %f,\n",cnt, lindices[cnt], cnt, llu[cnt], cnt, lbounds[cnt]);
#endif
	}
    
	status = change_bound(c->master, p->num->mast_cols, lindices, llu, lbounds);
	if (status)
	{
		fprintf(stderr, "Failed to change bounds in CPLEX.\n");
		exit(1);
	}
    
	mem_free(lbounds);
	mem_free(lindices);
	mem_free(llu);
	mem_free(ubounds);
	mem_free(uindices);
	mem_free(ulu);
    
#ifdef TRACE
	printf ("Exiting change bounds back.\n");
#endif
    
}
void add_flip_bounds(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s, double * incumbent_x_k, int *num_flip)
{
    int cnt,rhs_ind[1];
    double rhs_val[1];
    rhs_val[0] = 0;
    rhs_ind[0] = 2;
    for (cnt = 0; cnt < p->num->mast_cols; cnt++)
    {
        if (p->master->typex[cnt]=='B' && incumbent_x_k[cnt+1]==1) {
            CPXchgcoef(env, c->master->lp, 2, cnt, 1.0);
            rhs_val[0]++;
        }
        else{
            CPXchgcoef(env, c->master->lp, 2, cnt, 0.0);
        }
    }
    
    rhs_val[0] = rhs_val[0] - *num_flip;
    //rhs_val[0] = rhs_val[0] - 1;
    CPXchgrhs(env, c->master->lp, 1, rhs_ind, rhs_val);
    //print_problem(c->master, "fcssn.lp");
}

void add_box_bounds(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s)
{
	int status = 0;
	int cnt;
	double *lbounds;
	double *ubounds;
	int *lindices;
	int *uindices;
	char *llu;
	char *ulu;
    
#ifdef TRACE
	printf ("Inside add_box_bounds.\n");
#endif
    
	if (!(lbounds = arr_alloc(p->num->mast_cols, double)))
		err_msg("Allocation", "change_bounds", "lbounds");
    if (!(lindices = arr_alloc(p->num->mast_cols, int)))
        err_msg("Allocation", "change_bounds", "lindices");
    if (!(llu = arr_alloc(p->num->mast_cols, char)))
        err_msg("Allocation", "change_bounds", "llu");
    
    if (!(ubounds = arr_alloc(p->num->mast_cols, double)))
        err_msg("Allocation", "change_bounds", "ubounds");
    if (!(uindices = arr_alloc(p->num->mast_cols, int)))
        err_msg("Allocation", "change_bounds", "uindices");
    if (!(ulu = arr_alloc(p->num->mast_cols, char)))
        err_msg("Allocation", "change_bounds", "ulu");
    
#ifdef CAL_CHECK
    fprintf(g_FilePointer, "\n");
#endif
    
    for (cnt = 0; cnt < p->num->mast_cols; cnt++)
    {
        if (p->master->typex[cnt]=='I') {
            ubounds[cnt] = min(s->mip_x[cnt+1] + sd_global->config.BOX_SIZE, c->master->bdu[cnt]);
            uindices[cnt] = cnt;
            ulu[cnt] = 'U';
        }
        
#ifdef CAL_CHECK
        fprintf(g_FilePointer, "indices[%d] = %d, lu[%d] = %c, bounds[%d] = %f,\n",cnt, uindices[cnt], cnt, ulu[cnt], cnt, ubounds[cnt]);
#endif
    }
    
	status = change_bound(c->master, p->num->mast_cols, uindices, ulu, ubounds);
	if (status)
	{
		fprintf(stderr, "Failed to change bounds in CPLEX.\n");
		exit(1);
	}
    
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
	{
        if (p->master->typex[cnt]=='I') {
            lbounds[cnt] = max(s->mip_x[cnt+1] - sd_global->config.BOX_SIZE,c->master->bdl[cnt]);
            lindices[cnt] = cnt;
            llu[cnt] = 'L';
        }
        
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "indices[%d] = %d, lu[%d] = %c, bounds[%d] = %f,\n",cnt, lindices[cnt], cnt, llu[cnt], cnt, lbounds[cnt]);
#endif
	}
    
	status = change_bound(c->master, p->num->mast_cols, lindices, llu, lbounds);
	if (status)
	{
		fprintf(stderr, "Failed to change bounds in CPLEX.\n");
		exit(1);
	}
    
	mem_free(lbounds);
	mem_free(lindices);
	mem_free(llu);
	mem_free(ubounds);
	mem_free(uindices);
	mem_free(ulu);
    
#ifdef TRACE
	printf ("Exiting add_box_bounds.\n");
#endif
    
}

void update_cont_bounds(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s, double *incumbent_x_k,double alpha, int mip_idx)
{
    int cnt, status;
    char **r_names;
	double xrhs[1], coef[2];
	char xsense[1] = { 'G' };
	int rmatbeg[1] = { 0 };
    int coef_col[2];
    /* Give the lb a name */
	r_names = arr_alloc(2, string);
	r_names[0] = arr_alloc(NAME_SIZE, char);


    
    /* update existing capacity upper bound */
    for (cnt = 0; cnt < p->num->mast_cols; cnt++)
    {
        if (p->master->typex[cnt]=='B') {
            if (incumbent_x_k[cnt+1]==1) {
                CPXchgcoef(env, c->master->lp, 2 + cnt + 1 - 89, cnt, -(1+alpha) * incumbent_x_k[cnt+1-89]);
            }
            else{
                CPXchgcoef(env, c->master->lp, 2 + cnt + 1 - 89, cnt, -1008);
            }
        }
    }
    //print_problem(c->master, "after_update_cont_ubds.lp");
    
    /* add(update) capacity lower bound */
    if (mip_idx > 1) {
        for (cnt = 0; cnt < p->num->mast_cols; cnt++)
        {
            if (p->master->typex[cnt]=='B') {
                if (incumbent_x_k[cnt+1]==1) {
                    status = CPXchgcoef(env, c->master->lp, 2 + cnt + 1, cnt, -(1-alpha) * incumbent_x_k[cnt+1-89]);
                }
                else{
                    status = CPXchgcoef(env, c->master->lp, 2 + cnt + 1, cnt, 0);
                }
            }
        }
    }
    else{
        for (cnt = 0; cnt < p->num->mast_cols; cnt++)
        {
            if (p->master->typex[cnt]=='B') {
                if (incumbent_x_k[cnt+1]==1) {
                    strcpy(r_names[0], "ULB  ");
                    r_names[0][3] = '0' + (cnt - 89 + 1) / 10 % 10;
                    r_names[0][4] = '0' + (cnt - 89 + 1) / 1 % 10;
                    coef_col[0] = cnt - 89;
                    coef_col[1] = cnt;
                    coef[0] = 1;
                    coef[1] = -(1-alpha) * incumbent_x_k[cnt-89+1];
                    xrhs[0] = 0.0;
                    status = CPXaddrows(env, c->master->lp, 0, 1, 2, xrhs, xsense, rmatbeg,
                               coef_col, coef, NULL, r_names);
                    print_problem(c->master, "after_update_cont_bds.lp");
                }
                else{
                    strcpy(r_names[0], "ULB  ");
                    r_names[0][3] = '0' + (cnt - 89 + 1) / 10 % 10;
                    r_names[0][4] = '0' + (cnt - 89 + 1) / 1 % 10;
                    coef_col[0] = cnt - 89;
                    coef_col[1] = cnt ;
                    coef[0] = 1;
                    coef[1] = 0.0;
                    xrhs[0] = 0.0;
                    
                    status = CPXaddrows(env, c->master->lp, 0, 1, 2, xrhs, xsense, rmatbeg,
                               coef_col, coef, NULL, r_names);
                    print_problem(c->master, "after_update_cont_bds.lp");
                }
            }
        }
    }

    //print_problem(c->master, "after_update_cont_bds.lp");
    
}
