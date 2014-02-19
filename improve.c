/***********************************************************************\
**
 ** improve.c
 **
 ** This file contains the functions used in assessing any improvement
 ** in the solution.  In step 4 of the SD algorithm, the incumbent x
 ** may be updated if an improvement is shown.  Otherwise, the incumbent
 ** x remains the same, as does the incumbent cut.
 **
 **
 ** History:
 **   10 Mar 1992 - <Jason Mai> - created.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "utility.h"
#include "theta.h"
#include "improve.h"
#include "omega.h"
#include "sdglobal.h"

/***********************************************************************\ 
 ** This function determines whether the "stagewise descent property" is 
 ** satisified.  If the current approximation of f_k gives a lower difference
 ** between the candidate and incumbent x than the previous approximation
 ** gave, then the incumbent x is updated to the candidate x, and the
 ** reference to the incumbent cut is updated as well.  The function returns
 ** TRUE if the incumbent was updated; FALSE otherwise.
 \***********************************************************************/
BOOL check_improvement(sdglobal_type* sd_global, prob_type *prob,
		cell_type *cell, soln_type *soln, double *phi)
{
	double candid_est;
#ifdef CAL_CHECK
	int cnt;
#endif
	/*int		best;*/
#ifdef DEBUG
	FILE *gamma;
#endif
#ifdef TRACE
	printf("Inside check_improvement\n");
#endif

	/* Calculate height at new candidate x with newest cut included */
	/*
	 candid_est = f_k(soln->candid_x, prob->c, cell->cuts,
	 cell->theta, prob->num, cell, &best);*/

	/* Yifan 03/20/2012 Test for gamma issues*/

	candid_est = max_cut_height(sd_global, cell->cuts, soln->candid_x, cell,
			prob->num);
	candid_est += CxX(prob->c, soln->candid_x, prob->num->mast_cols);
	/*candid_est = soln->opt_value;*/

	/* Calculate height at current incumbent x with newest cut. */
	/*soln->incumb_est = f_k(soln->incumb_x, prob->c, cell->cuts,
	 cell->theta, prob->num, cell, &best);*/

	soln->incumb_est = max_cut_height(sd_global, cell->cuts, soln->incumb_x,
			cell, prob->num);
	soln->incumb_est += CxX(prob->c, soln->incumb_x, prob->num->mast_cols);

#ifdef CAL_CHECK
	printf("soln->incumb_est:%f\n",soln->incumb_est);
	print_vect(soln->incumb_x, prob->num->mast_cols, "incumb_x");
	for (cnt=0; cnt<cell->cuts->cnt; cnt++)
	{
		printf("alpha[%d]:%f\t",cnt,cell->cuts->val[cnt]->alpha);
		printf("cut_height[%d]:%f\t",cnt,cut_height(cell->cuts->val[cnt], soln->incumb_x, cell, prob->num));
		print_vect(cell->cuts->val[cnt]->beta, prob->num->mast_cols, "beta");
	}
	printf("max cut height: %f\n",soln->incumb_est);
#endif

	if (cell->k == 2)
	{ /*  SS messing around  */
		phi[0] = soln->incumb_est; /*  SS messing around  */
		phi[1] = soln->incumb_est; /*  SS messing around  */
	}

#ifdef DEBUG
	gamma = fopen("gamma.out", "a");
	printf("\n");
	print_vect(soln->candid_x, prob->num->mast_cols, "Candidate X");
	print_vect(soln->incumb_x, prob->num->mast_cols, "Incumbent X");
	printf("\nImprov: f_k(candid_x) = %f. f_k(incumb_x) = %f. G = %f\n",
			candid_est, soln->incumb_est, soln->gamma);
	fprintf(gamma, "%f\n", soln->gamma);
	fclose(gamma);
#endif

	/*
	 fprintf(g_FilePointer, "\n***** In improve.c *****\n");
	 fprintf(g_FilePointer, "candid_est = %f, incumb_est = %f, gamma = %f\n",
	 candid_est, soln->incumb_est, soln->gamma); 
	 */
	/* Print out d_norm's for the calculation check when updating quad_scalar. zl
	 fprintf(g_FilePointer, "\n*** in improve.c, check_improvement. ***\n");
	 fprintf(g_FilePointer, "norm_d_k = %f, norm_d_k_1 = %f, quad_scalar = %f\n", 
	 soln->norm_d_k, soln->norm_d_k_1, cell->quad_scalar);
	 */

	/* If we saw considerable improvement, then change the incumbent */
	if ((candid_est - soln->incumb_est) < (sd_global->config.R * soln->gamma))
	{

		new_incumbent(prob, cell, soln, candid_est);

		if (cell->k > 1 && soln->norm_d_k > sd_global->config.TOLERANCE)
			if (soln->norm_d_k >= sd_global->config.R3 * soln->norm_d_k_1)
			{
				/*printf("soln->norm_d_k is %f\n",soln->norm_d_k_1);
				 printf("soln->norm_d_k_1 is %f\n",soln->norm_d_k_1);*/
				cell->quad_scalar *= sd_global->config.R2 * sd_global->config.R3
						* soln->norm_d_k_1 / soln->norm_d_k;

				cell->quad_scalar =
						min(sd_global->config.MAX_QUAD_SCALAR, cell->quad_scalar );
				/* Yifan 03/19/2012 Make sure the quad sclar will never be smaller than 1.0*/
				cell->quad_scalar =
						max(sd_global->config.MIN_QUAD_SCALAR, cell->quad_scalar);

				if (sd_global->config.MASTER_TYPE == SDQP)
					construct_QP(prob, cell, cell->quad_scalar);
			}
		/*
		 fprintf(g_FilePointer, "new incumbent, quad_scalar = %f\n",
		 cell->quad_scalar);
		 */
		/* As mentioned above, in regularized QP method, we need to switch the 
		 master problem from x to d.  The rhs is to be changed as: 
		 new rhs = b - AxX, (for the original constraints) 
		 and 
		 new rhs = alpah - beta * X, (for cuts)
		 and 
		 new lower bound = -X, 
		 They are updated each time we change the incumbent. zl 
		 */
		if (sd_global->config.MASTER_TYPE == SDQP)
		{
			change_rhs(prob, cell, soln);
			change_bounds(prob, cell, soln);
		}

		phi[0] = soln->incumb_est; /*  SS messing around  */
		phi[1] = sd_global->config.SMOOTH_PARM * soln->incumb_est
				+ (1 - sd_global->config.SMOOTH_PARM) * phi[1]; /*  SS messing around  */
		phi[2] = 0.0;
		phi[3] = (double) cell->k;
		/*printf("*****   NEW INCUMBENT   *****\n");*//*comment out by Yifan to see impact ratio*/
		printf("+");
		fflush(stdout);

		/* Added to assure the validity of zero as subproblem LB. 
		 * Lei 09/13/05 */
		if (sd_global->config.SUB_LB_CHECK)
			printf("Subproblem LB = %lf\n", soln->sub_lb_checker);

		soln->norm_d_k_1 = soln->norm_d_k; /* zl. 06/10/02 */
		return TRUE;
	}
	else /* all the cuts are falling away a step */
	{
		/* Update quad_scalar when no incumbent is found. */
		cell->quad_scalar =
				min(sd_global->config.MAX_QUAD_SCALAR, cell->quad_scalar / sd_global->config.R2);

		if (sd_global->config.MASTER_TYPE == SDQP)
			construct_QP(prob, cell, cell->quad_scalar); //modified by Yifan 09/14/2011
		/*
		 fprintf(g_FilePointer, "not new incumbent, quad_scalar = %f\n",
		 cell->quad_scalar);
		 */
		phi[0] = sd_global->config.SMOOTH_PARM * soln->incumb_est
				+ (1 - sd_global->config.SMOOTH_PARM) * phi[0]; /*  SS messing around  */
		phi[1] = sd_global->config.SMOOTH_PARM * soln->incumb_est
				+ (1 - sd_global->config.SMOOTH_PARM) * phi[1]; /*  SS messing around  */
		soln->incumb_stdev *= (cell->k - 1) / (double) (cell->k);
		/*
		 soln->incumb_est *= (cell->k - 1) / (double) (cell->k);
		 */
		soln->norm_d_k_1 = soln->norm_d_k; /* zl. 06/10/02 */
		/*printf("cell->quad_scalar is : %f\n",cell->quad_scalar);*/
		return FALSE;
	}

}

BOOL new_check_improvement(sdglobal_type* sd_global, prob_type *prob,
		cell_type *cell, soln_type *soln, double *phi)
{
	double candid_est;

	/*int		best;*/
#ifdef DEBUG
	FILE *gamma;
#endif
#ifdef TRACE
	printf("Inside new_check_improvement\n");
#endif

	/* Calculate height at new candidate x with newest cut included */
	/*
	 candid_est = f_k(soln->candid_x, prob->c, cell->cuts,
	 cell->theta, prob->num, cell, &best);*/

	/* Yifan 03/20/2012 Test for gamma issues*/

	candid_est = max_cut_height(sd_global, cell->cuts, soln->candid_x, cell,
			prob->num);
	candid_est += CxX(prob->c, soln->candid_x, prob->num->mast_cols);
	/*candid_est = soln->opt_value;*/

	/* Calculate height at current incumbent x with newest cut. */
	/*soln->incumb_est = f_k(soln->incumb_x, prob->c, cell->cuts,
	 cell->theta, prob->num, cell, &best);*/

	soln->incumb_est = max_cut_height(sd_global, cell->cuts, soln->incumb_x,
			cell, prob->num);
	soln->incumb_est += CxX(prob->c, soln->incumb_x, prob->num->mast_cols);

	/* If we saw considerable improvement, then change the incumbent */
	if ((candid_est - soln->incumb_est) < (sd_global->config.R * soln->gamma))
	{
		calc_quad_scalar(sd_global, prob, cell, soln, candid_est);
		new_incumbent(prob, cell, soln, candid_est);

		if (sd_global->config.MASTER_TYPE == SDQP)
		{
			change_rhs(prob, cell, soln);
			change_bounds(prob, cell, soln);
		}

		/* printf("*****   NEW INCUMBENT   *****\n");*/
		printf("+");
		fflush(stdout);

		return TRUE;
	}
	else /* all the cuts are falling away a step */
	{
		/* Update quad_scalar when no incumbent is found. */
		cell->quad_scalar =
				min(sd_global->config.MAX_QUAD_SCALAR, 2 * cell->quad_scalar);

		if (sd_global->config.MASTER_TYPE == SDQP)
			construct_QP(prob, cell, cell->quad_scalar);
		return FALSE;
	}

}
/***********************************************************************\  
 ** If it has been decided that the incumbent solution should change, 
 ** several values must be updated, including the incumbent X, the 
 ** incumbent cut, the estimate at the incumbent and the standard 
 ** deviation of the estimate, as well as the last_update counter.
 ** This function performs these updates.  It assumes that the candidate 
 ** solution becomes the incumbent solution, with a height corresponding 
 ** to the _est_ parameter.
 \***********************************************************************/
void new_incumbent(prob_type *p, cell_type *c, soln_type *s, double est)
{
	vector X;
	double val, stdev;
	int obs, idx;
	int count;
	i_type i;

#ifdef TRACE
	printf("Inside new_incumbent\n");
#endif 

	/* Copy over information about the new incumbent */
	copy_arr(s->incumb_x, s->candid_x, p->num->mast_cols)
	s->incumb_cut = c->cuts->cnt - 1;
	s->incumb_est = est;
	s->last_update = c->k;
	s->incumb_k = c->k;

	/* 
	 ** The cut consists of an average of many pi X omega products.  We
	 ** use the mean and standard deviation of this collection of products
	 ** in the check for optimality, so they are updated here.  We loop 
	 ** through the istar indices of new incumbent cut to find the stdev.
	 */
	count = 0;
	stdev = 0.0;
	X = s->incumb_x;

	if (c->cuts->cnt > 0)
	{
		for (obs = 0; obs < c->cuts->val[s->incumb_cut]->omega_cnt; obs++)
			if (valid_omega_idx(s->omega, obs))
			{
				/* Find the pi from the cut's argmax for this observation */
				i.sigma = c->cuts->val[s->incumb_cut]->istar[obs];
				i.delta = c->sigma->lamb[i.sigma];

				/*
				 ** Couldn't this be done only in the optimality check, so we 
				 ** don't do it every time we change our minds about the incumbent ?
				 ** Which happens more often -- an optimality check or a change of
				 ** the incumbent x?
				 ** What about when the the incumbent cut is updated, but X is not?
				 */

				/* Calculate the height for this observation and dual vector */
				val = c->sigma->val[i.sigma].R + s->delta->val[i.delta][obs].R;
				for (idx = 0; idx < p->num->nz_cols; idx++)
					val -= c->sigma->val[i.sigma].T[idx]
							* X[c->sigma->col[idx]];
				for (idx = 0; idx < p->num->rv_cols; idx++)
					val -= s->delta->val[i.delta][obs].T[idx]
							* X[s->delta->col[idx]];

				stdev += s->omega->weight[obs] * SQR(s->incumb_est - val);
				count += s->omega->weight[obs];
			}
		s->incumb_stdev = sqrt(stdev / (double) count);
	}

	s->incumbent_change = TRUE;
#ifdef DEBUG
	printf("\nMade a new incumbent, since it qualified");
	print_cut(c->cuts, p->num, s->incumb_cut);
#endif

}

/* This function loops through a set of cuts and find the highest cut height
 at the specified position x */
double max_cut_height(sdglobal_type* sd_global, cut_type *cuts, vector X,
		cell_type *c, num_type *num)
{
	double Sm;
	double ht;
	int cnt;

	Sm = cut_height(sd_global, cuts->val[0], X, c, num);
	for (cnt = 1; cnt < cuts->cnt; cnt++)
	{
		ht = cut_height(sd_global, cuts->val[cnt], X, c, num);
		if (Sm < ht)
			Sm = ht;
	}
	return Sm;
}

