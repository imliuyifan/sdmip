/***********************************************************************\
**
 ** cuts.c
 **
 **
 ** A cut is formed according to the following algorithm:
 ** ----------------------------------------------------
 **  for all observations of omega
 **    for all dual vectors
 **      calculate pix(R - TxX)
 **      track the max
 **    record the max as an istar
 **    add the pi[istar]xR to alpha
 **    add the pi[istar]xT to beta
 **  divide alpha & beta by the # of observations made
 ** ----------------------------------------------------
 **
 **
 ** This file contains the following functions for forming a cut:
 ** ------------------------------------------------------------
 ** form_new_cut()
 ** form_incumbent_cut()
 ** stochastic_updates()
 ** SD_cut()
 ** compute_istar()
 ** new_cut()
 ** free_cut()
 ** print_cut()
 ** new_cuts()
 ** free_cuts()
 ** print_cut_info() // printing cuts for cut index checking. zl 
 **
 **
 ** History:
 **   20 Jan 1992 - <Jason Mai> - created.
 **   14 Dec 1992 - <Jason Mai> - fixed compilation errors.
 **   22 Feb 1992 - <Jason Mai> - updated with new structures (esp. num_type).
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "solver.h"
#include "utility.h"
#include "theta.h"
#include "testout.h"
#include "subprob.h"
#include "sigma.h"
#include "delta.h"
#include "omega.h"
#include "lambda.h"
#include "log.h"
#include "cuts.h"
#include "sdglobal.h"

#include <float.h>
#include <stdlib.h>
#include <time.h>

/*
 **  WAIT!  cuts->cnt is NOT the same as num_samples, or cell->k.
 ** you must use cuts->cnt when accessing the cuts array, since it
 ** tells you the place of the next cut, somewhere after all the
 ** ancestors' cuts.  But, when weighting cuts, need num_samples or
 ** cell->k, to get the number of samples during this cell.
 **
 ** Did you remember this?
 ** I think so.
 **
 ** Don't forget to multiply sigma by omega->weight also!
 */

/***********************************************************************\
** Once an observation of omega has been generated and the
 ** corresponding subproblem has been solved, a new cut may
 ** be formed for the master problem.  This function will
 ** update all the data (concerning omega and the dual vector)
 ** and form a new cut based upon it.  The new cut is added
 ** to the array of cuts, and the number of cuts is incremented.
 \***********************************************************************/
void form_new_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int omeg_idx, BOOL new_omega)
{
	one_cut *cut;

#ifdef DEBUG
	int idx;
#endif

#ifdef TRACE
	printf("Inside form_new_cut; %d cuts now\n", c->cuts->cnt);
#endif

	cut = new_cut(p->num->mast_cols, s->omega->most, c->k);

	stochastic_updates(sd_global, c, c->lambda, c->sigma, s->delta, s->omega,
			p->num, p->Rbar, p->Tbar, c->subprob, s->Pi, omeg_idx, new_omega);

	SD_cut(sd_global, c->sigma, s->delta, s->omega, p->num, cut, s->candid_x,
			s->pi_ratio, s->max_ratio, s->min_ratio, c->k,
			s->dual_statble_flag);

	add_cut(sd_global, cut, p, c, s);

	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl 

	 print_cut_info(c, p->num, "After adding a cut.");  
	 */

#ifdef RUN
	/*
	 if (c->k == 10)
	 {
	 print_vect(s->candid_x, p->num->mast_cols, "Candidate X");

	 for(idx = 0; idx < c->lambda->cnt; idx++)
	 for(obs = 0; obs < s->omega->most; obs++)
	 if (valid_omega_idx(s->omega, obs))
	 print_delta(s->delta, p->num, idx, obs);

	 for(idx = 0; idx < c->lambda->cnt; idx++)
	 print_lambda(c->lambda, p->num, idx);

	 for(obs = 0; obs < s->omega->most; obs++)
	 if (valid_omega_idx(omega, obs))
	 print_omega(s->omega, p->num, obs);

	 for(idx = 0; idx < c->sigma->cnt; idx++)
	 print_sigma(c->sigma, p->num, idx);

	 for(idx=0; idx<c->cuts->cnt; idx++)
	 print_cut(c->cuts, p->num, idx);
	 }
	 */
#endif

#ifdef DEBUG
	for(idx=0;idx<c->cuts->cnt;idx++)
	print_cut(c->cuts, p->num, idx);
	printf("\n");
#endif

}

/* modified by Yifan 2013.9.26 */
void remove_all_cuts(sdglobal_type* sd_global, prob_type *p, cell_type *c,
                     soln_type *s)
{
    int cut_idx;
#ifdef TRACE
	printf("Inside remove_all_cuts; %d cuts now\n", c->cuts->cnt);
#endif
    for (cut_idx = c->cuts->cnt - 1; cut_idx >= 0; cut_idx--) {
        if (c->cuts->val[cut_idx]->subfeaflag == TRUE)
        {
            drop_cut(cut_idx, p, c, s);
            decrease_feacut_rownum(c);
        }
    }
#ifdef TRACE
	printf("Exiting remove_all_cuts; %d cuts now\n", c->cuts->cnt);
#endif
}


void form_ip_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
                  soln_type *s, int omeg_idx, BOOL new_omega, int mip_idx)
{
	one_cut *cut;
    
#ifdef DEBUG
	int idx;
#endif
    
#ifdef TRACE
	printf("Inside form_ip_cut; %d cuts now\n", c->cuts->cnt);
#endif
    
    // push Benders' cut down also, but let's push them down gentally
	cut = new_cut(p->num->mast_cols, s->omega->most, c->k);
    
	stochastic_updates(sd_global, c, c->lambda, c->sigma, s->delta, s->omega,
                       p->num, p->Rbar, p->Tbar, c->subprob, s->Pi, omeg_idx, new_omega);
    
	SD_cut(sd_global, c->sigma, s->delta, s->omega, p->num, cut, s->mip_x,
           s->pi_ratio, s->max_ratio, s->min_ratio, c->k,
           s->dual_statble_flag);
    sd_global->config.MASTER_TYPE = SDLP;
    cut->beta[0] = 1.0;
	add_cut(sd_global, cut, p, c, s);
    sd_global->config.MASTER_TYPE = SDQP;
    
	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl
     
	 print_cut_info(c, p->num, "After adding a cut.");
	 */
#ifdef TRACE
	printf("Exiting form_ip_cut; %d cuts now\n", c->cuts->cnt);
#endif
    
}

/***********************************************************************\
 ** Once an observation of omega has been generated and the
 ** corresponding primal subproblem is infeasible. A feasibility cut will
 ** be generated for the master problem.  This function will
 ** update all the data (concerning omega and the extreme directions)
 ** and form new cuts based upon it.  The new cuts are added
 ** to the pool of feasibility cuts, and the number of cuts is incremented.
 \***********************************************************************/
void form_fea_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int omeg_idx, vector x_k, BOOL new_omega)
{
	BOOL new_sigma;
	int cut_cnt_added;
	int start, end;
	int idx;

#ifdef TRACE
	printf("Inside form_fea_cut; %d cuts now\n", c->cuts->cnt);
#endif    

	new_sigma = stochastic_updates(sd_global, c, c->feasible_lambda,
			c->feasible_sigma, s->feasible_delta, s->omega, p->num, p->Rbar,
			p->Tbar, c->subprob, s->Pi, omeg_idx, new_omega);

	FEA_cut(sd_global, c, s, c->feasible_sigma, s->feasible_delta, s->omega,
			p->num, s->omega->most, s->dual_statble_flag, new_omega, new_sigma);

//	cut_cnt_pool = c->feasible_cuts_pool->cnt;

#ifdef DEBUG
	printf("cut_cnt_pool is %d\n",cut_cnt_pool);
	for(idx=0; idx < c->feasible_cuts_pool->cnt; idx++)
	print_cut(c->feasible_cuts_pool, p->num, idx);
	printf("\n");
#endif

	start = c->feasible_cuts_added->cnt;
	cut_cnt_added = FEA_cut_check_add(sd_global, c, p, s, x_k);
	end = cut_cnt_added;

#ifdef DEBUG
	printf("cut_cnt_added is %d\n",cut_cnt_added);
	for(idx=0; idx < c->feasible_cuts_added->cnt; idx++)
	print_cut(c->feasible_cuts_added, p->num, idx);
	printf("\n");
#endif

	if (end > start)
	{
		for (idx = start; idx < end; idx++)
		{
			add_cut_to_master(sd_global, c->feasible_cuts_added->val[idx], p, c,
					s, idx);
		}
		update_dual_size(c, s, p);
	}

#ifdef TRACE
	printf("Exiting form_fea_cut; %d cuts now\n", c->cuts->cnt);
#endif

}

void update_dual_size(cell_type *c, soln_type *s, prob_type *p)
{

	mem_free(s->Master_pi);
	/* Yifan 03/11/2012 the zero-th element is saved for one-norm */
	/* Yifan 03/11/2012 new room in dual value is needed for new feasibility cut*/
	s->Master_pi =
			arr_alloc(p->num->mast_rows+1+p->num->max_cuts+c->feasible_cuts_added->cnt,double);
}

/***********************************************************************\
** As iterations proceed, the cut which was formed based on the
 ** incumbent X vector is periodically re-evaluated.  If tau iterations
 ** have passed since the last update, then the cut is reformed.  Or,
 ** if the value of the incumbent cut is less than the value of f_k
 ** at the incumbent X, the incumbent is reformed.  This function
 ** will solve an additional subproblem (based upon the incumbent X
 ** and the most recent observation of omega), add the dual solution
 ** to the data structures, and re-form the incumbent cut (thus
 ** replacing it in the array of cuts).  Note that the updates do
 ** NOT treat the omega as a new observation!
 \***********************************************************************/
BOOL form_incumb_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int omeg_idx)
{

	one_cut *cut;
	clock_t start, end; /* Recording time spent on argmax
	 procedures. zl, 06/30/04. */
	BOOL incumb_change = FALSE;

#ifdef TRACE
	printf("Inside form_incumb_cut, i_k=%d\n", s->incumb_cut);
#endif

	/*
	 ** Calculate the height of f_k at the incumbent X according to the
	 ** the cuts.  We expect the incumbent cut to be the one which has
	 ** the highest value at this point.
	 **
	 int		best;
	 double	object_est;
	 double	incumb_est;
	 
	 object_est = f_k(s->incumb_x, p->c, c->cuts, c->theta, p->num, c, &best);
	 incumb_est = cut_height(c->cuts->val[s->incumb_cut], s->incumb_x, c, p->num);
	 incumb_est += CxX(p->c, s->incumb_x, p->num->mast_cols);
	 **
	 #ifdef DEBUG
	 printf("object_est = %f.  incumb_est = %f.\n", object_est, incumb_est);
	 #endif
	 **
	 ** If the height of the objective function approximation, f_k, at 
	 ** the incumbent X is greater than the height of the incumbent cut
	 ** at the incumbent X, then the _best_ cut will not be the _incumb_cut_,
	 ** and we must reform the incumbent cut.  Also, we will reform 
	 ** the cut if tau iterations have passed since the last update.
	 **
	 ** What if two cuts tie, and incumb_cut is the 2nd one?  best != incumb
	 ** but we don't need to reform it!!!  Need height of incumbent...
	 **
	 if (object_est - incumb_est > sd_global->config.TOLERANCE ||
	 c->k - s->last_update == p->tau)
	 **
	 */

	if (c->k - s->last_update == p->tau)
	{

		solve_subprob(sd_global, p, c, s, s->incumb_x, omeg_idx);

		/* Record time spent on argmax procedures without counting the time
		 on solving the subproblem LP. zl, 06/30/04. */
		start = clock(); /* zl, 06/30/04. */

		if (c->subprob->feaflag == TRUE)
		{

			/*modified by Yifan to avoid dropping feasibility cuts 02/09/2012*/
			/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
			if (c->cuts->cnt >= p->num->max_cuts)
			{
				drop_cut(s->incumb_cut, p, c, s);
				decrease_feacut_rownum(c);
			}

		}

		/* zl, 06/30/04. */
		end = clock();

		s->run_time->argmax_iter += ((double) (end - start)) / CLOCKS_PER_SEC;

		if (c->subprob->feaflag == FALSE)
		{
			incumb_change = TRUE;
		}
		else
		{
			start = clock(); /* zl, 06/30/04. */
			cut = new_cut(p->num->mast_cols, s->omega->most, c->k);
			cut->is_incumbent = TRUE; /*added by Yifan 02/02/2012 indentify a new incumbent cut*/
			stochastic_updates(sd_global, c, c->lambda, c->sigma, s->delta,
					s->omega, p->num, p->Rbar, p->Tbar, c->subprob, s->Pi,
					omeg_idx, FALSE);
			SD_cut(sd_global, c->sigma, s->delta, s->omega, p->num, cut,
					s->incumb_x, s->pi_ratio, s->max_ratio, s->min_ratio, c->k,
					s->dual_statble_flag);
			s->incumb_cut = add_cut(sd_global, cut, p, c, s);
			s->last_update = c->k;
			/* zl, 06/30/04. */
			end = clock();
			s->run_time->argmax_iter += ((double) (end - start))
					/ CLOCKS_PER_SEC;

			/*added to record subproblem objective estimate by Yifan 02/01/12*/
			//s->subobj_est[c->k] = cut->alpha - CxX(cut->beta, s->incumb_x, p->num->mast_cols);
			s->subobj_est = cut_height(sd_global, cut, s->incumb_x, c, p->num);
		}

		//printf("s->subobj_est is : %f\n",s->subobj_est);

		/* Print all the cuts after adding a cut, for the purpose of cut index
		 checking. zl 

		 print_cut_info(c, p->num, "After adding a cut.");
		 */

#ifdef DEBUG
		printf("Changed incumbent cut to:\n");
		print_cut(c->cuts, p->num, s->incumb_cut);
#endif

	}
	return incumb_change;
}

void decrease_feacut_rownum(cell_type *c)
{
	int idx;
	for (idx = 0; idx < c->feasible_cuts_added->cnt; idx++)
		--c->feasible_cuts_added->val[idx]->row_num;
}

/***********************************************************************\
** This function updates all the structures necessary for forming
 ** a stochastic cut.  The latest observation of omega and the latest
 ** dual solution to the subproblem are added to their appropriate
 ** structures.  Then Pi x R and Pi x T are computed and for the latest
 ** omega and dual vector, and are added to the appropriate structures.
 **
 ** Note that the new column of delta is computed before a new row in
 ** lambda is calculated and before the new row in delta is completed,
 ** so that the intersection of the new row and new column in delta is
 ** only computed once (they overlap at the bottom, right-hand corner).
 \***********************************************************************/
BOOL stochastic_updates(sdglobal_type* sd_global, cell_type *c,
		lambda_type *lambda, sigma_type *sigma, delta_type *delta,
		omega_type *omega, num_type *num, sparse_vect *Rbar,
		sparse_matrix *Tbar, one_problem *subprob, vector Pi, int omeg_idx,
		BOOL new_omega)
{
	int lamb_idx;
	BOOL new_lamb = FALSE, new_sigma = FALSE;

#ifdef DEBUG
	int idx, obs;
#endif

#ifdef TRACE
	printf("\nInside stochastic_updates\n");
#endif

	/* Only need to calculate column if new observation of omega found */
	if (new_omega)
		calc_delta_col(sd_global, delta, lambda, omega, num, omeg_idx);

	/* Retrieve the dual solution from the latest subproblem */
	get_dual(Pi, subprob, num, num->sub_rows);

#ifdef RUN
	/*
	 print_vect(Pi, num->sub_rows, "Subproblem Dual, Pi");
	 */
#endif

#ifdef DEBUG
	printf("\n");
	print_vect(Pi, num->sub_rows, "Dual vector");
	printf("\n");
#endif

	lamb_idx = calc_lambda(sd_global, lambda, num, Pi, &new_lamb);

	// new_Mu???????????  added by Yifan to find out whether a new Mu has been generated!

	/* Do this afterwards because you need to have index to lambda */
	calc_sigma(sd_global, c, sigma, num, Pi, Rbar, Tbar, lamb_idx, new_lamb,
			&new_sigma);
	if (sd_global->MALLOC)
	{
		printf("After return from sigma\n");
		/* malloc_verify(); */
	}

	/* Only need to calculate row if a distinct lambda was found */
	/* We could use Pi, instead of lambda(Pi), for this calculation, */
	/* and save the time for expanding/reducing vector.  But for clarity... */

	/*Commented by Yifan: even though the lambda is the same, the current Pi might be a 
	 distinct one due to the variations in sigma*/
	if (new_lamb)
		calc_delta_row(sd_global, delta, lambda, omega, num, lamb_idx);

#ifdef DEBUG
	for(idx=0;idx<lambda->cnt;idx++)
	for(obs=0;obs<omega->most;obs++)
	if (valid_omega_idx(omega, obs))
	print_delta(delta, num, idx, obs);

	for(idx=0;idx<lambda->cnt;idx++)
	print_lambda(lambda, num, idx);

	for(obs=0;obs<omega->most;obs++)
	if (valid_omega_idx(omega, obs))
	print_omega(omega, num, obs);

	for(idx=0;idx<sigma->cnt;idx++)
	print_sigma(sigma, num, idx);
#endif

	return new_sigma;

}

/* Yifan 03/04/2012 Updated for Feasibility Cuts*/

/***********************************************************************\
** This function creates a new cut for the master program based on
 ** all of the observed outcomes of omega and on all the previous
 ** dual solutions to the subproblem.  It also uses some given
 ** solution to the master program (X).  For each outcome observed, the
 ** dual vector (Pi) is found which maximizes:  Pi x (R - T x X)
 ** The cut is then formed by calculating Pi x R and Pi x T for each
 ** of these maximal Pi's and averaging them.  If a given outcome
 ** has been repeated q times, then its maximal Pi is weighted
 ** by q in the average.  Note that the calculation will be divided
 ** into a determinisitic part, which need only be done once, and 
 ** a stochastic part, which must be computed for each observation.
 **
 ** Need to store Pi x Tbar x X for all Pi in a separate array ahead
 ** of time.  This way it isn't re-calculated for each omega...
 \***********************************************************************/
void SD_cut(sdglobal_type* sd_global, sigma_type *sigma, delta_type *delta,
		omega_type *omega, num_type *num, one_cut *cut, vector Xvect,
		double *pi_ratio, double max_ratio, double min_ratio, int num_samples,
		BOOL *dual_statble_flag)
{
	int c, cnt;
	int obs; /* Observation of omega being used */
	i_type istar; /* Index to optimizing Pi's */
	vector pi_Tbar_x; /* Array of PixTbarxX scalars for all Pi */
	BOOL pi_eval_flag = FALSE; /*TRUE for testing the impact of the new PI's */
	double *argmax_all; /*added by Yifan to calcuate argmax for new PI's and old PI's seperately*/
	double *argmax_new;
	double *argmax_old;
	double *beta;

	double argmax_dif_sum = 0;
	double argmax_all_sum = 0;
	double vari = 1.0;
	i_type istar_new;
	i_type istar_old;
	FILE *fptr;
#ifdef RECOURSE_OBJ
	FILE *subobj_ptr; /* by Yifan 02/02/12 */
#endif

#ifdef TRACE
	printf("Inside SD_cut\n");
#endif

	if (!(argmax_all = arr_alloc(1, double)))
		err_msg("Allocation", "SD_cut", "argmax_all");
	if (!(argmax_new = arr_alloc(1, double)))
		err_msg("Allocation", "SD_cut", "argmax_new");
	if (!(argmax_old = arr_alloc(1, double)))
		err_msg("Allocation", "SD_cut", "argmax_old");

	/* by Yifan 02/02/12 */
	if (cut->is_incumbent)
	{
		if (!(beta = arr_alloc(num->mast_cols+1, double))) /* Yifan 03/04/2012 Modified*/
			err_msg("Allocation", "SD_cut", "beta");
		if (!(cut->subobj_omega = arr_alloc(omega->most+1, double)))
			err_msg("Allocation", "SD_cut", "subobj_omega");
		if (!(cut->subobj_freq = arr_alloc(omega->most+1, int)))
			err_msg("Allocation", "SD_cut", "subobj_freq");

	}

	//added by Yifan 09/21/2011
	/*Calculate pi_eval_flag to determine the way of computing argmax*/
	if (num_samples > sd_global->config.PI_EVAL_START
			&& !(num_samples % sd_global->config.PI_CYCLE))
		pi_eval_flag = TRUE; //modified by Yifan for testing

	/* Need to store  Pi x Tbar x X independently of observation loop */
	if (!(pi_Tbar_x = arr_alloc(sigma->cnt, double)))
		err_msg("Allocation", "SD_cut", "pi_Tbar_x");

	/* Calculate (Pi x Tbar) x X by mult. each VxT by X, one at a time */
	for (cnt = 0; cnt < sigma->cnt; cnt++)
	{
		pi_Tbar_x[cnt] = 0;
		for (c = 1; c <= num->nz_cols; c++)
			pi_Tbar_x[cnt] += sigma->val[cnt].T[c] * Xvect[sigma->col[c]];
	}

	/* Assume the cut's fields were initialized to zero.  */

	/* Yifan 03/20/2012 Test for omega issues*/
	for (obs = 0; obs < omega->most; obs++)
		if (valid_omega_idx(omega, obs))
		{
			/* For each observation, find the Pi which maximizes height at X. */
			if (pi_eval_flag == TRUE)
			{
				//printf("\n This is iteration c->k: %d \n",num_samples);
				//istar = compute_istar(obs, cut, sigma, delta, Xvect, num, pi_Tbar_x, argmax_all, FALSE, num_samples);
				//printf("This is argmax OSD for obs %d : %f and istar(%d,%d)\n", obs, *argmax_all, istar.sigma, istar.delta);

				istar_old = compute_istar(obs, cut, sigma, delta, Xvect, num,
						pi_Tbar_x, argmax_old, pi_eval_flag, num_samples);
				istar_new = compute_new_istar(obs, cut, sigma, delta, Xvect,
						num, pi_Tbar_x, argmax_new, num_samples);
				if (*argmax_new > *argmax_old)
				{
					*argmax_all = *argmax_new;
					istar.sigma = istar_new.sigma;
					istar.delta = istar_new.delta;
				}
				else
				{
					*argmax_all = *argmax_old;
					istar.sigma = istar_old.sigma;
					istar.delta = istar_old.delta;
				}
				//printf("argmax:%f and finalistar(%d,%d) and ck %d\n",*argmax_all, istar.sigma, istar.delta, sigma->ck[istar.sigma]);
				/* modified by Yifan 2013.02.15 */

				/* argmax_dif_sum += *argmax_all*omega->weight[obs] - max(*argmax_old*omega->weight[obs],0);
				 
				 argmax_all_sum += *argmax_all*omega->weight[obs];*//* modified by Yifan 2013.05.12 */
				/* If Eta0=0, the above code works fine. But if Eta0<0, then we need to calculate the following way */
				argmax_dif_sum += max((*argmax_old-sd_global->Eta0),0)
						* omega->weight[obs];
				argmax_all_sum += max((*argmax_all-sd_global->Eta0),0)
						* omega->weight[obs];

				// printf("This is argmax NSD for obs %d : %f and istar(%d,%d) and sigma from %d iteration \n", obs, *argmax_all, istar.sigma, istar.delta, sigma->ck[istar.sigma]);
			}
			else
				istar = compute_istar(obs, cut, sigma, delta, Xvect, num,
						pi_Tbar_x, argmax_all, pi_eval_flag, num_samples);

			cut->istar[obs] = istar.sigma;

			/* by Yifan 02/02/12 */
			if (cut->is_incumbent)
			{
				cut->subobj_freq[obs] = omega->weight[obs];

				cut->subobj_omega[obs] = sigma->val[istar.sigma].R
						+ delta->val[istar.delta][obs].R;

				for (c = 1; c <= num->nz_cols; c++)
					beta[sigma->col[c]] = sigma->val[istar.sigma].T[c];

				for (c = 1; c <= num->rv_cols; c++)
					beta[delta->col[c]] += delta->val[istar.delta][obs].T[c];

				cut->subobj_omega[obs] -= CxX(beta, Xvect, num->mast_cols);

			}

			// printf("c->k=%d  STABLE_FLAG: %d, max: %f, min: %f, diff: %f, This is new pi's impact ratio # %d: %f \n", num_samples, *dual_statble_flag,max_ratio, min_ratio, max_ratio-min_ratio, num_samples %sd_global->config.SCAN_LEN,  pi_ratio[num_samples %sd_global->config.SCAN_LEN]);

#ifdef LOOP
			printf("\nistar.sigma=%d.  istar.delta=%d", istar.sigma, istar.delta);
#endif

			/* Average using these Pi's to calculate the cut itself. */
			cut->alpha += sigma->val[istar.sigma].R * omega->weight[obs];
			cut->alpha += delta->val[istar.delta][obs].R * omega->weight[obs];

			for (c = 1; c <= num->nz_cols; c++)
				cut->beta[sigma->col[c]] += sigma->val[istar.sigma].T[c]
						* omega->weight[obs];
			for (c = 1; c <= num->rv_cols; c++)
				cut->beta[delta->col[c]] += delta->val[istar.delta][obs].T[c]
						* omega->weight[obs];

		}

	if (pi_eval_flag == TRUE)
	{
		pi_ratio[num_samples % sd_global->config.SCAN_LEN] = argmax_dif_sum
				/ argmax_all_sum;
		if (num_samples - sd_global->config.PI_EVAL_START
				> sd_global->config.SCAN_LEN)
			vari = calc_var(sd_global, pi_ratio, NULL, NULL, 0); /*added by Yifan return vari*/
		if (DBL_ABS(vari) >= .000002
				|| (pi_ratio[num_samples % sd_global->config.SCAN_LEN]) < 0.95)
			*dual_statble_flag = FALSE;
		else
			*dual_statble_flag = TRUE;
	}

	if (0)
	{
		/*printf("c->k=%d  STABLE_FLAG: %d, max: %f, min: %f, diff: %f, This is new pi's impact ratio # %d: %f \n", num_samples, *dual_statble_flag,max_ratio, min_ratio, max_ratio-min_ratio, num_samples %sd_global->config.SCAN_LEN,  pi_ratio[num_samples %sd_global->config.SCAN_LEN]);*/
		fptr = fopen("pi_ratio.log", "a");
		fprintf(fptr, "c->k=%d, %d:%f, vari:%.9f, dual_stable_flag: %d\n",
				num_samples, num_samples % sd_global->config.SCAN_LEN,
				pi_ratio[num_samples % sd_global->config.SCAN_LEN], vari,
				*dual_statble_flag);
		fclose(fptr);
	}

	/*added by Yifan 02/02/12 output subproblem objective with a fixed incumbent_x and different omegas*/
#ifdef RECOURSE_OBJ
	if (cut->is_incumbent)
	{
		subobj_ptr = fopen("subobj_omega.txt","a");
		fprintf(subobj_ptr,"New incumbent found.\n");
		for (obs=0; obs<omega->most; obs++)
		{
			fprintf(subobj_ptr,"%f\t%d\n",cut->subobj_omega[obs],cut->subobj_freq[obs]);
		}
		fclose(subobj_ptr);
	}
#endif

	cut->alpha /= num_samples;
	for (c = 1; c <= num->mast_cols; c++)
		cut->beta[c] /= num_samples;

	mem_free(pi_Tbar_x);
	mem_free(argmax_all);
	mem_free(argmax_new);
	mem_free(argmax_old);
	if (cut->is_incumbent)
	{
		mem_free(beta);
	}
#ifdef TRACE
	printf("Exiting SD_cut\n");
#endif
}

/***************************************************************************\
 ** This function adds new feasibility cuts. It first adds feasibility cuts from old pi's
 ** associated with the new omega generated. Cuts from a new dual extreme ray(new pi) and
 ** all omegas generated so far are added to the feasible_cuts_pool structure afterwards.
 \***************************************************************************/

int FEA_cut(sdglobal_type* sd_global, cell_type *cell, soln_type *soln,
		sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num,
		int num_samples, BOOL *dual_statble_flag, BOOL new_omega,
		BOOL new_sigma)
{
	int c, cut_cnt, sig_pi, del_pi;
	int obs, last_omega;
	double alpha;
	double *beta;
	last_omega = omega->most - 1; /* The most recent omega generated */
	cut_cnt = cell->feasible_cuts_pool->cnt;

	if (!(beta = arr_alloc(num->mast_cols+1, double)))
		err_msg("Allocation", "SD_cut", "beta");

	/* Adds feasibility cuts using old Pi's */
	if (new_omega)
	{
		for (sig_pi = 0; sig_pi < sigma->cnt - 1; sig_pi++)
		{
			del_pi = sigma->lamb[sig_pi];
			alpha = 0.0;
			for (c = 0; c <= num->mast_cols; c++)
				beta[c] = 0.0;
			alpha = sigma->val[sig_pi].R + delta->val[del_pi][last_omega].R;
			for (c = 1; c <= num->nz_cols; c++)
				beta[sigma->col[c]] += sigma->val[sig_pi].T[c];
			for (c = 1; c <= num->rv_cols; c++)
				beta[delta->col[c]] += delta->val[del_pi][last_omega].T[c];
			cut_cnt = add_to_cutpool(sd_global, &alpha, beta, cell, soln,
					num->mast_cols);
		}
	}
	/* Adds feasibility cuts using the new Pi */
	if (new_sigma)
	{
		for (obs = 0; obs < omega->most; obs++)
		{
			if (valid_omega_idx(omega, obs))
			{
				sig_pi = sigma->cnt - 1;
				del_pi = sigma->lamb[sig_pi];
				alpha = 0.0;
				for (c = 0; c <= num->mast_cols; c++)
					beta[c] = 0.0;
				alpha = sigma->val[sig_pi].R + delta->val[del_pi][obs].R;
				for (c = 1; c <= num->nz_cols; c++)
					beta[sigma->col[c]] += sigma->val[sig_pi].T[c];
				for (c = 1; c <= num->rv_cols; c++)
					beta[delta->col[c]] += delta->val[del_pi][obs].T[c];
				cut_cnt = add_to_cutpool(sd_global, &alpha, beta, cell, soln,
						num->mast_cols);
			}
		}
	}

	if (new_omega && !new_sigma)
	{
		/* add cut with data located at buttom right corner */
		sig_pi = sigma->cnt - 1;
		del_pi = sigma->lamb[sig_pi];
		alpha = 0.0;
		for (c = 0; c <= num->mast_cols; c++)
			beta[c] = 0.0;
		del_pi = sigma->lamb[sig_pi];
		alpha = sigma->val[sig_pi].R + delta->val[del_pi][last_omega].R;
		for (c = 1; c <= num->nz_cols; c++)
			beta[sigma->col[c]] += sigma->val[sig_pi].T[c];
		for (c = 1; c <= num->rv_cols; c++)
			beta[delta->col[c]] += delta->val[del_pi][last_omega].T[c];
		cut_cnt = add_to_cutpool(sd_global, &alpha, beta, cell, soln,
				num->mast_cols);
	}mem_free(beta);
#ifdef TRACE
	printf("Exiting FEA_cut\n");
#endif
	return cut_cnt;

}

/***********************************************************************\
 ** This function add a new feasibility cut to the cut pool using
 ** alpha and beta provided.
 \***********************************************************************/
int add_to_cutpool(sdglobal_type* sd_global, double *alpha, double *beta,
		cell_type *c, soln_type *s, int mast_cols)
{
	one_cut *cut;
	double cut_alpha;
	int cnt;
#ifdef TRACE
	printf("Inside add_to_cutpool\n");
#endif
	for (cnt = 0; cnt < c->feasible_cuts_pool->cnt; cnt++)
	{
		cut_alpha = c->feasible_cuts_pool->val[cnt]->alpha;
		if (DBL_ABS(*alpha-cut_alpha) < sd_global->config.TOLERANCE)
		{
			if (equal_arr(beta, c->feasible_cuts_pool->val[cnt]->beta,
					mast_cols, sd_global->config.TOLERANCE))
			{
				return c->feasible_cuts_pool->cnt;
			}
		}
	}
	cut = (one_cut *) mem_malloc (sizeof(one_cut));
	cut->cut_obs = c->k;
	cut->omega_cnt = s->omega->most;
	cut->slack_cnt = 0;
	cut->is_incumbent = FALSE;

	if (!(cut->istar = arr_alloc(s->omega->most, int)))
		err_msg("Allocation", "add_to_cutpool", "istar");

	if (!(cut->beta = arr_alloc(mast_cols+1, double)))
		err_msg("Allocation", "add_to_cutpool", "beta");

	cut->subfeaflag = FALSE;

	cut->alpha = *alpha;

	for (cnt = 0; cnt <= mast_cols; cnt++)
		cut->beta[cnt] = beta[cnt];

	/* Yifan 03/04/2012 The following is commented out by Yifan*/
	/* One_norm might cuase the same vector become different*/
	/* cut->beta[0] = one_norm(cut->beta, mast_cols); */

	c->feasible_cuts_pool->val[c->feasible_cuts_pool->cnt] = cut;

#ifdef TRACE
	printf("Exiting add_to_cutpool\n");
#endif
	return ++c->feasible_cuts_pool->cnt;
}

/***********************************************************************\
 ** Feasibility cuts which are violated by the current candidate solution
 ** are added to the c->feasibile_cuts_added structure.
 ** If Beta x X  >=  Alpha can not be satisfied, then it  will be added
 ** as a new feasibility cut to the master problem temporarily in
 ** c->feasibile_cuts_added structure.

 \***********************************************************************/
int FEA_cut_check_add(sdglobal_type* sd_global, cell_type *cell, prob_type *p,
		soln_type *soln, vector x_k)
{
	int cnt;
	int idx, c;
	double Beta_x, Alpha;
	BOOL duplic_cut;

	for (idx = 0; idx < cell->feasible_cuts_pool->cnt; idx++)
	{
		duplic_cut = FALSE;
		Alpha = cell->feasible_cuts_pool->val[idx]->alpha;
		for (c = 0; c < cell->feasible_cuts_added->cnt; c++)
		{
			if (DBL_ABS(Alpha - cell->feasible_cuts_added->val[c]->alpha)
					< sd_global->config.TOLERANCE)
			{
				if (equal_arr(cell->feasible_cuts_pool->val[idx]->beta,
						cell->feasible_cuts_added->val[c]->beta,
						p->num->mast_cols, sd_global->config.TOLERANCE))
				{
					duplic_cut = TRUE;
					break;
				}
			}
		}

		/* Yifan 03/23/2012 Add those cuts in cut pool that will be violated by incumbent solution*/
		/* Check if the cut will be violated by the incumbent solution*/
		Beta_x = CxX(cell->feasible_cuts_pool->val[idx]->beta, soln->incumb_x,
				p->num->mast_cols);
		/* If Beta_x is extremetly close to Alpha, we need to detect that. Yifan 04/26/2013 */
		//if (Beta_x < Alpha || DBL_ABS(Beta_x-Alpha) <sd_global->config.TOLERANCE) {
		if (Beta_x < Alpha)
		{
			cell->incumb_infea = TRUE;
			if (duplic_cut == TRUE)
			{
				printf(
						"Incumbent violates one old cut from fea cuts pool (this cut also exists in fea_cut_added)\n");
				continue;
			}
			else
			{
				printf(
						"Incumbent violates one new cut from fea cuts pool (this cut is not in fea_cut_added but will be)");
				print_cut(cell->feasible_cuts_pool, p->num, idx);
			}
			printf("\n");
			cell->feasible_cuts_added->val[cell->feasible_cuts_added->cnt] =
					cell->feasible_cuts_pool->val[idx];
			cell->feasible_cuts_added->cnt++;
			printf("cut added to master...due to Incumbent violation");
			print_cut(cell->feasible_cuts_pool, p->num, idx);
			printf("\n");
		}
		else
		{
			/* Check if the cut will be violated by the candidate solution*/
			if (duplic_cut == TRUE)
				continue;
			Beta_x = CxX(cell->feasible_cuts_pool->val[idx]->beta, x_k,
					p->num->mast_cols);

			/* If Beta_x is extremetly close to Alpha, we need to detect that. Yifan 04/26/2013 */
			if (Beta_x < Alpha)
			{
				printf(
						"Candidate violates one cut from fea cuts pool (this cut is not in fea_cut_added but will be)");
				cell->feasible_cuts_added->val[cell->feasible_cuts_added->cnt] =
						cell->feasible_cuts_pool->val[idx];
				cell->feasible_cuts_added->cnt++;
				printf("cuts added to master...due to Candidate violation");
				print_cut(cell->feasible_cuts_pool, p->num, idx);
				printf("\n");
			}
		}

	}

	cnt = cell->feasible_cuts_added->cnt;
	return cnt;
}

/***********************************************************************\
** This function loops through all the dual vectors found so far
 ** and returns the index of the one which satisfies the expression:
 **
 **	argmax { Pi x (R - T x X) | all Pi }
 **
 ** where X, R, and T are given.  It is calculated in this form:
 **
 ** 	Pi x Rbar + Pi x Romega + (Pi x Tbar) x X + (Pi x Tomega) x X
 **
 ** Since the Pi's are stored in two different structures (sigma and delta),
 ** the index to the maximizing Pi is actually a structure containing
 ** two indices.  (While both indices point to pieces of the dual vectors,
 ** sigma and delta may not be in sync with one another due to elimination
 ** of non-distinct or redundant vectors.
 \***********************************************************************/
i_type compute_istar(int obs, one_cut *cut, sigma_type *sigma,
		delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
		double *argmax, BOOL pi_eval, int ictr)
{
	double arg;
//  double	argmax;             //modified by Yifan to return argmax value 09/22/2011
	int sig_pi, del_pi;
	int c, new_pisz;
	i_type ans;

#ifdef LOOP
	printf("Inside compute_istar\n");
#endif
	//printf("\n***IN COMPUTE ISTAR ...\n");
	if (pi_eval == TRUE)
	{
		//new_pisz =sd_global->config.PI_EVAL_START;
		new_pisz = ictr / 10 + 1;
	}
	else
		new_pisz = 0;

	ictr -= new_pisz; /*evaluate the pi's generated in the first 90% iterations */

	*argmax = -DBL_MAX;

	/*added by Yifan to enable parallel process*/
// #pragma omp for private(sig_pi,del_pi, arg)
	for (sig_pi = 0; sig_pi < sigma->cnt; sig_pi++)
	{
		if (sigma->ck[sig_pi] <= ictr)
		{
			/* Find the row in delta corresponding to this row in sigma */
			del_pi = sigma->lamb[sig_pi];

			/* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
			arg = sigma->val[sig_pi].R + delta->val[del_pi][obs].R
					- Pi_Tbar_X[sig_pi];

			/* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */

			for (c = 1; c <= num->rv_cols; c++)
				arg -= delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];

#ifdef LOOP
			print_sigma(sigma, num, sig_pi);
			print_delta(delta, num, del_pi, obs);
			printf("\nResulting arg=%f", arg);
#endif

			if (arg > (*argmax))
			{
				*argmax = arg;
				ans.sigma = sig_pi;
				ans.delta = del_pi;
				//printf("argmax:%f and istar(%d,%d)\n", arg, ans.sigma, ans.delta);
			}
		}
	}

	/* 
	 if (pi_eval==FALSE) {
	 printf("argmax:%f and istar(%d,%d) and ck %d\n", *argmax, ans.sigma, ans.delta, sigma->ck[ans.sigma]);
	 }
	 else
	 printf("argmax:%f and oldistar(%d,%d) and ck %d\n", *argmax, ans.sigma, ans.delta, sigma->ck[ans.sigma]);
	 */
#ifdef LOOP
	printf("Exiting compute_istar\n");
#endif

	return ans;
}

i_type compute_new_istar(int obs, one_cut *cut, sigma_type *sigma,
		delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
		double *argmax, int ictr)
{
	double arg;
	int sig_pi, del_pi;
	int c, new_pisz;
	i_type ans;
	ans.sigma = 0;
	ans.delta = 0;

#ifdef LOOP
	printf("Inside compute_new_istar\n");
#endif

	new_pisz = ictr / 10 + 1;
	ictr -= new_pisz; /*evaluate the pi's generated in the last 10% iterations */

	*argmax = -DBL_MAX;
	for (sig_pi = 0; sig_pi < sigma->cnt; sig_pi++)
	{
		if (sigma->ck[sig_pi] > ictr)
		{
			/* Find the row in delta corresponding to this row in sigma */
			del_pi = sigma->lamb[sig_pi];

			/* Start with (Pi x Rbar) + (Pi x Romega) + (Pi x Tbar) x X */
			arg = sigma->val[sig_pi].R + delta->val[del_pi][obs].R
					- Pi_Tbar_X[sig_pi];

			/* Subtract (Pi x Tomega) x X. Multiply only non-zero VxT values */
			for (c = 1; c <= num->rv_cols; c++)
				arg -= delta->val[del_pi][obs].T[c] * Xvect[delta->col[c]];

#ifdef LOOP
			print_sigma(sigma, num, sig_pi);
			print_delta(delta, num, del_pi, obs);
			printf("\nResulting arg=%f", arg);
#endif

			if (arg > (*argmax))
			{
				*argmax = arg;
				ans.sigma = sig_pi;
				ans.delta = del_pi;

			}
		}
	}

//printf("argmax:%f and newistar(%d,%d) and ck %d\n",*argmax, ans.sigma, ans.delta, sigma->ck[ans.sigma]);

#ifdef LOOP
	printf("Exiting compute_new_istar\n");
#endif
	return ans;
}

/***********************************************************************\
** This function will remove the oldest cut whose corresponding dual
 ** variable is zero (thus, a cut which was slack in last solution).
 \***********************************************************************/
void reduce_cuts(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	int oldest_cut, min_cut_obs;
	int idx;
	double *dual;
	double height, min_height;

#ifdef TRACE
	printf("Inside reduce_cuts\n");
#endif

    min_cut_obs = c->k;
	oldest_cut = c->cuts->cnt;
    dual = arr_alloc(p->num->mast_rows+c->cuts->cnt+1, double);
    
    /* modified by Yifan 2013.9.16  Later, we need to update SDLP as SDIP */
    if (sd_global->config.MASTER_TYPE == SDQP) {
        /* Original rows and all the added cuts.  The eta row isn't needed */
        get_dual(dual, c->master, p->num, p->num->mast_rows + c->cuts->cnt);
        
        
        for (idx = 0; idx < c->cuts->cnt; idx++)
        {
            /* avoid dropping incumbent cut*/
            if (idx == s->incumb_cut)
            {
                continue;
            }
            
            if (c->cuts->val[idx]->cut_obs < min_cut_obs
				&& DBL_ABS(dual[c->cuts->val[idx]->row_num + 1])
                <= sd_global->config.TOLERANCE)
            {
                min_cut_obs = c->cuts->val[idx]->cut_obs;
                oldest_cut = idx;
            }
        }
    }
	

    /* If dual variables cannot be used to decide whcih cut to drop
     then we check each cut's height and drop the most loose cut */
	if (oldest_cut == c->cuts->cnt)
	{
		min_height = cut_height(sd_global, c->cuts->val[0], s->candid_x, c,
				p->num);
		oldest_cut = 0;
		for (idx = 1; idx < c->cuts->cnt; idx++)
		{
			if (idx == s->incumb_cut && sd_global->config.MASTER_TYPE == SDQP)
			{
				continue;
			}
            if (sd_global->config.MASTER_TYPE == SDQP) {
                height = cut_height(sd_global, c->cuts->val[idx], s->candid_x, c,
                                    p->num);
            }
            else{
                height = cut_height(sd_global, c->cuts->val[idx], s->mip_x, c,
                                    p->num);
            }
			
			if (height < min_height)
			{
				min_height = height;
				oldest_cut = idx;
			}
		}
		/* Modified by Yifan to drop the loose cut except incumbent */
		//err_msg("No_cut_dropped", "reduce_cuts", "found");
	}

#ifdef DEBUG
	printf("Dropped cut %d as the oldest\n", c->cuts->val[oldest_cut]->row_num);
#endif

	if (c->cuts->val[oldest_cut]->subfeaflag == TRUE)
	{
		drop_cut(oldest_cut, p, c, s);
		decrease_feacut_rownum(c);
	}

	mem_free(dual);
}

/***********************************************************************\
** When invoked, this function will drop all cuts which have been slack
 ** in more than a given number of consecutive solutions of the master problem.
 ** It first obtains a dual solution to the most recent master, and 
 ** increments the slack_cnt of each cut whose dual variable is zero.
 ** It then scans the list of cuts for potential cuts to drop.
 \***********************************************************************/
void thin_cuts(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	int cnt;
	double *pi, *dual;

#ifdef TRACE
	printf("Inside thin_cuts\n");
#endif

	/* 
	 ** Get the dual solution to master, particularly for cut constraints 
	 */
	if (!(dual = arr_alloc(p->num->mast_rows+c->cuts->cnt+1, double)))
		err_msg("Allocation", "thin_cuts", "dual");
	get_dual(dual, c->master, p->num,
			p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt);
	pi = dual + p->num->mast_rows;

	/* 
	 ** For every cut, check its dual.  If it is zero, increment slack_cnt
	 ** and check whether it has exceeded the limit.  If so, drop it.
	 */
	for (cnt = c->cuts->cnt - 1; cnt >= 0; cnt--)
		if (DBL_ABS(pi[cnt]) <= sd_global->config.TOLERANCE
				&& cnt != s->incumb_cut
				&& c->cuts->val[cnt]->subfeaflag == TRUE)
			if ((++c->cuts->val[cnt]->slack_cnt) >= sd_global->config.DROP_TIME)
			{
#ifdef RUN
				printf("\n**Dropped cut %d for being slack\n\n", cnt);
#endif
				drop_cut(cnt, p, c, s);
			}
			else
				/* nothing */;
		else
			c->cuts->val[cnt]->slack_cnt = 0;

	/* What if the incumbent cut has been slack for DROP_TIME iterations ? */
	/* It better not be... re-form it at least! */
	/* Here, you reset the slack_cnt of the incumbent to zero every time */

	mem_free(dual);
}

/***********************************************************************\
** This function removes a cut from both the cut_type structure and the
 ** master problem constraint matrix.  In the cuts->val array, the last
 ** cut is swapped into the place of the exiting cut.  In the constraint
 ** matrix, the row is deleted, and the row numbers of all constraints
 ** below it are decremented.  Note that it is possible for the incumbent
 ** cut to be swapped from its current position to that of the dropped
 ** cut.  In this case, we must update s->incumb_cut appropriately.
 ** 
 ** Cuts can be added and removed from the cuts structure only via
 ** this function and add_cut below.
 **
 ** Should there be a corresponding add_cut, so that the only way to get
 ** a new cut is through add_cut, and the only way to get rid of one 
 ** is through drop_cut... control of data...
 \***********************************************************************/
void drop_cut(int cut_idx, prob_type *p, cell_type *c, soln_type *s)
{
	int row; /* Row number in constraint matrix */
	int idx;

#ifdef TRACE
	printf("Inside drop_cut, cut_idx=%d\n", cut_idx);
#endif

	/* Print all the cuts before droping a cut, for the purpose of cut index
	 checking. zl 
	 
	 print_cut_info(c, p->num, "Before Dropping a cut");
	 */

	/* Get rid of the old cut */
	row = c->cuts->val[cut_idx]->row_num;
	if (!remove_row(c->master, row))
	{
		print_contents(c->master, "contents.out");
		err_msg("remove_row", "drop_cut", "returned FALSE");
	}
	free_cut(c->cuts->val[cut_idx]);

#ifdef DEBUG
	printf("cuts->cnt=%d, row=%d, s->incumb_cut=%d\n",
			c->cuts->cnt, row, s->incumb_cut);
#endif

	/* Update the surviving cuts */
	c->cuts->val[cut_idx] = c->cuts->val[--c->cuts->cnt];
	for (idx = 0; idx < c->cuts->cnt; idx++)
		if (c->cuts->val[idx]->row_num > row)
			--c->cuts->val[idx]->row_num;

	/* Worry about swapping down the incumbent cut */
	if (s->incumb_cut == c->cuts->cnt)
		s->incumb_cut = cut_idx;

	/* What if the incumbent cut was the one that was swapped? */
	/* But it might be the incumbent cut being dropped ! */

	/* Print all the cuts after dropping a cut, for the purpose of cut index
	 checking. zl 

	 print_cut_info(c, p->num, "After dropping a cut");
	 */
}

/***********************************************************************\
** This function will add a new cut to the master problem.  If 
 ** This function adds a cut to the master problem, in the form
 ** of a constraint.  When the objective function is: Min cx + eta
 ** a cut is specified as:   eta >= Alpha + Beta x X.
 ** In this function, it is entered in the master problem in this form:
 **
 **      Beta x X  +  Eta  >=  Alpha
 **
 ** For now, Eta does not have a coefficient.  All the eta
 ** coefficients (not just the one for this cut) will be updated in
 ** change_eta_col(), before each master program is solved.
 ** The function returns TRUE if the cut was successfully added;
 ** FALSE otherwise, in case anyone cares.
 **
 ** Doesn't work for cells!!! cuts->cnt may be larger than p->num->max_cuts
 ** if there were many ancestor cuts!!!! 
 \***********************************************************************/
int add_cut(sdglobal_type* sd_global, one_cut *cut, prob_type *p, cell_type *c,
		soln_type *s)
{
	int beg_col; /* column where beta coefficients begin */
	int end_col; /* column where beta coefficients finish */
	int *coef_col; /* column number of each beta coefficient */
	int cnt;
	double rhs; /* rhs value in regularized QP method. */
//	BOOL reach_max_cuts = TRUE;
	int row;
	/*int idx;*/
#ifdef TRACE
	printf("Inside add_cut\n");
#endif

	/* Print all the cuts before adding a cut, for the purpose of cut index
	 checking. zl 

	 print_cut_info(c, p->num, "Before adding a cut");
	 */

	/* Insure that there is room to add another cut */
	if (c->cuts->cnt >= p->num->max_cuts && sd_global->config.MASTER_TYPE == SDQP)
		reduce_cuts(sd_global, p, c, s);
//	else
//		reach_max_cuts = FALSE;

	/*deleate all feasibility cuts so that the new optimality cut can be added to the end of the optimality cuts  Yifan /03/09/2012*/
	for (cnt = c->feasible_cuts_added->cnt - 1; cnt >= 0; cnt--)
	{
		row = c->feasible_cuts_added->val[cnt]->row_num;
		if (!remove_row(c->master, row))
		{
			print_problem(c->master, "master_errors_notice");
			print_contents(c->master, "contents.out");
			err_msg("remove_row", "drop_cut", "returned FALSE");
		}
	}

	/*
	 ** Initialize an array to specify columns of each coefficient in
	 ** beta.  The one-norm of beta is temporarily used as the coefficeint on
	 ** eta (it is assumed to be replaced in the next step in solve_master()).
	 */
	beg_col = 0;
	end_col = p->num->mast_cols;

	/*
	 ** Assign the column values for each element in the beta vector.
	 ** The zeroth position will be set to the eta column, which will
	 ** temporarily receive a coefficient equal to the 1-norm of beta.
	 ** (or whatever is stored in beta[0])
	 **
	 ** This ought to be reworked... shouldn't have to alloc every time...
	 */
	/*
	 if (cut->subfeaflag==FALSE) {
	 for (idx=0; idx<c->cuts->cnt; idx++) {
	 if (idx<c->cuts->val[idx]->subfeaflag==FALSE) {
	 if (DBL_ABS(cut->alpha-c->cuts->val[idx]->alpha)<sd_global->config.TOLERANCE) {
	 if (equal_arr(cut->beta, c->cuts->val[idx]->beta, p->num->mast_cols,sd_global->config.TOLERANCE)) {
	 new_cut = FALSE;
	 }
	 }
	 }
	 }
	 }
	 */

	/*added by Yifan 02/27/2012*/

	if (!(coef_col = arr_alloc(p->num->mast_cols+1, int)))
		err_msg("Allocation", "add_cut", "coef_col");
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
		coef_col[cnt + 1] = cnt;
	coef_col[0] = p->num->mast_cols;

#ifdef DEBUG
	printf("Adding the row:\n");
	print_vect(cut->beta, end_col - beg_col, "c->beta");
#endif

	/*
	 ** Add the cut (it's a ">=" constraint) to the master, with coefficeints 
	 ** as specified in beta, and right hand side as specified by alpha.
	 */

	/* In the regularized QP method, we need to shift the rhs of the cut from
	 'x' to 'd' each time we add a cut. (We do not need to worry about it
	 when dropping a cut.) That is, in regularized QP method, the rhs will
	 become    
	 alpha - beta * incumb_x
	 instead of
	 alpha
	 as in the LP method.
	 */

	if (sd_global->config.MASTER_TYPE == SDLP)
		rhs = cut->alpha;
	else
		rhs = cut->alpha - CxX(cut->beta, s->incumb_x, p->num->mast_cols);

	cut->alpha_incumb = rhs;

	/* print out information for checking purpose. zl */
	/*  fprintf (g_FilePointer, "\n***** In cuts.c -- add_cut() *****\n");
	 fprintf (g_FilePointer, "mast_cols = %d\n", p->num->mast_cols);
	 fprintf (g_FilePointer, "alpha = %f,  rhs = %f.\n", cut->alpha, rhs);
	 for (cnt=0; cnt<p->num->mast_cols; cnt++)
	 fprintf (g_FilePointer, "inc_x[%d] = %f,  beta[%d] = %f\n", 
	 cnt+1, s->incumb_x[cnt+1], cnt+1, cut->beta[cnt+1]); 
	 */
	if (!add_row(c->master, beg_col, end_col, coef_col, cut->beta, GE, rhs))
		err_msg("LP solver", "add_cut", "ans");

	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	cut->row_num = p->num->mast_rows + c->cuts->cnt;

	/* Since a new optimality cut is added, the row number of every feasibility cut will increase by 1*/
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		++c->feasible_cuts_added->val[cnt]->row_num;
	}

	/* add all the removed feasibility cuts back to the problem after optimality cuts*/
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		if (sd_global->config.MASTER_TYPE == SDLP)
			rhs = c->feasible_cuts_added->val[cnt]->alpha;
		else
			rhs = c->feasible_cuts_added->val[cnt]->alpha
					- CxX(c->feasible_cuts_added->val[cnt]->beta, s->incumb_x,
							p->num->mast_cols) + sd_global->config.FEA_TOLER;
		/* putshing the solution a little inside the boundary 04/30/2013 Yifan*/

		if (!add_row(c->master, beg_col, end_col, coef_col,
				c->feasible_cuts_added->val[cnt]->beta, GE, rhs))
			err_msg("LP solver", "add_cut_to_master", "ans");
	}

	mem_free(coef_col);

	c->cuts->val[c->cuts->cnt] = cut;

	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl  
	 
	 print_cut_info(c, p->num, "After adding a cut");
	 */

	return c->cuts->cnt++;
}

void add_cut_to_master(sdglobal_type* sd_global, one_cut *cut, prob_type *p,
		cell_type *c, soln_type *s, int idx)
{
	int beg_col; /* column where beta coefficients begin */
	int end_col; /* column where beta coefficients finish */
	int *coef_col; /* column number of each beta coefficient */
	int cnt;
	double rhs; /* rhs value in regularized QP method. */

#ifdef TRACE
	printf("Inside add_cut_to_master\n");
#endif

	/*
	 ** Initialize an array to specify columns of each coefficient in
	 ** beta.  The one-norm of beta is temporarily used as the coefficeint on
	 ** eta (it is assumed to be replaced in the next step in solve_master()).
	 */
	beg_col = 0;
	end_col = p->num->mast_cols;

	/*
	 ** Assign the column values for each element in the beta vector.
	 ** The zeroth position will be set to the eta column, which will
	 ** temporarily receive a coefficient equal to the 1-norm of beta.
	 ** (or whatever is stored in beta[0])
	 **
	 ** This ought to be reworked... shouldn't have to alloc every time...
	 */
	if (!(coef_col = arr_alloc(p->num->mast_cols+1, int)))
		err_msg("Allocation", "add_cut", "coef_col");
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
		coef_col[cnt + 1] = cnt;
	coef_col[0] = p->num->mast_cols;

#ifdef DEBUG
	printf("Adding the row:\n");
	print_vect(cut->beta, end_col - beg_col, "c->beta");
#endif

	/*
	 ** Add the cut (it's a ">=" constraint) to the master, with coefficeints 
	 ** as specified in beta, and right hand side as specified by alpha.
	 */

	/* In the regularized QP method, we need to shift the rhs of the cut from
	 'x' to 'd' each time we add a cut. (We do not need to worry about it
	 when dropping a cut.) That is, in regularized QP method, the rhs will
	 become alpha - beta * incumb_x instead of alpha as in the LP method.
	 */

	if (sd_global->config.MASTER_TYPE == SDLP)
		rhs = cut->alpha;
	else
		rhs = cut->alpha - CxX(cut->beta, s->incumb_x, p->num->mast_cols)
				+ sd_global->config.FEA_TOLER;
	/* putshing the solution a little inside the boundary 04/30/2013 Yifan*/

	/* print out information for checking purpose. zl */
	/*  fprintf (g_FilePointer, "\n***** In cuts.c -- add_cut() *****\n");
	 fprintf (g_FilePointer, "mast_cols = %d\n", p->num->mast_cols);
	 fprintf (g_FilePointer, "alpha = %f,  rhs = %f.\n", cut->alpha, rhs);
	 for (cnt=0; cnt<p->num->mast_cols; cnt++)
	 fprintf (g_FilePointer, "inc_x[%d] = %f,  beta[%d] = %f\n", 
	 cnt+1, s->incumb_x[cnt+1], cnt+1, cut->beta[cnt+1]); 
	 */

	if (!add_row_to_master(c->master, beg_col, end_col, coef_col, cut->beta, GE,
			rhs))
		err_msg("LP solver", "add_cut", "ans");
	/* THIS IS A SERIOUS PROBLEM, ROW_NUM IS NOT CORRECT, ONLY ONE CUT IN MASTER!!!*/
	/* Yifan 03/11/2012 Be careful of this row number*/
	cut->row_num = p->num->mast_rows + c->cuts->cnt + idx;

	mem_free(coef_col);

	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl  
	 
	 print_cut_info(c, p->num, "After adding a cut");
	 */
#ifdef TRACE
	printf("Exiting add_cut_to_master\n");
#endif

}

/***********************************************************************\
** This function allocates memory for the arrays inside a single
 ** cut, and initializes its values accordingly.  The cut structure 
 ** itself is assumed to be already allocated.  Note, each beta 
 ** vector contains room for its one-norm, thought it just gets
 ** filled with zero anyway.
 \***********************************************************************/
one_cut *new_cut(int num_x, int num_istar, int num_samples)
{
	one_cut *cut;
	int cnt;

#ifdef TRACE
	printf("Inside new_cut\n");
#endif

	cut = (one_cut *) mem_malloc (sizeof(one_cut));
	cut->cut_obs = num_samples;
	cut->omega_cnt = num_istar;
	cut->slack_cnt = 0;
	cut->is_incumbent = FALSE; /*added by Yifan 02/02/2012 new cut is by defalut not incumbent*/

	if (!(cut->istar = arr_alloc(num_istar, int)))
		err_msg("Allocation", "new_cut", "istar");

	if (!(cut->beta = arr_alloc(num_x+1, double)))
		err_msg("Allocation", "new_cut", "beta");

	cut->subfeaflag = TRUE;

	cut->alpha = 0.0;
	cut->alpha_incumb = 0.0;

	for (cnt = 0; cnt <= num_x; cnt++)
		cut->beta[cnt] = 0.0;

	return cut;
}

/*added by Yifan to generate the feasibility cut*/
one_cut *new_fea_cut(int num_x, int num_istar, int num_samples)
{
	one_cut *cut;
	int cnt;

#ifdef TRACE
	printf("Inside new_cut\n");
#endif

	cut = (one_cut *) mem_malloc (sizeof(one_cut));
	cut->cut_obs = num_samples;
	cut->omega_cnt = num_istar;
	cut->slack_cnt = 0; /*make sure the cut won't be dropped Yifan /08/22/2011*/
	cut->is_incumbent = FALSE; /*added by Yifan 02/02/2012 new cut is by defalut not incumbent*/

	if (!(cut->istar = arr_alloc(num_istar, int)))
		err_msg("Allocation", "new_fea_cut", "istar");

	if (!(cut->beta = arr_alloc(num_x+1, double)))
		err_msg("Allocation", "new_fea_cut", "beta");

	cut->subfeaflag = FALSE;

	cut->alpha = 0.0;
	for (cnt = 0; cnt <= num_x; cnt++)
		cut->beta[cnt] = 0.0;

	return cut;
}

/***********************************************************************\
** This function frees the two arrays allocated for a single cut.
 ** When the incumbent cut is re-evaluated, its original value
 ** is erased and freed, and replaced by a new one.  Also, before
 ** exiting, the program must free all the cuts.
 \***********************************************************************/
void free_cut(one_cut *cut)
{

#ifdef LOOP
	printf("Inside free_cut\n");
#endif

	if (cut)
	{
		/*
		 printf("zl_free_cut ~1\n");
		 */
		if (cut->istar)
			mem_free(cut->istar);
		if (cut->beta)
			mem_free(cut->beta);
		if (cut->is_incumbent)
		{
			mem_free(cut->subobj_omega);
			mem_free(cut->subobj_freq);
		} /*added by Yifan 02/02/12 */
		mem_free(cut);
	}
}

/***********************************************************************\
** This function allocates memory for a new cut structure.  This entails
 ** the structure itself, and the _val_ array of one_cut pointers inside
 ** the structure.  The actual one_cut structures are allocated 
 ** according to the _num_betas_ parameter, via calls to new_cut().
 ** Note that the allocated cuts are NOT initialized.
 \***********************************************************************/
cut_type *new_cuts(int num_cuts, int num_x, int num_betas)
{
	cut_type *cuts;
	int cnt;

#ifdef TRACE
	printf("Inside new_cuts\n");
#endif

	if (!(cuts = (cut_type *) mem_malloc (sizeof(cut_type))))
		err_msg("Allocation", "new_cuts", "cuts");

	if (!(cuts->val = arr_alloc(num_cuts, cut_ptr)))
		err_msg("Allocation", "new_cuts", "cuts->val");

	cuts->cnt = num_betas;
	for (cnt = 0; cnt < num_betas && cnt < num_cuts; cnt++)
		cuts->val[cnt] = new_cut(num_x, 0, 0);

	return cuts;
}

/***********************************************************************\
**
 **
 \***********************************************************************/
void free_cuts(cut_type *cuts)
{
	int cnt;

#ifdef TRACE
	printf("Inside free_cuts\n");
#endif

	for (cnt = 0; cnt < cuts->cnt; cnt++)
		free_cut(cuts->val[cnt]);
	mem_free(cuts->val);
	mem_free(cuts);
}

/***********************************************************************\
 ** This function allocates memory for a new batch cuts structure.  This entails
 ** the structure itself, and the _batch_ array of _cuts_ pointers inside
 ** the structure.  The actual _cuts_ structures are allocated
 ** according to the bsize parameter, via calls to new_cuts().
 ** Note that not _one cut_ structure are allocated! Later we will point the
 ** _cuts_ structure here to those _one cut_  in compromise problem.
 \***********************************************************************/
batch_cut_type *new_bcuts(prob_type *p, int bsize, batch_cut_type *batch_cuts)
{
	/* modified by Yifan 2013.05.05 */
#ifdef TRACE
	printf("Inside new_bcuts\n");
#endif

	if (!(batch_cuts = (batch_cut_type *) mem_malloc (sizeof(batch_cut_type))))
		err_msg("Allocation", "batch_cut_type", "bcuts");

	if (!(batch_cuts->batch = arr_alloc(bsize, cuts_ptr)))
		err_msg("Allocation", "new_cuts", "cuts->val");

	batch_cuts->b_size = bsize;

	return batch_cuts;
}
/* modified by Yifan 2013.05.05 */
/* This function frees up all meomory allocated for batch cuts structure */
void free_bcuts(batch_cut_type *batch_cuts)
{
	int idx;

#ifdef TRACE
	printf("Inside free_bcuts\n");
#endif
	/* Free all cuts structure saved for the compromise problem Yifan 2013/01/17 */
	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
      if (batch_cuts->batch[idx] != NULL) {
        free_cuts(batch_cuts->batch[idx]);
      }
    }

	/*
	 for (idx = 0; idx < BATCH_SIZE; idx++) {
	 mem_free(sd_global->bcuts->batch[idx]);
	 }
	 */
	mem_free(batch_cuts->batch);
	mem_free(batch_cuts);

}

/***********************************************************************\
** This function prints the relevant information in a cut.
 ** It is meant to be used for debugging.
 \***********************************************************************/
void print_cut(cut_type *cuts, num_type *num, int idx)
{
	int cnt;

	printf("\nCut #%d:: c:%d o:%d\n  a:%f B:", idx, cuts->val[idx]->cut_obs,
			cuts->val[idx]->omega_cnt, cuts->val[idx]->alpha);
	for (cnt = 0; cnt <= num->mast_cols; cnt++)
		printf("%f ", cuts->val[idx]->beta[cnt]);
	printf("\nistar: ");
	for (cnt = 0; cnt < cuts->val[idx]->omega_cnt; cnt++)
		printf("%d ", cuts->val[idx]->istar[cnt]);
	printf("\n");
}

/**********************************************************************\
** This function prints some information of a cut for the purpose of
 ** cut index checking. zl
 \**********************************************************************/
void print_cut_info(cell_type *c, num_type *num, char *phrase)
{
	int cnt;
	int idx;
	int j;
	int num_rows; /* # of rows in CPLEX constraints. */
	int num_cuts; /* # of cuts in our data structure. */
	int m;
	int n;

	int status; /* indicating the CPLEX status. */
	int nzcnt;
	int rmatbeg[20];
	int rmatind[200];
	double rmatval[200];
	int surplus;

	FILE *g_FilePointer;

	g_FilePointer = fopen("cuts_idx.out", "a");
	fprintf(g_FilePointer, "-------------%s-------------- \n\n", phrase);

	num_rows = get_numrows(c->master); /* 2011.10.30 */

	if (0 == num_rows)
	{
		printf("Fail to get # of rows in CPLEX");
		exit(1);
	}
	fprintf(g_FilePointer, "\nCPLEX: # of rows = %d.", num_rows);

	num_cuts = num_rows - c->master->mar;
	fprintf(g_FilePointer, "\nMaster: # of rows = %d.", c->master->mar);
	fprintf(g_FilePointer, "\nOur structure: # of cuts = %d.\n", num_cuts);

	/*
	 int rmatspace; 
	 status = get_rows(c->master, &nzcnt, NULL, NULL, NULL,
	 0, &surplus, c->master->mar, num_rows-1);  
	 rmatspace = -surplus; 
	 if (status) {
	 fprintf(g_FilePointer, "\nFail to get rmatspace. Error #%d", status);
	 exit(1);
	 }  
	 
	 */

	status = get_rows(c->master, &nzcnt, rmatbeg, rmatind, rmatval, 200,
			&surplus, c->master->mar, num_rows - 1); /* 2011.10.30 */
	/*
	 fprintf(g_FilePointer, "\nsurplus = %d", surplus);
	 */
	if (status)
	{
		fprintf(g_FilePointer, "\nFail to get rows in CPLEX, error #%d.\n",
				status);
		exit(1);
	}

	for (idx = 0; idx < c->cuts->cnt; idx++)
	{

		/* Print cuts in our structure. */

		m = c->cuts->val[idx]->row_num; /*  the cut's row # in master constraint matrix. zl */

		fprintf(g_FilePointer,
				"\n***OurCut[#%d], Row # in master constraint matrix: %d.***\n",
				idx, m);
		for (cnt = 0; cnt <= num->mast_cols; cnt++)
			fprintf(g_FilePointer, "::%d: %f ", cnt,
					c->cuts->val[idx]->beta[cnt]);

		/* Print cuts in CPLEX structure. */

		fprintf(g_FilePointer, "\n+++CPLEX[#%d]+++\n", m);

		n = m - c->master->mar; /* the row # among the rows got from get_rows(). zl */

		if (n == c->cuts->cnt - 1)
		{
			for (j = rmatbeg[n]; j < nzcnt; j++)
				fprintf(g_FilePointer, "::%d: %f ", rmatind[j], rmatval[j]);
		}
		else
		{
			for (j = rmatbeg[n]; j < rmatbeg[n + 1]; j++)
			{
				fprintf(g_FilePointer, "::%d: %f ", rmatind[j], rmatval[j]);
			}
		}
		fprintf(g_FilePointer, "\n\n");
	}
	fclose(g_FilePointer);
}

