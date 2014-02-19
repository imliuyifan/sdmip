/***********************************************************************\
**
 ** optimal.c
 **
 ** This file contains the functions used in determining the optimality
 ** of a given solution.
 **
 **
 ** History:
 **   10 Mar 1992 - <Jason Mai> - created.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include <time.h>
#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "theta.h"
#include "optimal.h"
#include "omega.h"
#include "log.h"
#include "master.h"
#include "cuts.h"
#include "rvgen.h"
#include "sdconstants.h"
#include "sdglobal.h"

/***********************************************************************  
 ** This function determines whether or not the current incumbent 
 ** solution is considered to be optimal.  In the early iterations
 ** (less than MIN_ITER iterations), the function automatically assumes
 ** that the incumbent is NOT optimal.  Provided that MIN_ITER iterations
 ** have passed, it performs an optimality pre-test, to determine whether
 ** a full optimality check is a worthwhile pursuit.  If so, it does
 ** the full test, which involves reforming cuts and solving new master
 ** problems.  The function returns TRUE if the incumbent is optimal;
 ** FALSE otherwise.
 \***********************************************************************/
BOOL optimal(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s,
		double *phi)
{
#ifdef TRACE
	printf("Inside optimal\n");
#endif

	/*
	 printf("In optimal.c, full_test(), MIN = %d\n", sd_global->config.MIN_ITER);
	 */
	if (c->k > sd_global->config.MIN_ITER && *s->dual_statble_flag)
	{
		if (pre_test_1(sd_global, s))
		{ /* if (pre_test_2(p,c,s) )  JH */
			{
				if (sd_global->config.TEST_TYPE == 0)
					return TRUE;
				else if (full_test(sd_global, p, c, s))
				{
					printf(">");
					fflush(stdout);
					return TRUE;
				}
				else
				{
					printf("<");
					fflush(stdout);
				}
			}

		}

	}
#ifdef TRACE
	printf("Exiting optimal\n");
#endif
	return FALSE;
}

/***********************************************************************\  
 ** Because checking optimality is an arduous task, we first do a pre-
 ** check to determine if the full test is worthwhile.  This function
 ** determines whether the height at the candidate is close enough to 
 ** the height at the incumbent to warrant an optimality test.
 \***********************************************************************/
BOOL pre_test_1(sdglobal_type* sd_global, soln_type *s)
{
#ifdef TRACE
	printf("Inside pre_test_1\n");
#endif

	/* modified by zl, 08/17/04. */
	/* The candidate must be within some small percentage of incumbent cut */
	/*rare situation for s->candid_est < 0 and s->incumb_est > 0*/

	if (s->candid_est > 0)
	{
		s->smpl_test_flag = (s->candid_est
				> (1 - sd_global->config.PRE_EPSILON) * s->incumb_est);
	}
	else
	{
		s->smpl_test_flag = (s->candid_est
				> (1 + sd_global->config.PRE_EPSILON) * s->incumb_est);
	}

	/* Record if the simple test (pre_test) is ever passed. zl, 08/20/04. */
	if (s->smpl_test_flag)
		s->smpl_ever_flag = s->smpl_test_flag;

	return s->smpl_test_flag;

	/*
	 return 
	 (s->candid_est - s->incumb_est > 
	 -sd_global->config.CONFID_LO * s->incumb_stdev - sd_global->config.EPSILON);
	 */
}

/***********************************************************************\  
 ** If we pass the first pre_test, we undertake a second pre_test to 
 ** check to determine if the full test is worthwhile.  This function
 ** does bootstrapped based calculations to check whether or not the 
 ** full test is likely to be passed.  If this second pretest fails, the
 ** full test would fail as well, so it need not be performed. 
 **
 ** JH NOTE: Save bootstrap seed and re-initialize it if pre_test_2 passes
 \***********************************************************************/
BOOL pre_test_2(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	cut_type *T;
	double ULm;
	int *cdf, *observ;
	int m, sum;
	sd_long seed;

	int j;
	double ht, Sm;
	double eps_factor = 1.0; /* JH maybe 4? */
	double p_pass_factor = 1.0; /* JH maybe 0.8? */

#ifdef TRACE
	printf("Inside pre_test_2\n");
#endif

#ifdef RUN
	printf("\nPerforming the second pre-test \n\n");
#endif

#ifdef OPT
	printf("omega->cnt=%d. sigma->cnt=%d. lambda->cnt=%d.\n",
			s->omega->cnt, c->sigma->cnt, c->lambda->cnt);
	print_num(p->num);
#endif

	seed = sd_global->config.EVAL_SEED1;
	sum = 0;
	T = choose_cuts(p, c, s);
	observ = arr_alloc(c->k, int);
	cdf = arr_alloc(s->omega->most, int);

	empirical_distrib(s->omega, cdf);
	/*
	 printf("in pre_test2, seed = %d\n", sd_global->config.EVAL_SEED1);
	 fprintf(g_FilePointer, "\n*** in pre-test2: seed =  %d*** \n", sd_global->config.EVAL_SEED1);
	 */
	/* Find out how many of the resamplings satisfy optimality requirement */
	for (m = 0; m < sd_global->config.M; m++)
	{

		/* 
		 ** ??????????????????
		 ** What do we use here? Should we resample k omegas, or omega->cnt omegas?
		 ** As long as omega->weight is incremented when an omega is dropped,
		 ** we can use c->k.  (Jan 93 - They are, so we can).
		 ** ??????????????????
		 */
		sample_omega(sd_global, cdf, observ, c->k);
#ifdef OPT
		printf("SS-PRINT: passed sample_omega \n");
#endif
		reform_cuts(sd_global, c->sigma, s->delta, s->omega, p->num, T, observ,
				c->k);
#ifdef OPT
		printf("SS-PRINT: passed reform_cuts \n");
#endif

		/* Find the highest reformed cut at the incumbent solution */
		Sm = T->val[0]->alpha
				- CxX(T->val[0]->beta, s->incumb_x, p->num->mast_cols);
		for (j = 1; j < T->cnt; j++)
		{
			ht = T->val[j]->alpha
					- CxX(T->val[j]->beta, s->incumb_x, p->num->mast_cols);
			if (Sm < ht)
				Sm = ht;
		}
		Sm += CxX(p->c, s->incumb_x, p->num->mast_cols);

		ULm = calculate_ULm(p, c, T, s);

#ifdef OPT
		printf("ULm=%lf, Sm=%lf\n", ULm, Sm);
#endif
		/*
		 printf("ULm=%lf, Sm=%lf inc = %lf Sm-ULm/inc = %lf \n", 
		 ULm, Sm, s->incumb_est, (Sm-ULm)/s->incumb_est);
		 */

		/* JH: What value for epsilon here?  */
		if ((Sm - ULm) / s->incumb_est
				<= eps_factor * sd_global->config.PRE_EPSILON)
			sum++;

		/*  Skip out of the loop if there's no hope of meeting condition. - SS  
		 *** JH deleting for now 6/10/02 ***
		 if (sum + sd_global->config.M - m - 1 < sd_global->config.PERCENT_PASS*sd_global->config.M) return FALSE;
		 if (sum + sd_global->config.M - m - 1 < 0.75*sd_global->config.M)
		 {      sd_global->config.MIN_ITER = c->k + 20;
		 return (FALSE); 
		 } 
		 */

	}

#ifdef OPT
	printf("\nsum=%d\n",sum);
#endif
#ifdef RUN
	printf("\nc->k = %d, sum_pre2 = %d\n", c->k, sum);
#endif

	mem_free(cdf);
	mem_free(observ);
	free_cuts(T);

	/* JH: What value for percent_pass here?  */
	/*  return (sum >= p_pass_factor*sd_global->config.PERCENT_PASS * sd_global->config.M); */
	if (sum
			>= p_pass_factor * sd_global->config.PERCENT_PASS
					* sd_global->config.M)
	{
		sd_global->config.EVAL_SEED1 = seed;
		/*
		 scanf("n = %d");
		 */
		return (TRUE);
	}
	/* this else shouldn't be executed, given SS's pass above */
	/*  jh whoopsie 6/10/02 
	 else { if (sum <= 0.75*p_pass_factor*sd_global->config.PERCENT_PASS*sd_global->config.M)
	 sd_global->config.MIN_ITER = c->k + 20;
	 return (FALSE); 
	 } 
	 */
	else
		return (FALSE); /* part of whoopsie */

}

/***********************************************************************\  
 ** This function performs a complete statistical test of optimality.
 ** First, it selects cuts whose height at the incumbent is "close" to
 ** incumbent cut's height.  Then, it performs M resamplings of the
 ** observations in omega, and reforms the selected cuts with respect
 ** to these observations (as if *they* were observed instead of the actual
 ** omega).  For each of the M resamplings, a master program containing the
 ** reformed cuts is solved, and if almost all of the solutions to these
 ** master programs 
 \***********************************************************************/
BOOL full_test(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	cut_type *T;
	double Lm;
	int *cdf, *observ;
	int m, sum;
#if 0
	int i; /* for printing cuts information. zl */
#endif
	int j;
	double ht, Sm;
#ifdef CAL_CHECK
	double temp;
#endif
	int num_failed = 0; /* Number of failed replications. */
	double error_sum = 0.0;/* The sum of errors of failed replications. 
	 added by zl, 0817/04. */
	clock_t start, end; /* Recording time on performing full
	 optimality test. zl, 06/29/04. */

#ifdef TRACE
	printf("Inside full_test\n");
#endif

#ifdef RUN
	/* printf("\nPerforming the full optimality test at iter %d\n\n", c->k);*/
#endif

#ifdef OPT
	printf("omega->cnt=%d. sigma->cnt=%d. lambda->cnt=%d.\n",
			s->omega->cnt, c->sigma->cnt, c->lambda->cnt);
	print_num(p->num);
#endif

#ifdef CAL_CHECK
	/* calculate the lower bound based on the original cuts, and compare it
	 with the incumbant estimate, for checking purpose only. zl */

	temp = cal_temp_lb(p, c, s, c->cuts);
	printf("temp is : %f\n",temp);
	printf("incumbent est is: %f\n",s->incumb_est);
	fprintf(g_FilePointer, "temp = %f\n", temp);

	fprintf(g_FilePointer, "*** Enter full_test(), num_cuts = %d ***\n",
			c->cuts->cnt);
#endif

	Sm = c->cuts->val[0]->alpha
			- CxX(c->cuts->val[0]->beta, s->incumb_x, p->num->mast_cols);

	for (j = 1; j < c->cuts->cnt; j++)
	{
		ht = c->cuts->val[j]->alpha
				- CxX(c->cuts->val[j]->beta, s->incumb_x, p->num->mast_cols);
		if (Sm < ht)
			Sm = ht;
	}

#ifdef CAL_CHECK
	printf("temp LM is : %f\n",temp);
	printf("SM regular is around %f\n",Sm);
	printf("incumbent est is: %f\n",s->incumb_est);
#endif

	Sm = cut_height(sd_global, c->cuts->val[0], s->incumb_x, c, p->num);

	for (j = 1; j < c->cuts->cnt; j++)
	{
		ht = cut_height(sd_global, c->cuts->val[j], s->incumb_x, c, p->num);

		if (Sm < ht)
			Sm = ht;
	}

#ifdef CAL_CHECK
	printf("temp LM is : %f\n",temp);
	printf("SM cut_height is around %f\n",Sm);
	printf("incumbent est is: %f\n",s->incumb_est);
#endif

	/* Record the time spent on performing full optimality test. 
	 zl, 06/29/04. */
	start = clock();

	sum = 0;
	T = choose_cuts(p, c, s);
	observ = arr_alloc(c->k, int);
	cdf = arr_alloc(s->omega->most, int);

	empirical_distrib(s->omega, cdf);
	/*
	 printf(" In full_test seed = %d\n", sd_global->config.EVAL_SEED1);
	 */
	/* Find out how many of the resamplings satisfy optimality requirement */

	for (m = 0; m < sd_global->config.M; m++)
	{

		/* 
		 ** ??????????????????
		 ** What do we use here? Should we resample k omegas, or omega->cnt omegas?
		 ** As long as omega->weight is incremented when an omega is dropped,
		 ** we can use c->k.  (Jan 93 - They are, so we can).
		 ** ??????????????????
		 */

		sample_omega(sd_global, cdf, observ, c->k);
#ifdef OPT
		printf("SS-PRINT: passed sample_omega \n");
#endif
		reform_cuts(sd_global, c->sigma, s->delta, s->omega, p->num, T, observ,
				c->k);
#ifdef OPT
		printf("SS-PRINT: passed reform_cuts \n");
#endif

#if 0
		/* Print out the reformed cuts for checking purpose only.  zl */
		fprintf(g_FilePointer, "\n*** after reform_cuts ***\n");
		fprintf(g_FilePointer, "m = %d, num_cuts (T) = %d\n", m, T->cnt);
		for (i=0; i<T->cnt; i++)
		{
			fprintf(g_FilePointer, "++ T ++:  i = %d, alpha = %f, row_num = %d, theta = %f\n", i,
					T->val[i]->alpha,
					T->val[i]->row_num,
					s->Master_pi[T->val[i]->row_num]);
			for (j=0; j<=p->num->mast_cols; j++)
			fprintf(g_FilePointer, "beta[%d] = %f, ", j,
					T->val[i]->beta[j]);
			fprintf(g_FilePointer, "\n");
		}
#endif
		/* Find the highest reformed cut at the incumbent solution */

		Sm = T->val[0]->alpha
				- CxX(T->val[0]->beta, s->incumb_x, p->num->mast_cols);

		/* Yifan 04/01/2012 No t/k for the cut height of previous cuts */
		/* modified by Yifan 2013.05.06 */
		/*Sm = Sm * (double)T->val[0]->cut_obs / (double)c->k;*/

		/* Yifan 03/17/2012 Updated for optimality cut height*/
		/*
		 if (sd_global->config.LB_TYPE == 1) {
		 Sm += (1- (double)T->val[0]->cut_obs / (double)c->k) * sd_global->Eta0;
		 }*/

		for (j = 1; j < T->cnt; j++)
		{
			ht = T->val[j]->alpha
					- CxX(T->val[j]->beta, s->incumb_x, p->num->mast_cols);

			/* modified by Yifan 2013.05.06 */
			/* ht = ht * (double)T->val[j]->cut_obs / (double)c->k;*/

			/*
			 if (sd_global->config.LB_TYPE == 1) {
			 ht += (1- (double)T->val[j]->cut_obs / (double)c->k) * sd_global->Eta0;
			 }*/

			if (Sm < ht)
				Sm = ht;
		}

		/* In QP approach, we don't include the incumb_x * c in Sm. zl */
		if (sd_global->config.MASTER_TYPE == SDLP)
			Sm += CxX(p->c, s->incumb_x, p->num->mast_cols);

		if (sd_global->config.MASTER_TYPE == SDLP)
			Lm = solve_temp_master(sd_global, p, T, c);
		else
			Lm = cal_temp_lb(p, c, s, T);

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "m = %d, Sm = %f, Lm = %f\n", m+1, Sm, Lm);
		printf("m = %d, Sm = %f, Lm = %f", m+1, Sm, Lm);
#endif

#ifdef OPT
		printf("Lm=%lf, Sm=%lf\n", Lm, Sm);
#endif

		if (DBL_ABS((Sm - Lm) / s->incumb_est)
						<= sd_global->config.EPSILON)
		{
			sum++;
		}
		else
		{
			/* Sum up errors of failed replications in full test.
			 Only at Max_iter. zl, 08/17/04. */
			if (c->k >= p->num->iter)
			{
				num_failed++;
				error_sum += (Sm - Lm) / s->incumb_est;
			}
		}

#ifdef OPT_TEST
        printf("m = %d, UB = %f, LB = %f, sum = %d\n", m+1, Sm, Lm, sum);
		printf("s->incumb_est is: %f\n",s->incumb_est);
		printf("(Sm-Lm)/incumb_est is : %f\n",(Sm - Lm) / s->incumb_est);
		printf("**********************************************\n");
#endif

		if (m + 1 - sum
				>= (1 - sd_global->config.PERCENT_PASS) * sd_global->config.M)
		{
			/* Record average error of failed replications in full test.
			 Only at Max_iter. zl, 08/17/04. */
			if (c->k >= p->num->iter)
				s->full_test_error = error_sum / num_failed;

			/* Record the time spent on performing full optimality test. 
			 zl, 10/20/04. */
			end = clock();
			s->run_time->full_test_iter = ((double) (end - start))
					/ CLOCKS_PER_SEC;
			s->run_time->full_test_accum += s->run_time->full_test_iter;
#if 0
			printf("iter %d: full_test_iter = %lf, full_test_accum = %lf\n",
					c->k,
					s->run_time->full_test_iter,
					s->run_time->full_test_accum);
#endif

			/* Yifan 03/19/2012 Avoid memory leaks*/
			mem_free(observ);
			mem_free(cdf);
			free_cuts(T);
			return FALSE;
		}

		/*  Skip out of the loop if there's no hope of meeting condition. - SS  
		 if (sum + sd_global->config.M - m - 1 < sd_global->config.PERCENT_PASS*sd_global->config.M) return FALSE;
		 */

		/*
		 if ((s->incumb_est - Lm) / Lm <= sd_global->config.EPSILON)
		 sum++;
		 */
	}

	/* Record the number of bootstraped replications that satisfy the optimality 
	 requirement. zl, 06/30/04. */
	s->passed = sum;

	/* Write down all the relevant information on the modification on the
	 stopping rule. zl */

#ifdef OPT
	printf("\nsum=%d\n",sum);
#endif
#ifdef RUN
	/* printf("\nc->k = %d, sum_full = %d\n", c->k, sum);*/
#endif

	/* Record the time spent on performing full optimality test. 
	 zl, 06/29/04. */
	end = clock();
	s->run_time->full_test_iter = ((double) (end - start)) / CLOCKS_PER_SEC;
	s->run_time->full_test_accum += s->run_time->full_test_iter;
#if 0
	printf("iter %d: full_test_iter = %lf, full_test_accum = %lf\n",
			c->k,
			s->run_time->full_test_iter,
			s->run_time->full_test_accum);
#endif

	/* Yifan 03/19/2012 Avoid memory leaks*/
	mem_free(observ);
	mem_free(cdf);
	free_cuts(T);

	return (sum >= sd_global->config.PERCENT_PASS * sd_global->config.M);
}

/***********************************************************************\  
 ** This function selects all cuts whose height at the incumbent x
 ** is close to the height of the incumbent cut.  These cuts together
 ** are likely to provide good approximations of f at incumb_x, when
 ** they are reformed with new observations.  The function returns a
 ** new cut structure which contains room for cuts to be reformed.
 ** Only the _istar_ and _cut_obs_ fields of each cut have been initialized.
 \***********************************************************************/
cut_type *choose_cuts(prob_type *p, cell_type *c, soln_type *s)
{
	cut_type *T;
	int cnt;

	T = new_cuts(p->num->max_cuts, p->num->mast_cols, 0);

	/* Define what it means to be "close" to the incumbent */
//	sd_global->config.CONFID_HI * s->incumb_stdev;
	CxX(p->c, s->incumb_x, p->num->mast_cols);

	/* removed by Yifan 05/26/2013 */
	/*
	 printf("\nc->cuts->cnt is :%d\n", c->cuts->cnt);
	 printf("c->feasible_cuts_added->cnt is :%d\n",c->feasible_cuts_added->cnt);
	 printf("c->feasible_cuts_pool->cnt is :%d\n\n",c->feasible_cuts_pool->cnt);*/

#ifdef OPT
	printf("Cutoff=%lf, incumb_est=%lf.\n", cutoff, s->incumb_est);
#endif

	/* Loop through the cuts available for reforming; pick close ones */
	for (cnt = 0; cnt < c->cuts->cnt; cnt++)
	{ /*added by Yifan to avoid select feasibility cut for opt test 02/26/2012*/

		/* Yifan 04/03/2012 Choosing cuts with nonzero dual multipliers */
		/*if(c->cuts->val[cnt]->subfeaflag == FALSE) continue;
		 height = cut_height(c->cuts->val[cnt], s->incumb_x, c, p->num);*/

#ifdef OPT
		printf("%d height=%lf\n", cnt, height);
#endif

		/*if( DBL_ABS(incumb_ht - height) <= 0.25*DBL_ABS(incumb_ht) )*/

		/* Yifan 04/03/2012 Choosing cuts with nonzero dual multipliers */
		if (s->Master_pi[c->cuts->val[cnt]->row_num + 1] > 0.00001)
		{
			T->val[T->cnt] = new_cut(p->num->mast_cols,
					c->cuts->val[cnt]->omega_cnt, c->cuts->val[cnt]->cut_obs);
			copy_arr(T->val[T->cnt]->istar, c->cuts->val[cnt]->istar,
					T->val[T->cnt]->omega_cnt - 1)
			T->val[T->cnt]->row_num = c->cuts->val[cnt]->row_num; /* JH 5/98 */

#ifdef OPT
			printf("Taking the cut:");
			print_cut(c->cuts, p->num, cnt);
#endif

			T->cnt++;
		}
	}

#ifdef OPT
	printf("Selected %d cuts\n", T->cnt);
#endif

	return T;
}

/***********************************************************************\  
 ** This function forms an empirical distribution on the observations
 ** stored in omega, and calculates an integer cdf to represent the
 ** distribution.  An observation which has been seen n times will have n 
 ** times the probability of being chosen as an observation seen only once.  
 \***********************************************************************/
void empirical_distrib(omega_type *omega, int *cdf)
{
	int cnt;

	/* Calculate an integer cdf distribution for observations */
	/* If the cnt is not a valid omega idx, we know that weight[cnt] is 0 */
	cdf[0] = omega->weight[0];
	for (cnt = 1; cnt < omega->most; cnt++)
		cdf[cnt] = cdf[cnt - 1] + omega->weight[cnt];
}

/***********************************************************************\  
 ** This function randomly selects a new set of observations from the 
 ** old set of observations stored in omega.  Entries in omega which 
 ** have been observed multiple times have a proportionally higher 
 ** chance of being selected for the new set.  The function fills an
 ** array, assumed to be of a size equal to the number of iterations 
 ** which have passed, with the new set of observations.
 \***********************************************************************/
void sample_omega(sdglobal_type* sd_global, int *cdf, int *observ, int k)
{
	int cnt, obs;
	int sample;

#ifdef LOOP
	printf("Inside sample_omega\n");
#endif

#ifdef OPT
	printf("Sample=");
#endif

	/*
	 ** ???????????????????
	 ** The k here better be less than the number of omegas!
	 ** Since omega->weight is incremented when another is dropped, the
	 ** sum stays the same...
	 ** ???????????????????
	 */

	/* Choose k observations according to cdf (k = number of iterations) */
	for (obs = 0; obs < k; obs++)
	{
		sample = randfun(k, &(sd_global->config.EVAL_SEED1));
		/*
		 fprintf(g_FilePointer, "obs = %d, sample = %d\n", obs, sample);    
		 */
		for (cnt = 0; sample > cdf[cnt]; cnt++)
			/* Loop until sample falls below cdf */;
		observ[obs] = cnt;

#ifdef OPT
		printf("%d ", cnt);
#endif
	}

#ifdef OPT
	printf(" ***\n");
#endif

}

/*********************************************************\
** This function calculates an upper bound on a lower 
 ** bound on the value of the resampled master problem, 
 ** using the dual solution ** to the actual master problem 
 ** and an exact penalty function.    
 ** JH 5/98  WORK
 \*********************************************************/
double calculate_ULm(prob_type *p, cell_type *c, cut_type *T, soln_type *s)
{
	int j, t;
	double zbar, dbar, value;

	/* modified.  zl 
	 lpwrite(c->master->lp,"master.lp");
	 lpwrite(c->subprob->lp,"sub.lp"); 
	 */
	write_prob(c->master, "master.lp"); /* 2011.10.30 */
	write_prob(c->subprob, "sub.lp"); /* 2011.10.30 */

	zbar = s->candid_est;
	for (t = 0; t < c->cuts->cnt; ++t)
	{ /* printf(" row = %d ; alpha = %lf ; lambda = %lf \n", 
	 c->cuts->val[t]->row_num, c->cuts->val[t]->alpha, 
	 s->Master_pi[c->cuts->val[t]->row_num+1]); */
		zbar -= c->cuts->val[t]->alpha
				* s->Master_pi[c->cuts->val[t]->row_num + 1];

		/*printf(" zbar was %lf, and is now %lf \n", s->candid_est, zbar); */
	}

	value = zbar;
	for (t = 0; t < T->cnt; ++t)
	{ /* printf(" row = %d ; alpha_hat = %lf ; lambda = %lf \n", 
	 T->val[t]->row_num, T->val[t]->alpha, 
	 s->Master_pi[T->val[t]->row_num+1]); */
		value += T->val[t]->alpha * s->Master_pi[T->val[t]->row_num + 1];
	}
	/* printf(" Value is initialized at %lf \n", value); */

	for (j = 1; j <= p->num->mast_cols; ++j)
	{
		dbar = s->Master_dj[j];
		for (t = 0; t < c->cuts->cnt; ++t)
		{ /* printf(" lambda[%d] = %lf  beta[%d,%d] = %lf \n",
		 c->cuts->val[t]->row_num, s->Master_pi[c->cuts->val[t]->row_num+1],
		 c->cuts->val[t]->row_num, j, c->cuts->val[t]->beta[j]); */
			dbar += c->cuts->val[t]->beta[j]
					* s->Master_pi[c->cuts->val[t]->row_num + 1];
		}
		/*printf("\n dbar[%d] was %lf and is now %lf \n", j, s->Master_dj[j],dbar);
		 */

		for (t = 0; t < T->cnt; ++t)
		{ /* printf(" lambda[%d] = %lf  beta[%d,%d] = %lf \n",
		 T->val[t]->row_num, s->Master_pi[T->val[t]->row_num+1],
		 T->val[t]->row_num, j, T->val[t]->beta[j]);  */
			dbar -= T->val[t]->beta[j] * s->Master_pi[T->val[t]->row_num + 1];
		}
		/* printf("\n OK, dbar[%d] is now %lf and x[%d] is %lf \n", 
		 j, dbar, j, s->candid_x[j]); */

		/*if (dbar < 0 )*/value += dbar * s->candid_x[j];
		/* printf("\n value is now %lf \n\n", value); */

	}
	return (value);
}

/***********************************************************************\
** This function loads, solves, and frees a master program with cuts
 ** specified by _T_.  It returns the value of the objective function
 ** at the optimal solution. 
 \***********************************************************************/
double solve_temp_master(sdglobal_type* sd_global, prob_type *p, cut_type *T,
		cell_type *c)
{
	one_problem *temp;
	double ans;

	/*
	 malloc_verify();
	 printf("Before new_master\n");
	 */
	temp = orig_new_master(p->master, T, 0);

	/*
	 malloc_verify();
	 printf("Before solve_problem\n");
	 */

	solve_problem(sd_global, temp);
	c->LP_test++;

	/*
	 malloc_verify();
	 printf("Before get_objective\n");
	 */

	ans = get_objective(temp);

	/*
	 malloc_verify();
	 printf("Before remove_problem\n");
	 */

	remove_problem(temp);

	/*
	 malloc_verify();
	 printf("Before free_one_prob\n");
	 */

	free_one_prob(temp);

	/*
	 malloc_verify();
	 printf("Before returning\n");
	 */

	return ans;
}

/***********************************************************************\  
 ** This function will calculate a new set of cuts based on the 
 ** observations of omega passed in as _observ_, and the istar's 
 ** which have already been stored in the _istar_ field of each cut.  
 ** If an istar field does not exist for a given observation, then 
 ** a value of zero is averaged into the calculation of alpha & beta.
 \***********************************************************************/
void reform_cuts(sdglobal_type* sd_global, sigma_type *sigma, delta_type *delta,
		omega_type *omega, num_type *num, cut_type *T, int *observ, int k)
{
	int cnt, obs, idx, count; /* modified by Yifan 2013.05.06 */
	i_type istar;

#ifdef TRACE
	printf("Inside reform_cuts\n");
#endif
#ifdef OPT
	printf("SS-PRINT: Inside reform_cuts\n");
	printf("SS-PRINT: T->cnt=%d\n",T->cnt);
#endif

	/* Loop through all the cuts and reform them */
	for (cnt = 0; cnt < T->cnt; cnt++)
	{

#ifdef OPT
		printf("\nReforming #%d\n", cnt);
		malloc_verify();
#endif

		/* Begin with cut coefficients of zero */
		for (idx = 0; idx <= num->mast_cols; idx++)
			T->val[cnt]->beta[idx] = 0.0;
		T->val[cnt]->alpha = 0.0;

		count = 0; /* modified by Yifan 2013.05.06 */
		/* Reform this cut based on resampled observations */
		for (obs = 0; obs < k; obs++)
		{
			/* Only compute if this observation actually exists in omega */
			if (valid_omega_idx(omega, observ[obs]))
			{
				/* Only sum values if the cut has an istar for this observation */
				if (observ[obs] < T->val[cnt]->omega_cnt&&
				T->val[cnt]->istar[observ[obs]] != DROPPED){
				istar.sigma = T->val[cnt]->istar[observ[obs]];
				istar.delta = sigma->lamb[istar.sigma];

				T->val[cnt]->alpha += sigma->val[istar.sigma].R +
				delta->val[istar.delta][observ[obs]].R;

				for (idx = 1; idx <= num->nz_cols; idx++)
				T->val[cnt]->beta[sigma->col[idx]] +=
				sigma->val[istar.sigma].T[idx];

				for (idx = 1; idx <= num->rv_cols; idx++)
				T->val[cnt]->beta[delta->col[idx]] +=
				delta->val[istar.delta][observ[obs]].T[idx];
				count++; /* modified by Yifan 2013.05.06 */
			}
		}
	}

				/* Take the average of the alpha and beta values */
				/* modified by Yifan 2013.05.06 */
		for (idx = 0; idx <= num->mast_cols; idx++)
			T->val[cnt]->beta[idx] /= (double) k;
		T->val[cnt]->alpha /= (double) k;
		if (sd_global->config.LB_TYPE == 1)
		{
			T->val[cnt]->alpha += (1 - (double) count / (double) k)
					* sd_global->Eta0;
		}

#ifdef OPT
		print_cut(T, num, cnt);
#endif

	}
}

/*
 ** This function returns a uniform random number between [0, greatest-1]
 ** using our own random number generator.  This is not so good, since
 ** we are all running off the same seed... we ought to have different 
 ** streams of random numbers.
 */
int randfun(int greatest, sd_long *seed)
{
	return (int) (randUniform(seed) * greatest);
}

/****************************************************************************\
  This function is to calculate the lower bound on the optimal value which
 is used in stopping rule in full_test() in optimal.c in the case of
 regularized approach. 
 \****************************************************************************/
double cal_temp_lb(prob_type *p, cell_type *c, soln_type *s, cut_type *T)
{
	double *bk; /* vector: b - A*incumb_x. */
	double *lambda; /* vector: the dual of the primal constraints. */
	double bk_lambda; /* scalar: bk*lambda. */

	sparse_matrix *A_Trans; /* sparse_matrix: the transpose of A. */
	double *A_Trans_lambda; /* vector: - A_Trans * lambda. */

	double *theta; /* the dual of the reformed cut constraints. */
	double *Vk; /* the vector of scalars of cut constraints. 
	 Vk = alpha - beta * incumb_x  */
	double Vk_theta; /* scalar: Vk*theta. */

	double *Bk_theta; /* vector: Bk_Transpose * theta, where Bk_Transpose is
	 the matrix of cut coefficients. */
	double *Bk_col; /* vector: Store one column of Bk while calculating
	 Bk_theta. */
	double *q_vec; /* vector: c + Bk_theta - A_Trans_lambda. */
	double q_term; /* scalar: q_vec * q_vec. */
	double ht;
	double *lb; /* vector: store the lower bounds.  */
	double *ub; /* vector: store the upper bounds.  */

	double Sm; /* The upper bound of the optimal value. */
	double Lm; /* The calculated lower bound of the optimal value. */
	/* double eta_zero; */
	int cnt;
#ifdef CAL_CHECK
	int idx;
#endif
	int i;

#ifdef TRACE
	printf ("Inside cal_temp_lb.\n");
#endif

	if (!(bk =
			arr_alloc(p->num->mast_rows + c->feasible_cuts_added->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "bk");
	if (!(lambda =
			arr_alloc(p->num->mast_rows + c->feasible_cuts_added->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "lambda");
	/* Yifan 03/19/2012 Test update A and B matrix*/
	if (!(A_Trans_lambda = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "cal_temp_lb", "A_lambda");
	if (!(A_Trans = (sparse_matrix *) mem_malloc(sizeof(sparse_matrix))))
		err_msg("Allocation", "cal_temp_lb", "A_Trans");
	if (!(A_Trans->val = arr_alloc(p->A->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "A_Trans->val");
	if (!(A_Trans->row = arr_alloc(p->A->cnt+1, int)))
		err_msg("Allocation", "cal_temp_lb", "A_Trans->row");
	if (!(A_Trans->col = arr_alloc(p->A->cnt+1, int)))
		err_msg("Allocation", "cal_temp_lb", "A_Trans->col");
	if (!(theta = arr_alloc(T->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "theta");
	if (!(Vk = arr_alloc(T->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "Vk");
	/* Yifan 03/19/2012 Test update A and B matrix*/
	if (!(Bk_theta = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "cal_temp_lb", "Bk_theta");
	if (!(Bk_col = arr_alloc(T->cnt+1, double)))
		err_msg("Allocation", "cal_temp_lb", "Bk_col");
	/* Yifan 03/19/2012 Test update A and B matrix*/
	if (!(q_vec = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "cal_temp_lb", "q_vec");
	if (!(lb = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "cal_temp_lb", "lb");
	if (!(ub = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "cal_temp_lb", "ub");

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n****** Entering full_test()---cal_temp_lb ******\n");
#endif

	/* Calculate upper bound, Sm. */

	/*modified by Yifan to avoid evaluating feasibility cut 02/26/2011*/
	if (T->val[0]->subfeaflag == TRUE)
	{
		Sm = T->val[0]->alpha
				- CxX(T->val[0]->beta, s->incumb_x, p->num->mast_cols);
	}
	else
	{
		Sm = -INFBOUND;
	}

	/* The adjustment when using the original cuts instead of the reformed cuts
	 to check the correctness of the lb calculation. 
	 Sm *= (double)T->val[0]->cut_obs / (double)c->k; 
	 */
	for (cnt = 1; cnt < T->cnt; cnt++)
	{
		/*added by Yifan to avoid evaluating feasibility cut 02/26/2011*/
		if (T->val[cnt]->subfeaflag == FALSE)
			continue;

		ht = T->val[cnt]->alpha
				- CxX(T->val[cnt]->beta, s->incumb_x, p->num->mast_cols);

		/* The adjustment when using the original cuts instead of the reformed cuts
		 to check the correctness of the lb calculation.
		 ht *= (double)T->val[cnt]->cut_obs / (double)c->k; 
		 */
		if (Sm < ht)
			Sm = ht;
	}
	/*  
	 Sm += CxX(p->c, s->incumb_x, p->num->mast_cols);
	 */
	/*
	 ** Calculate bk, which is A*incumb_x - b. Note: in fact, we are
	 ** calculating -bk here, due to the way function TxX works. Also be aware
	 ** of the one-norm. 
	 */

	/* Later used for updating rhs of Lower Bound constraint */

	/* modified by Yifan 2013.05.06 */
	/*
	 if (sd_global->config.LB_TYPE == 1) {
	 eta_zero = sd_global->Eta0;
	 }
	 else {
	 eta_zero = 0.0;
	 }*/

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n+++ bk*lambda +++\n");
#endif
	/* 1. */
	/* Print out rhs. */
	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		bk[cnt + 1] = p->master->rhsx[cnt];

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "rhs[%d] = %f, ", cnt+1, bk[cnt+1]);
#endif
	}

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");

	/* Print out incumb_x. */
	for (cnt=0; cnt<p->num->mast_cols+1; cnt++)
	fprintf(g_FilePointer, "incumb_x[%d] = %f, ", cnt, s->incumb_x[cnt]);
	fprintf(g_FilePointer, "\n");
#endif

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	/* Calculate bk = b - A * incumb_x. */
	TxX(p->A, s->incumb_x, bk);

	/* Yifan 03/11/2012 Updated for Feasibility Cuts*/
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		bk[p->num->mast_rows + cnt + 1] =
				c->feasible_cuts_added->val[cnt]->alpha
						- CxX(c->feasible_cuts_added->val[cnt]->beta,
								s->incumb_x, p->num->mast_cols);
	}

#ifdef CAL_CHECK
	for (cnt=0; cnt<p->num->mast_rows; cnt++)
	fprintf(g_FilePointer, "bk[%d] = %f, ", cnt+1, bk[cnt+1]);
	fprintf(g_FilePointer, "\n");

	/* Print out s->Master_pi. */
	for (cnt=0; cnt<=p->num->mast_rows+c->cuts->cnt; cnt++)
	fprintf(g_FilePointer, "Master_pi[%d] = %f ,", cnt,
			s->Master_pi[cnt]);
	fprintf(g_FilePointer, "\n");
#endif

	/* 2. */

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	/* Obtain lambda from s->Master_pi of original master problem constraints. */
	/* Dual values' sign need to be flipped here before assigning to lambda*/
	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		lambda[cnt + 1] = -s->Master_pi[cnt + 1];

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "lambda[%d] = %f, ", cnt+1, lambda[cnt+1]);
#endif
	}

	/* Yifan 03/11/2012 Updated for Feasibility Cuts NOTICEHERE*/
	/* Obtain lambda from s->Master_pi of  master problem feasibility cuts. */
	/* Dual values' sign need to be flipped here before assigning to lambda*/
	/* Yifan 03/19/2012 Test removing eta0*/
	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		lambda[p->num->mast_rows + cnt + 1] =
				-s->Master_pi[c->feasible_cuts_added->val[cnt]->row_num + 1];
	}

#ifdef OPT_TEST
	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		printf("lambda_m[%d]: %f\t",cnt+1, lambda[cnt+1]);
		printf("bk_m[%d]: %f\n",cnt+1,bk[cnt+1]);
	}

	for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
	{
		printf("lambda_f[%d]: %f\t", p->num->mast_rows+cnt+1, lambda[ p->num->mast_rows+cnt+1]);
		printf("bk_f[%d]: %f\n", p->num->mast_rows+cnt+1,bk[ p->num->mast_rows+cnt+1]);
	}
#endif

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");

	for (i=1; i<=p->num->mast_cols; i++)
	fprintf(g_FilePointer, "dj[%d] = %f, ", i, s->Master_dj[i]);
	fprintf(g_FilePointer, "\n");
#endif

	/* 
	 ** Calculate bk_lambda = bk * lambda. 
	 */

	/* 3. */

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	/* Yifan 03/11/2012 Updated for Feasibility Cuts*/
	bk_lambda = CxX(bk, lambda,
			p->num->mast_rows + c->feasible_cuts_added->cnt);

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "bk_lambda = %f\n", bk_lambda);
#endif

	/* Be careful here!!! We should also include the bound constraints of the 
	 variables (this is in fact the negative of identity matrix) into the 
	 constraint matrix (A) and the lower bounds into the rhs vector (bk). 
	 - I * x <= -bk
	 The dual of these constraints are the reduced costs, dj. */

	/*modified by Yifan to get bk_lambda correctlly considering lower and upper bounds*/
	//get_lb(lb, c->master, p->num->mast_cols);
	//  get_ub(ub, c->master, p->num->mast_cols);
	/*
	 #ifdef CAL_CHECK
	 for (i=1; i<=p->num->mast_cols; i++) 
	 fprintf(g_FilePointer, "lb[%d] = %f, ", i, lb[i]);
	 fprintf(g_FilePointer, "\n");
	 #endif
	 */
	/*!!!!!! THIS PART IS PROBLEMATIC */
	/*THINK AGAIN, NO NEED FOR CORRECTION */
	/*  lb <= X_k <=ub */
	/* If it is touching the lb, then X_k - lb is 0 */
	/* If it is touching the ub, then X_k - ub is 0 */
	/* If it is touching no bound, then dj is 0 */
	/* THEREFORE, ( X_k - b ) * dj is 0 in any case */

	/* Yifan 03/14/2012 Make sure addition or subtraction*/
	/*
	 for (i=1; i<p->num->mast_cols+1; i++) 
	 bk_lambda += s->Master_dj[i] * s->incumb_x[i];*/

	//bk_lambda -= CxX(lb, s->Master_dj, p->num->mast_cols);
	//bk_lambda -= CxX(ub, s->Master_dj, p->num->mast_cols);
	/*modified by Yifan to get bk_lambda correctlly considering lower and upper bounds*/

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "bk_lambda = %f\n", bk_lambda);
#endif

	/* 
	 ** Calculate A_Trans*lambda. Note: in fact, we are calculating
	 ** - A_Trans*lambda here, due to the way function TxX works. 
	 */

	/* Print out the A matrix. */
	/*     
	 fprintf(g_FilePointer, "+++ A matrix +++\n");
	 for (cnt=1; cnt<=p->A->cnt; cnt++) 
	 fprintf(g_FilePointer, "cnt = %d, val = %f, row = %d, col = %d\n",
	 cnt, p->A->val[cnt], p->A->row[cnt], p->A->col[cnt]);
	 */

	/* Get A Transpose from A. */
	/*
	 printf("p->num->mast_rows: %d\n",p->num->mast_rows);
	 printf("cell->cuts->cnt is:%d\n",c->cuts->cnt);
	 printf("p->A->cnt is %d\n",p->A->cnt);
	 debug Yifan*/

	/* 4. */

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	A_Trans->cnt = p->A->cnt;

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\nBefore A_Trans: cnt = %d, %d.\n", A_Trans->cnt,
			p->A->cnt);
#endif

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	for (cnt = 1; cnt <= A_Trans->cnt; cnt++)
	{
		A_Trans->val[cnt] = p->A->val[cnt];
		A_Trans->row[cnt] = p->A->col[cnt];
		A_Trans->col[cnt] = p->A->row[cnt];
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "A:       cnt = %d, val = %f, row = %d, col = %d\n",
				cnt, p->A->val[cnt], p->A->row[cnt], p->A->col[cnt]);
		fprintf(g_FilePointer, "A_Trans: cnt = %d, val = %f, row = %d, col = %d\n",
				cnt, A_Trans->val[cnt], A_Trans->row[cnt], A_Trans->col[cnt]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "After A_Trans\n");
#endif

	/* 5. */

	/* Yifan 03/09/2012 Possible changes in optimality test for added feasibility Cuts*/
	/* Calculate - A_Trans * lambda. */
	TxX(A_Trans, lambda, A_Trans_lambda);

	/* Yifan 03/11/2012 Updated for Feasibility Cuts*/
	for (i = 0; i < p->num->mast_cols; i++)
	{
		for (cnt = 0; cnt < c->feasible_cuts_added->cnt; cnt++)
		{
			A_Trans_lambda[i + 1] -= c->feasible_cuts_added->val[cnt]->beta[i
					+ 1] * lambda[p->num->mast_rows + cnt + 1];
		}
	}

	/* Be careful!! Since we included the bound constraints, we also need the
	 adjustment on A_Trans * lambda here. */
	/* wrong !!!!
	 for (cnt=1; cnt<=A_Trans->cnt; cnt++) {
	 A_Trans_lambda[A_Trans->row[cnt]] += s->Master_dj[A_Trans->row[cnt]]; 
	 }
	 */

	/* Yifan 03/14/2012 Make sure addition or subtraction*/
	/* we need pi*A_j - c_j = Mu_L_j - Mu_U_j here */
	/* But A_trans_lambda is the negative of the value we do need*/
	for (i = 0; i < p->num->mast_cols; i++)
	{
		A_Trans_lambda[i + 1] += s->Master_dj[i + 1];
	}

	/* 
	 ** Calculate Vk_theta = Vk * theta.
	 */
	/*added by Yifan for debugging*/
	/*print_vect(s->Master_pi,p->num->mast_rows+p->num->max_cuts+1,"master_pi");*/

	for (cnt = 0; cnt < T->cnt; cnt++)
	{

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "alpha[%d] = %f\n", cnt, T->val[cnt]->alpha);
		for (idx=1; idx<=p->num->mast_cols; idx++)
		fprintf(g_FilePointer, "beta[%d, %d] = %f, ", cnt, idx,
				T->val[cnt]->beta[idx]);
		fprintf(g_FilePointer, "\n");
#endif

		/*
		 ** Obtain theta from s->Master_pi. Obtain Vk from 
		 **       T->val[i]->alpha - T->val[i]->beta * incumb_x. 
		 ** Note: be careful of the corresponding row number of each cut. And also 
		 ** be careful that the first element of s->Master_pi was reserved for one 
		 ** norm. 
		 */

		/* 6. */
		/* Yifan 03/14/2012 Make sure the index is correct*/
		/* Yifan 03/19/2012 Test removing eta0*/
		theta[cnt + 1] = ((double) c->k / (double) T->val[cnt]->cut_obs)
				* s->Master_pi[T->val[cnt]->row_num + 1];

#ifdef OPT_TEST
		printf("theta[%d]:%f\t",cnt+1,theta[cnt+1]);
#endif

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "theta[%d] = %f, rom_num = %d\n", cnt+1,
				theta[cnt+1], T->val[cnt]->row_num);
#endif

		/* 7. */

		/* Yifan 03/14/2012 Updated for optimality cut height*/
		Vk[cnt + 1] = T->val[cnt]->alpha
				- CxX(T->val[cnt]->beta, s->incumb_x, p->num->mast_cols);

		/* Yifan 05/03/13 Update Vk */
		/* modified by Yifan 2013.05.06 */
		/*Vk[cnt+1] = Vk[cnt+1] * (double) T->val[cnt]->cut_obs / (double) c->k;*/

		/* Yifan 04/08/2012 Update rhs including new lower bound */

		/* modified by Yifan 2013.05.06 */
		//if (sd_global->config.LB_TYPE ==1) {
		/* Vk[cnt+1] += ((double) c->k / (double) T->val[cnt]->cut_obs - 1) * sd_global->Eta0;*/
		//Vk[cnt+1] += (1 - (double) T->val[cnt]->cut_obs / (double) c->k ) * sd_global->Eta0;
		//}
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "Vk'[%d] = %f, ", cnt+1, Vk[cnt+1]);
#endif

		/* Also be careful to update alpha and beta based on the change of 
		 coefficients of eta. */
		/* Yifan 03/14/2012 Updated for optimality cut height*/

		/* NOTICE THE CALCULATION CHANGES!!!!! 03/29/2012*/
		/*Vk[cnt+1] *= (double)T->val[cnt]->cut_obs / (double)c->k;*/

		/*
		 if (sd_global->config.LB_TYPE == 1) {
		 Vk[cnt+1] +=(1- (double)T->val[cnt]->cut_obs / (double)c->k) * eta_zero;
		 }
		 else {
		 Vk[cnt+1] +=(1- (double)T->val[cnt]->cut_obs / (double)c->k) * sd_global->Eta0;
		 }*/

#ifdef OPT_TEST
		printf("Vk[%d]: %f\n",cnt+1,Vk[cnt+1]);
#endif

		/* Yifan 03/14/2012 Updated for optimality cut height*/
		// Vk[cnt+1] = cut_height(T->val[cnt], s->incumb_x, c, p->num);
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "Vk[%d] = %f\n", cnt+1, Vk[cnt+1]);
#endif
		/*added by Yifan for debugging*/
		/*printf("theta[%d] = %f, rom_num = %d \n", cnt+1, 
		 theta[cnt+1], T->val[cnt]->row_num);*/

	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n");
#endif

	/* 8. */
	/* Calculate Vk_theta. */
	Vk_theta = CxX(Vk, theta, T->cnt);

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "Vk_theta = %f\n", Vk_theta);
#endif

	/* 
	 ** Calculate Bk_theta. Bk is the transpose of the matrix of cut
	 ** coefficients {T->beta}, where T are the reformed cuts. 
	 */
	for (i = 1; i <= p->num->mast_cols; i++)
	{

		/* 9. */
		for (cnt = 0; cnt < T->cnt; cnt++)
		{

			/* Yifan 03/14/2012 migth need updates here!!!*/
			/* NOTICE THE CALCULATION CHANGES!!!!! 03/29/2012*/
			Bk_col[cnt + 1] = T->val[cnt]->beta[i];
			/*Bk_col[cnt+1] *= (double)T->val[cnt]->cut_obs / (double)c->k;*/

			/* Yifan 03/14/2012 Updated for optimality cut height*/

#ifdef CAL_CHECK
			fprintf(g_FilePointer, "Bk_col[%d] = %f, ", cnt+1, Bk_col[cnt+1]);
#endif
		}
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\n");
#endif

		/* 10. */
		Bk_theta[i] = CxX(Bk_col, theta, T->cnt);
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "Bk_theta[%d] = %f\n", i, Bk_theta[i]);
#endif
	}

	/* Calculate the quadratic vector, q_vec. 
	 Note: A_Trans_lambda = - A_Trans * lambda. */

	/* 11. */
	for (i = 1; i <= p->num->mast_cols; i++)
	{
		/*
		 q_vec[i] = Bk_theta[i] - A_Trans_lambda[i];
		 */
		q_vec[i] = p->c[i] - Bk_theta[i] - A_Trans_lambda[i];

#ifdef OPT_TEST
		printf("c: %f - B_t: %f - A_T_l: %f   = q_vec: %f\n", p->c[i], Bk_theta[i], A_Trans_lambda[i], q_vec[i]);
#endif
	}

	q_term = CxX(q_vec, q_vec, p->num->mast_cols);

#ifdef OPT_TEST
	printf("\n***  Vk_theta= %f\n",Vk_theta);
	printf("*** bk_lambda = %f\n",bk_lambda);
	printf("*** q_term = %f\n",q_term);
	printf("*** c->quad_scalar = %f\n\n",c->quad_scalar);
#endif

	/* update_quad_scalar(p, c, s, q_term, eta_zero);*/

	/*
	 ** Now it is the time to calculate the lower bound, Lm. 
	 */

	/* 12. */
	Lm = Vk_theta + bk_lambda - q_term / c->quad_scalar / 2.0;

	/*Lm = Vk_theta + bk_lambda - q_term / 2.0;*/

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "V_t: %f - b_l: %f - q_t: %f / %f / 2.0 = %f\n",
			Vk_theta, bk_lambda, q_term, c->quad_scalar, Lm);

	fprintf(g_FilePointer, "s->opt_value (d_obj) = %f\n", s->opt_value);
	fprintf(g_FilePointer, "c * incumb_x = %f\n", CxX(p->c, s->incumb_x,
					p->num->mast_cols));
	fprintf(g_FilePointer, "incumb_est = %f\n", s->incumb_est);

	fprintf(g_FilePointer, "Sm = %f, Lm = %f.\n", Sm, Lm);

	fprintf(g_FilePointer, "\n****** Leaving full_test()---cal_temp_lb ******\n");
#endif

	mem_free(bk);
	mem_free(lambda);
	mem_free(A_Trans_lambda);
	mem_free(A_Trans->col);
	mem_free(A_Trans->row);
	mem_free(A_Trans->val);
	mem_free(A_Trans);
	mem_free(theta);
	mem_free(Vk);
	mem_free(Bk_theta);
	mem_free(Bk_col);
	mem_free(q_vec);
	mem_free(lb);
	mem_free(ub);
	return Lm;
}

