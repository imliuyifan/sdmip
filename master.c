/***********************************************************************\
**
 ** master.c
 **
 ** This file contains the functions needed for solving the master
 ** problem in the SD algorithm, given the cuts and all the coefficients.
 **
 **
 **
 ** solve_master()
 ** solve_QP_master()   added for regularized QP method. zl 
 ** add_cut()
 ** change_cut()
 ** change_eta_col()
 **
 ** History:
 **   02 Feb 1992 - <Jason Mai> - created.
 **   11 Mar 1992 - <Jason Mai> - revised structures and parameters.
 **   20 Mar 1992 - <Jason Mai> - debugged.
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
#include "testout.h"
#include "log.h"
#include "master.h"
#include "sdglobal.h"

/***********************************************************************\
** This function adds the newest cut to the master problem, updates
 ** the incumbent cut if necessary, updates the coefficients on all
 ** the cuts, and finally solves the master problem.
 ** It returns TRUE if the problem was solved; FALSE otherwise.
 \***********************************************************************/
BOOL solve_master(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	BOOL ans;
	/*
	 int	best;
	 double  candid_est;
	 int   i; 
	 */
	clock_t start, end; /* Recording solution time for solving master LPs. 
	 added by zl, 06/29/04. */

#ifdef SAVE
	char fname[20] = "master   .lp";
	static int mnum = 0;
#endif

#ifdef TRACE
	printf("Inside solve_master\n");
#endif

	/* Update eta coefficient on all cuts, based on cut_obs */
	change_eta_col(c->master, c->cuts, c->k, s, p->num);

#ifdef SAVE
	fname[6] = '0' + mnum / 100 % 10;
	fname[7] = '0' + mnum / 10 % 10;
	fname[8] = '0' + mnum / 1 % 10;
	++mnum;
	print_problem(c->master, fname);
	printf("Saving file: %s\n", fname);
#endif

	/* Recording the time for solving master LPs. zl, 06/29/04. */
	start = clock();
	ans = solve_problem(sd_global, c->master);
	end = clock();
	s->run_time->soln_master_iter = ((double) (end - start)) / CLOCKS_PER_SEC;
	s->run_time->soln_master_accum += s->run_time->soln_master_iter;

	c->LP_cnt++; /* # of LPs solved increase by 1. zl 06/30/02 */

	/* Get the most recent optimal solution to master program */
	get_primal(s->candid_x, c->master, p->num->mast_cols);
	/* Get the dual solution too, JH 4/8/98 */

	/* Yifan 03/12/2012 Updated for Feasibility Cuts, make sure the total number is correct here*/
	get_dual(s->Master_pi, c->master, p->num,
			p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt);
	get_dual_slacks(s->Master_dj, c->master, p->num, p->num->mast_cols + 1);

	/** 
	 if (c->k == config.MAX_ITER)
	 { printf(" In solve_master \n");
	 for (i=0; i<p->num->mast_rows+c->cuts->cnt+1; i++)
	 printf(" Master_pi[%d] = %lf \n", i, s->Master_pi[i]); 
	 printf("\n"); 
	 for (i=0; i<p->num->mast_cols+2; i++)
	 printf(" Master_dj[%d] = %lf \n", i, s->Master_dj[i]); 
	 } 
	 **/

	s->candid_est = get_objective(c->master);

	/* Calculate gamma for next improvement check on incumbent x */
	s->gamma = s->candid_est - s->incumb_est;

	return ans;
}

/***********************************************************************\
** This function is the regularized QP version of solve_master(). It does 
 ** the same things as solve_master() does in the LP version. zl 
 ** It adds the newest cut to the master problem, updates
 ** the incumbent cut if necessary, updates the coefficients on all
 ** the cuts, and finally solves the master problem.
 ** It returns TRUE if the problem was solved; FALSE otherwise.
 \***********************************************************************/
BOOL solve_QP_master(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	BOOL ans;
//	double gamma;
	/*  int	best; */
//	double candid_est;
	double ht; /* height at the candidate solution. */
	double Sm;
	double d2 = 0.0;
	double eta[1];
	int i;
	clock_t start, end; /* Recording solution time for solving master QPs. 
	 added by zl, 06/29/04. */

#ifdef SAVE
	char fname[20] = "master   .lp";
	static int mnum = 0;
#endif

#ifdef TRACE
	printf("Inside solve_QP_master\n");
#endif

	/* In regularized QP method, we switch from dual simplex optimizer to
	 barrier optimizer for solving quadratic master problem. zl */

	change_solver_barrier(c->master);

	/* Update eta coefficient on all cuts, based on cut_obs */
	change_eta_col(c->master, c->cuts, c->k, s, p->num);

	if (sd_global->config.LB_TYPE == 1)
	{
		update_rhs(sd_global, p, c, s);
	}

#ifdef SAVE
	fname[6] = '0' + mnum / 100 % 10;
	fname[7] = '0' + mnum / 10 % 10;
	fname[8] = '0' + mnum / 1 % 10;
	++mnum;
	print_problem(c->master, fname);
	printf("Saving file: %s\n", fname);
#endif

	/* Recording the time for solving master QPs. zl, 06/29/04. */
	start = clock();
	ans = solve_problem(sd_global, c->master);
	end = clock();
	s->run_time->soln_master_iter = ((double) (end - start)) / CLOCKS_PER_SEC;
	s->run_time->soln_master_accum += s->run_time->soln_master_iter;

	c->LP_cnt++; /* # of LPs solved increase by 1. zl 06/30/02 */

	s->opt_value = get_objective(c->master);
	/*
	 fprintf(g_FilePointer, "\n****** In solve_QP_master ******\n");
	 */
	/* Get the most recent optimal solution to master program */
	get_primal(s->candid_x, c->master, p->num->mast_cols);

	if (s->opt_value > s->incumb_est)
	{
		//printf("!!!\n");
	}
#ifdef CAL_CHECK
	/* Print out d[] for checking purpose. */
	for (i=1; i<=p->num->mast_cols+1; i++)
	fprintf(g_FilePointer, "d[%d] = %f, ", i, s->candid_x[i]);
	fprintf(g_FilePointer, "\n");

	/* Print out incumb_x for checking purpose. */
	for (i=1; i<=p->num->mast_cols; i++)
	fprintf(g_FilePointer, "inc_x[%d] = %f, ", i, s->incumb_x[i]);
	fprintf(g_FilePointer, "\n");
#endif

	/* Since we are solving for 'd' now, we need to switch back to 'x'.
	 x = xbar + d.  zl */

	for (i = 1; i <= p->num->mast_cols; i++)
	{
		d2 += s->candid_x[i] * s->candid_x[i];
		s->candid_x[i] += s->incumb_x[i];
#ifdef CAL_CHECK
		fprintf(g_FilePointer, "can_x[%d] = %f, ", i, s->candid_x[i]);
#endif
	}
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n\nd2 = %f, c*inc_x = %f, d_obj = %f, c*can_x = %f",
			d2, CxX(p->c, s->incumb_x, p->num->mast_cols), s->opt_value,
			CxX(p->c, s->candid_x, p->num->mast_cols));
#endif

	/* update d_norm_k in soln_type. */
	if (c->k == 1)
		s->norm_d_k_1 = d2;
	s->norm_d_k = d2;

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\n***in master.c, solve_QP_master. ***\n");
	fprintf(g_FilePointer, "c->k = %d, d2 = %f, norm_d_k_1 = %f, norm_d_k = %f\n",
			c->k, d2, s->norm_d_k_1, s->norm_d_k);
#endif

	/* Calculating the one_norm of 'x'.  zl */
	s->candid_x[0] = one_norm(s->candid_x + 1, p->num->mast_cols);

	/* Get the dual solution too, JH 4/8/98 */
	/* Yifan 03/12/2012 Updated for Feasibility Cuts, make sure the total number is correct here*/
	get_dual(s->Master_pi, c->master, p->num,
			p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt);

	/* Yifan 03/19/2012 Test duals*/
#if 0
	printf("c->k:%d\n",c->k);

	printf("\nTESTING DUALS\n");
	for (i=0; i< p->num->mast_rows + c->cuts->cnt + c->feasible_cuts_added->cnt; i++)
	{
		printf("dual[%d]:%f\n",i+1, s->Master_pi[i+1]);
	}

	printf("\nTESTING PRIMALS\n");
	for (i=1; i<=p->num->mast_cols; i++)
	{
		printf("can_x[%d] = %f\n", i, s->candid_x[i]);
	}
	printf("\n");
#endif

	/*
	 for (i=0 ; i< p->num->mast_rows + c->cuts->cnt+c->feasible_cuts_added->cnt; i++) {
	 printf("s->Master_pi[%d]:%f\n",i,s->Master_pi[i]);
	 }
	 */

	get_dual_slacks(s->Master_dj, c->master, NULL, p->num->mast_cols + 1);

#ifdef CAL_CHECK
	write_prob(c->master, "quad_cut.lp");
#endif

	/* Find the highest cut at the candidate solution. 
	 ** where
	 ** cut_height = alpha - beta(xbar+d)  
	 */

	/*
	 printf("c->cuts->val->cnt is %d\n",c->cuts->cnt);
	 printf("c->cuts->val[0]->alpha is what???\n");
	 */

	if (c->cuts->cnt > 0)
	{

		/*modified by Yifan to avoid evaluating feasibility cut 02/26/2011*/
		if (c->cuts->val[0]->subfeaflag == TRUE)
		{
			//printf("c->cuts->val[0]->alpha is %f\n",c->cuts->val[0]->alpha);
			/*  Sm = c->cuts->val[0]->alpha - CxX(c->cuts->val[0]->beta, s->candid_x,
			 p->num->mast_cols);*/
			/* Yifan 03/14/2012 Updated for optimality cut height*/
			Sm = cut_height(sd_global, c->cuts->val[0], s->candid_x, c, p->num);

			/* Be careful to update alpha and beta based on the change of the 
			 coefficients of eta column.  zl */
			/* Sm *= (double)c->cuts->val[0]->cut_obs / (double)c->k;*/
		}
		else
			Sm = -INFBOUND;

#ifdef CAL_CHECK
		fprintf(g_FilePointer, "\nht[0] = %f, row_num = %d :: ", Sm,
				c->cuts->val[0]->row_num);
#endif

		for (i = 1; i < c->cuts->cnt; i++)
		{

			/*ht = c->cuts->val[i]->alpha - CxX(c->cuts->val[i]->beta, s->candid_x,
			 p->num->mast_cols);
			 ht *= (double)c->cuts->val[i]->cut_obs / (double)c->k;*/
			/* Yifan 03/14/2012 Updated for optimality cut height*/
			ht = cut_height(sd_global, c->cuts->val[i], s->candid_x, c, p->num);
#ifdef CAL_CHECK
			fprintf(g_FilePointer, "ht[%d] = %f, row_num = %d :: ", i, ht,
					c->cuts->val[i]->row_num);
#endif
			if (Sm < ht)
				Sm = ht;
		}
	}
	else
	{
		if (sd_global->config.LB_TYPE == 0)
		{
			Sm = 0.0;
		}
		else
		{
			Sm = sd_global->Eta0;
		}
	}

	get_x(c->master, eta, p->num->mast_cols, p->num->mast_cols); /* 2011.10.30 */

#ifdef CAL_CHECK
	fprintf(g_FilePointer, "\nSm = %f, eta = %f \n", Sm, eta[0]);
#endif

	/* s->candid_est = c(xbar+d) + max{cut_height} */

	/*s->candid_est = s->opt_value + CxX(p->c, s->incumb_x, p->num->mast_cols)- d2 * c->quad_scalar / 2.0;*/
	s->candid_est = Sm + CxX(p->c, s->candid_x, p->num->mast_cols);
	/*s->candid_est = s->opt_value;*/

	/* Calculate gamma for next improvement check on incumbent x */
	/* if it is not solved in opt mode, the gamma will not change*/
	/* gamma will be updated with 0.0 when incumbent sub is infeasibile*/
	if (c->opt_mode == TRUE)
	{
		s->gamma = s->candid_est - s->incumb_est;
	}

	/*
	 printf("Sm : candid_est = %f, incumb_est = %f, gamma = %f\n", 
	 s->candid_est, s->incumb_est, s->gamma);*/
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "Sm : candid_est = %f, incumb_est = %f, gamma = %f\n",
			s->candid_est, s->incumb_est, s->gamma);
#endif

	/* s->candid_est = d_obj + c*xbar - quad_scalar*d2/2. */
//	candid_est = s->opt_value + CxX(p->c, s->incumb_x, p->num->mast_cols)
//			- d2 * c->quad_scalar / 2.0;
	/* Calculate gamma for next improvement check on incumbent x */
//	gamma = candid_est - s->incumb_est;
#ifdef CAL_CHECK
	fprintf(g_FilePointer, "eta: candid_est = %f, incumb_est = %f, gamma = %f\n",
			candid_est, s->incumb_est, gamma);
#endif

	/* In regularized QP method, we need to switch back from barrier
	 optimizer to primal simplex optimizer each time after solving the 
	 quadratic master problem. zl */

	change_solver_primal(c->master);

	//mem_free(eta);  //Deleted by Yifan since _eta_ is not allocated above April 27 2011

#ifdef TRACE
	printf("Exiting solve_QP_master\n");
#endif

	return ans;
}

/***********************************************************************\
** This function performs the updates on all the coefficients of eta
 ** in the master problem constraint matrix.  During every iteration,
 ** each of the coefficieints on eta are increased, so that the effect
 ** of the cut on the objective function is decreased.
 **
 ** THIS ONLY WORKS FOR A SINGLE CELL !!!!!!!!   FIX IT !!!!!
 ** Multiply them properly, as you do in f_k in theta.c.
 **
 ** And what is going to happen when you start to drop cuts ???
 ** The location of the cut is no longer valid...
 **
 ** Each cut needs to memorize its row.  Then this function will
 ** loop through cuts, and change correctly... if no cuts are dropped...
 ** 
 ** Or each cut can memorize its row name (you'll have to give it one, 
 ** like "0" and "1" and "2"), then request the row number of this
 ** row name from the solver, then change the coefficient at that row.
 \***********************************************************************/
void change_eta_col(one_problem *p, cut_type *cuts, int k, soln_type *soln,
		num_type *num)
{
	double *coef;
	int c;
	int eta_col[1]; /* added by Yifan to change eta coefficient in objective function 08/26/11*/
	double lb_0_neginf[1];/*added by Yifan to change eta coefficient in objective function 08/26/11*/
	char lu[1];
	double eta_coef0[1];
	double eta_coef1[1];
	BOOL etaFlag = FALSE;

#ifdef TRACE
	printf("Inside change_eta_col\n");
#endif

	eta_col[0] = num->mast_cols;
	eta_coef0[0] = 0.0;
	eta_coef1[0] = 1.0;
	lu[0] = 'L';

	/*
	 ** Calculate an array of coefficients for the eta column.
	 */
	if (!(coef = arr_alloc(cuts->cnt, double)))
		err_msg("Allocation", "change_eta_col", "coef");

	for (c = 0; c < cuts->cnt; c++)
	{
		coef[cuts->val[c]->row_num - num->mast_rows] = (double) k
				/ (double) cuts->val[c]->cut_obs;
		etaFlag = TRUE;
	}

	/*
	 ** Change the eta column in the master program, which corresponds
	 ** to the eta variable, starting with the row after the A matrix
	 ** and ending with the row of the last cut.
	 */
	if (!(change_col(p, num->mast_cols, coef, num->mast_rows,
			num->mast_rows + cuts->cnt)))
	{
		print_contents(p, "contents.out");
		err_msg("change_col", "change_eta_col", "returned FALSE");
	}

	/*added by Yifan to update the coefficient of eta in the objective 08/24/2011*/
	if (etaFlag == TRUE)
	{
		change_objective(p, 1, eta_col, eta_coef1);
		lb_0_neginf[0] = -INFBOUND;
	}
	else
	{
		change_objective(p, 1, eta_col, eta_coef0);
		lb_0_neginf[0] = 0.0;
	}

	change_bound(p, 1, eta_col, lu, lb_0_neginf); /* 2011.10.30 */
	/* if feasibility cut is added without any general cut, eta's lower bd sholud be 0. Yifan 08/26/2011*/

	mem_free(coef);
}

/***********************************************************************\
** This function will allocate and initialize memory for a copy of the 
 ** master problem passed to it as _master_.  In addition to the original
 ** master, it will add a new column, called "eta", whose coefficients
 ** are all zero and whose cost is one.  It will also include all cuts
 ** which exist in the _cuts_ structure as additional rows in the problem.
 ** Finally, it will allocate enough room in all relevant problem arrays
 ** to hold up to _extra_cuts_ more constraints, assuming these constraints
 ** to be non-sparse (non-zero coefficeints in every column).  It returns
 ** a pointer to the newly created master problem.
 **
 ** Note how carefully the coefficients are placed in the matval array.
 ** CPLEX requires them to be ordered by columns, but we know we will
 ** be adding rows.  So, each segment of the array corresponding to a
 ** given column contains enough empty locations at its end to hold all
 ** present and future cuts (rows).  This way, when we add a row, CPLEX 
 ** does not have to reorder and shift all the coefficients down.
 \***********************************************************************/
one_problem *new_master(one_problem *master, cut_type *cuts, int extra_cuts,
		vector x_k)
{
	one_problem *copy;
	int r, i, j, idx, cnt, len;
	sd_long col_offset, row_offset; //modified by Yifan to avoid memory leaks July 26 2011
	char cut_name[NAME_SIZE] =
	{ "Old    " }; /* diff. from add_cut */
	char *q;

#ifdef TRACE
	printf("Inside new_master\n");
#endif

	if (!(copy = (one_problem *) mem_malloc (sizeof(one_problem))))
		err_msg("Allocation", "new_master", "copy");

	/* Initialize unused fields in a master copy's CPLEX data. */
	copy->mae = 0;
	copy->maesz = 0;
	copy->estorsz = 0;
	copy->rngcol = NULL;
	copy->nrowind = NULL;
	copy->etype = NULL;
	copy->ename = NULL;
	copy->estore = NULL;
	copy->enzbeg = NULL;
	copy->enzcnt = NULL;
	copy->enzind = NULL;
	copy->enzval = NULL;
	copy->dataname = NULL;
	copy->rhsname = NULL;
	copy->rngname = NULL;
	copy->bndname = NULL;

	/* Initialize dimensions of copy based on master and new cuts. */
	copy->matsz = master->matsz + (cuts->cnt + extra_cuts) * (master->mac + 1);
	copy->marsz = master->mar + cuts->cnt + extra_cuts;
	copy->mar = master->mar + cuts->cnt;
	copy->macsz = master->mac + 1;
	copy->mac = master->mac + 1;
	copy->cstorsz = master->cstorsz + NAME_SIZE;
	copy->rstorsz = master->rstorsz + NAME_SIZE * (cuts->cnt + extra_cuts);
	copy->objsen = master->objsen;

	/* Make all allocations of known sizes, as calculated above */
	if (!(copy->name = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "new_master", "copy->name");
	if (!(copy->objname = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "new_master", "copy->objname");
	if (!(copy->objx = arr_alloc(copy->macsz, double)))
		err_msg("Allocation", "new_master", "copy->objx");
	if (!(copy->bdl = arr_alloc(copy->macsz, double)))
		err_msg("Allocation", "new_master", "copy->bdl");
	if (!(copy->bdu = arr_alloc(copy->macsz, double)))
		err_msg("Allocation", "new_master", "copy->bdu");
	if (!(copy->rhsx = arr_alloc(copy->marsz, double)))
		err_msg("Allocation", "new_master", "copy->rhsx");
	if (!(copy->senx = arr_alloc(copy->marsz, char)))
		err_msg("Allocation", "new_master", "copy->senx");
    if (!(copy->typex = arr_alloc(copy->macsz, char)))
		err_msg("Allocation", "new_master", "copy->typex");
	if (!(copy->matbeg = arr_alloc(copy->macsz, int)))
		err_msg("Allocation", "new_master", "copy->matbeg");
	if (!(copy->matcnt = arr_alloc(copy->macsz, int)))
		err_msg("Allocation", "new_master", "copy->matcnt");
	if (!(copy->cname = arr_alloc(copy->macsz, string)))
		err_msg("Allocation", "new_master", "copy->cname");
	if (!(copy->cstore = arr_alloc(copy->cstorsz, char)))
		err_msg("Allocation", "new_master", "copy->cstore");
	if (!(copy->rname = arr_alloc(copy->marsz, string)))
		err_msg("Allocation", "new_master", "copy->rname");
	if (!(copy->rstore = arr_alloc(copy->rstorsz, char)))
		err_msg("Allocation", "new_master", "copy->rstore");
	if (!(copy->matval = arr_alloc(copy->matsz, double)))
		err_msg("Allocation", "new_master", "copy->matval");
	if (!(copy->matind = arr_alloc(copy->matsz, int)))
		err_msg("Allocation", "new_master", "copy->matind");

	/*
	 ** First copy information directly from the original master problem.
	 */

	/* Copy the master problem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	i = 0;
	for (q = master->cname[0]; q < master->cname[0] + master->cstorsz; q++)
		copy->cstore[i++] = *q;

	i = 0;
	for (q = master->rname[0]; q < master->rname[0] + master->rstorsz; q++)
		copy->rstore[i++] = *q;

	strcpy(copy->name, master->name);
	strcpy(copy->objname, master->objname);

	/* Calculate difference in pointers for master/copy row and column names */
	col_offset = copy->cstore - master->cname[0];
	row_offset = copy->rstore - master->rname[0];

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < master->mac; j++)
	{
		copy->objx[j] = master->objx[j];
		copy->bdu[j] = master->bdu[j];
		copy->bdl[j] = master->bdl[j];
		copy->cname[j] = master->cname[j] + col_offset;
		copy->matbeg[j] = cnt;
		copy->matcnt[j] = master->matcnt[j];
		for (idx = master->matbeg[j];
				idx < master->matbeg[j] + master->matcnt[j]; idx++)
		{
			copy->matval[cnt] = master->matval[idx];
			copy->matind[cnt] = master->matind[idx];
			cnt++;
		}
		cnt += cuts->cnt + extra_cuts;
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < master->mar; r++)
	{
		copy->rhsx[r] = master->rhsx[r];
		copy->senx[r] = master->senx[r];
		copy->rname[r] = master->rname[r] + row_offset;
	}

	/*
	 ** Initialize information for the extra column in the new master.
	 */

	strcpy(copy->cstore + master->cstorsz, "eta");
	copy->cname[master->mac] = copy->cstore + master->cstorsz;
	copy->objx[master->mac] = 1.0;
	/* Change from cplex70 to cplex81. zl 04/13/05. */
	copy->bdu[master->mac] = INFBOUND;
	copy->bdl[master->mac] = -INFBOUND;
	copy->matbeg[master->mac] = cnt;
	copy->matcnt[master->mac] = cuts->cnt;

	for (idx = 0; idx < cuts->cnt; idx++)
	{
		copy->matval[idx + cnt] = 1.0;
		copy->matind[idx + cnt] = master->mar + idx;
	}

	/*
	 ** Now copy information from the cuts into the new master problem.
	 */

	if (cuts->cnt)
	{
		/* Copy the constraint coefficients */
		for (j = 0; j < master->mac; j++)
		{
			cnt = copy->matbeg[j] + copy->matcnt[j];
			copy->matcnt[j] += cuts->cnt;

			for (i = 0; i < cuts->cnt; i++)
			{
				copy->matval[cnt + i] = cuts->val[i]->beta[j + 1];
				copy->matind[cnt + i] = master->mar + i;
			}
		}

		/* Give names to the cut constraints, and add rhs values */
		len = master->rstorsz;
		for (cnt = 0; cnt < cuts->cnt; cnt++)
		{
			cut_name[3] = '0' + cnt / 1000 % 10;
			cut_name[4] = '0' + cnt / 100 % 10;
			cut_name[5] = '0' + cnt / 10 % 10;
			cut_name[6] = '0' + cnt / 1 % 10;
			copy->rname[master->mar + cnt] = copy->rstore + len;
			strcpy(copy->rname[master->mar + cnt], cut_name);
			len += strlen(cut_name) + 1;
			copy->rhsx[master->mar + cnt] = cuts->val[cnt]->alpha;
			copy->senx[master->mar + cnt] = 'G';
		}
	}

#ifdef DEBUG
	print_contents(copy, "copy.out");
#endif

	/* Load the copy into CPLEX now */
	if (!(setup_problem(copy)))
		err_msg("Problem Setup", "new_master", "copy");

	//print_problem(copy, "copy.lp"); /* added by Yifan to test copy.mps structure */

#ifdef DEBUG
	print_problem(copy, "copy.mps");
#endif

	/*
	 ** We're done, and we have room for extra_cuts more constraints.
	 */
#ifdef TRACE
	printf("Exiting new_master\n");
#endif

	return copy;
}

one_problem *orig_new_master(one_problem *master, cut_type *cuts,
		int extra_cuts)
{
	one_problem *copy;
	int r, i, j, idx, cnt, len;
	sd_long col_offset, row_offset;
	char cut_name[NAME_SIZE] =
	{ "Old    " }; /* diff. from add_cut */
	char *q;

#ifdef TRACE
	printf("Inside orig_new_master\n");
#endif

	if (!(copy = (one_problem *) mem_malloc (sizeof(one_problem))))
		err_msg("Allocation", "orig_new_master", "copy");

	/* Initialize unused fields in a master copy's CPLEX data. */
	copy->mae = 0;
	copy->maesz = 0;
	copy->estorsz = 0;
	copy->rngcol = NULL;
	copy->nrowind = NULL;
	copy->etype = NULL;
	copy->ename = NULL;
	copy->estore = NULL;
	copy->enzbeg = NULL;
	copy->enzcnt = NULL;
	copy->enzind = NULL;
	copy->enzval = NULL;
	copy->dataname = NULL;
	copy->rhsname = NULL;
	copy->rngname = NULL;
	copy->bndname = NULL;

	/* Initialize dimensions of copy based on master and new cuts. */
	copy->matsz = master->matsz + (cuts->cnt + extra_cuts) * (master->mac + 1);
	copy->marsz = master->mar + cuts->cnt + extra_cuts;
	copy->mar = master->mar + cuts->cnt;
	copy->macsz = master->mac + 1;
	copy->mac = master->mac + 1;
	copy->cstorsz = master->cstorsz + NAME_SIZE;
	copy->rstorsz = master->rstorsz + NAME_SIZE * (cuts->cnt + extra_cuts);
	copy->objsen = master->objsen;

	/* Make all allocations of known sizes, as calculated above */
	if (!(copy->name = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "orig_new_master", "copy->name");
	if (!(copy->objname = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "orig_new_master", "copy->objname");
	if (!(copy->objx = arr_alloc(copy->macsz, double)))
		err_msg("Allocation", "orig_new_master", "copy->objx");
	if (!(copy->bdl = arr_alloc(copy->macsz, double)))
		err_msg("Allocation", "orig_new_master", "copy->bdl");
	if (!(copy->bdu = arr_alloc(copy->macsz, double)))
		err_msg("Allocation", "orig_new_master", "copy->bdu");
	if (!(copy->rhsx = arr_alloc(copy->marsz, double)))
		err_msg("Allocation", "orig_new_master", "copy->rhsx");
	if (!(copy->senx = arr_alloc(copy->marsz, char)))
		err_msg("Allocation", "orig_new_master", "copy->senx");
    if (!(copy->typex = arr_alloc(copy->macsz, char)))
		err_msg("Allocation", "new_master", "copy->typex");
	if (!(copy->matbeg = arr_alloc(copy->macsz, int)))
		err_msg("Allocation", "orig_new_master", "copy->matbeg");
	if (!(copy->matcnt = arr_alloc(copy->macsz, int)))
		err_msg("Allocation", "orig_new_master", "copy->matcnt");
	if (!(copy->cname = arr_alloc(copy->macsz, string)))
		err_msg("Allocation", "orig_new_master", "copy->cname");
	if (!(copy->cstore = arr_alloc(copy->cstorsz, char)))
		err_msg("Allocation", "orig_new_master", "copy->cstore");
	if (!(copy->rname = arr_alloc(copy->marsz, string)))
		err_msg("Allocation", "orig_new_master", "copy->rname");
	if (!(copy->rstore = arr_alloc(copy->rstorsz, char)))
		err_msg("Allocation", "orig_new_master", "copy->rstore");
	if (!(copy->matval = arr_alloc(copy->matsz, double)))
		err_msg("Allocation", "orig_new_master", "copy->matval");
	if (!(copy->matind = arr_alloc(copy->matsz, int)))
		err_msg("Allocation", "orig_new_master", "copy->matind");

	/*
	 ** First copy information directly from the original master problem.
	 */

	/* Copy the master problem's column and row names */
	/* Assume uninitialized elements are zero, or '\0', from calloc */
	i = 0;
	for (q = master->cname[0]; q < master->cname[0] + master->cstorsz; q++)
		copy->cstore[i++] = *q;
	i = 0;
	for (q = master->rname[0]; q < master->rname[0] + master->rstorsz; q++)
		copy->rstore[i++] = *q;

#ifdef TRACE
	printf("Inside orig_new_master3\n");
#endif

	strcpy(copy->name, master->name);
	strcpy(copy->objname, master->objname);

	/* Calculate difference in pointers for master/copy row and column names */
	col_offset = copy->cstore - master->cname[0];
	row_offset = copy->rstore - master->rname[0];

#ifdef TRACE
	printf("Inside orig_new_master4\n");
#endif

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < master->mac; j++)
	{
		copy->objx[j] = master->objx[j];
		copy->bdu[j] = master->bdu[j];
		copy->bdl[j] = master->bdl[j];
		copy->cname[j] = master->cname[j] + col_offset;

		copy->matbeg[j] = cnt;
		copy->matcnt[j] = master->matcnt[j];
		for (idx = master->matbeg[j];
				idx < master->matbeg[j] + master->matcnt[j]; idx++)
		{
			copy->matval[cnt] = master->matval[idx];
			copy->matind[cnt] = master->matind[idx];
			cnt++;
		}
		cnt += cuts->cnt + extra_cuts;
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < master->mar; r++)
	{
		copy->rhsx[r] = master->rhsx[r];
		copy->senx[r] = master->senx[r];
		copy->rname[r] = master->rname[r] + row_offset;
	}

#ifdef TRACE
	printf("Inside orig_new_master5\n");
#endif
	/*
	 ** Initialize information for the extra column in the new master.
	 */

	strcpy(copy->cstore + master->cstorsz, "eta");
	copy->cname[master->mac] = copy->cstore + master->cstorsz;
	copy->objx[master->mac] = 1.0;
	/* Change from cplex70 to cplex81. zl 04/13/05. */
	copy->bdu[master->mac] = INFBOUND;
	copy->bdl[master->mac] = -INFBOUND;
	copy->matbeg[master->mac] = cnt;
	copy->matcnt[master->mac] = cuts->cnt;
	for (idx = 0; idx < cuts->cnt; idx++)
	{
		copy->matval[idx + cnt] = 1.0;
		copy->matind[idx + cnt] = master->mar + idx;
	}
#ifdef TRACE
	printf("Inside orig_new_master6\n");
#endif

	/*
	 ** Now copy information from the cuts into the new master problem.
	 */

	if (cuts->cnt)
	{
		/* Copy the constraint coefficients */
		for (j = 0; j < master->mac; j++)
		{
			cnt = copy->matbeg[j] + copy->matcnt[j];
			copy->matcnt[j] += cuts->cnt;

			for (i = 0; i < cuts->cnt; i++)
			{
				copy->matval[cnt + i] = cuts->val[i]->beta[j + 1];
				copy->matind[cnt + i] = master->mar + i;
			}
		}

		/* Give names to the cut constraints, and add rhs values */
		len = master->rstorsz;
		for (cnt = 0; cnt < cuts->cnt; cnt++)
		{
			cut_name[3] = '0' + cnt / 1000 % 10;
			cut_name[4] = '0' + cnt / 100 % 10;
			cut_name[5] = '0' + cnt / 10 % 10;
			cut_name[6] = '0' + cnt / 1 % 10;
			copy->rname[master->mar + cnt] = copy->rstore + len;
			strcpy(copy->rname[master->mar + cnt], cut_name);
			len += strlen(cut_name) + 1;
			copy->rhsx[master->mar + cnt] = cuts->val[cnt]->alpha;
			copy->senx[master->mar + cnt] = 'G';
		}
	}

#ifdef TRACE
	printf("Inside orig_new_master7\n");
#endif

#ifdef DEBUG
	print_contents(copy, "copy.out");
#endif

	/* Load the copy into CPLEX now */
	if (!(setup_problem(copy)))
		err_msg("Problem Setup", "orig_new_master", "copy");

#ifdef DEBUG
	print_problem(copy, "copy.mps");
#endif
	print_problem(copy, "copy.mps");

	/*
	 ** We're done, and we have room for extra_cuts more constraints.
	 */
#ifdef TRACE
	printf("Exiting orig_new_master\n");
#endif

	return copy;
}

/***********************************************************************\
** The master being used by each cell is just a copy of the real
 ** master which is stored in prob.  Once a cell is done, its copy
 ** (which includes all the cell's cuts) may be freed.  This 
 ** function unloads the problem from the CPLEX and frees the memory.
 \***********************************************************************/
void free_master(one_problem *copy)
{

#ifdef TRACE
	printf("Inside free_master\n");
#endif

	remove_problem(copy);
	free_one_prob(copy);
}

void update_rhs(sdglobal_type* sd_global, prob_type *prob, cell_type *cell,
		soln_type *soln)
{
	int cnt;
	double *rhs;
	int *indices;

#ifdef TRACE
	printf("Inside change_rhs.\n");
#endif

	if (!(rhs = arr_alloc(cell->cuts->cnt, double)))
		err_msg("Allocation", "change_rhs", "rhs");
	if (!(indices = arr_alloc(cell->cuts->cnt, int)))
		err_msg("Allocation", "change_rhs", "indices");

	for (cnt = 0; cnt < cell->cuts->cnt; cnt++)
	{
		rhs[cnt] = cell->cuts->val[cnt]->alpha_incumb;
		rhs[cnt] += ((double) cell->k / (double) cell->cuts->val[cnt]->cut_obs
				- 1) * sd_global->Eta0;
		indices[cnt] = cell->cuts->val[cnt]->row_num;
	}

	/* Now we change the rhs of the master problem. */
	change_rhside(cell->master, cell->cuts->cnt, indices, rhs);
	mem_free(rhs);
	mem_free(indices);
}

