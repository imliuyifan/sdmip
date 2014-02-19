/***********************************************************************\
**
 ** prob.c
 **
 ** This file contains the routines needed to solve a stochastic LP 
 ** using SD, given the original combined problem (in CPLEX format).
 **
 **    solve_SD()
 **    new_prob()
 **    free_prob()
 **    free_one_prob()
 **
 ** History:
 **   16 Mar 1992 - <Jason Mai> - created.
 **   22 Mar 1992 - <Jason Mai> - debugged.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "supomega.h"
#include "utility.h"
#include "testout.h"
#include "log.h"
#include "solver.h"
#include "sdglobal.h"
#include <string.h>

/***********************************************************************\
** This function solves an LP using stochastic decomposition.  The
 ** original problem is separated into a master problem and a subproblem,
 ** and also into cells for parallel processing.  Each cell is then
 ** successively solved and married until there is but one cell remaining,
 ** and it is solved.
 \***********************************************************************/
void solve_SD(sdglobal_type* sd_global, one_problem *original, vector x_k,
		int num_rv, int num_cipher, int row, int col, char *fname, int batch_id)
{
	prob_type *prob;
	cell_type *cell;

#ifdef TRACE
	printf("Inside solve_SD\n");
#endif

	prob = new_prob(sd_global, original, num_rv, num_cipher, row, col);
	prob->current_batch_id = batch_id;
	init_param(sd_global, prob);
	cell = new_cell(sd_global, prob, 0);
#ifdef RUN
	if (batch_id == 0)
	{
		print_num(sd_global, prob->num);
	}
#endif

#ifdef SAVE
	print_contents(prob->master, "master.out");
	print_contents(prob->subprob, "subprob.out");
#endif
	printf("\n\n-------------------- Replication No.%d --------------------\n",
			batch_id);
	solve_cell(sd_global, cell, prob, x_k, fname);
	print_contents(prob->master, "final_master.lp");

	free_cell(cell, prob->num);
	free_prob(sd_global, prob);
}

/***********************************************************************\
** This function allocates memory for the fields of the problem
 ** data structure, separates the original problem into two stages,
 ** and initializes the rest of the fields in the problem.
 \***********************************************************************/
prob_type *new_prob(sdglobal_type* sd_global, one_problem *original, int num_rv,
		int num_cipher, int row, int col)
{
	prob_type *p;
	int r, c;

#ifdef TRACE
	printf("Inside new_prob\n");
#endif

	if (!(p = (prob_type *) mem_malloc (sizeof(prob_type))))
		err_msg("Allocation", "new_prob", "prob");

	if (!(p->master = (one_problem *) mem_malloc (sizeof(one_problem))))
		err_msg("Allocation", "new_prob", "prob->master");

	if (!(p->subprob = (one_problem *) mem_malloc (sizeof(one_problem))))
		err_msg("Allocation", "new_prob", "prob->subprob");

	if (!(p->Rbar = (sparse_vect *) mem_malloc (sizeof(sparse_vect))))
		err_msg("Allocation", "new_prob", "prob->Rbar");

	if (!(p->Tbar = (sparse_matrix *) mem_malloc (sizeof(sparse_matrix))))
		err_msg("Allocation", "new_prob", "prob->Tbar");

	/* In the regularized QP method, we need master's A matrix too. zl */
	if (sd_global->config.MASTER_TYPE == SDQP)
	{
		if (!(p->A = (sparse_matrix *) mem_malloc (sizeof(sparse_matrix))))
			err_msg("Allocation", "new_prob", "prob->A");
	}

	if (!(p->coord = (coord_type *) mem_malloc (sizeof(coord_type))))
		err_msg("Allocation", "new_prob", "prob->coord");

	if (!(p->num = (num_type *) mem_malloc (sizeof(num_type))))
		err_msg("Allocation", "new_prob", "prob->num");

	p->tau = sd_global->config.TAU;
	p->num->iter = sd_global->config.MAX_ITER;
	p->num->rv = num_rv;
	p->num->cipher = num_cipher;
	p->eval_seed = sd_global->config.EVAL_SEED1;

	/* 
	 ** Obtain the master and subproblem from the original LP.  This will
	 ** initialize num, Rbar, Tbar, and c, as well as master and subprob.
	 */
	decompose(sd_global, original, p, row, col);

	/* 
	 ** Initialize the coord structure according to the locations
	 ** of the random variables in omega.  Also initialize some
	 ** some fields of num, having to do with random rows and columns.
	 */
	if (!(p->coord->omega_col = arr_alloc(num_rv+1, int)))
		err_msg("Allocation", "new_prob", "p->coord->omega_col");
	p->coord->omega_col[0] = get_omega_col(sd_global, p->coord->omega_col + 1);

	if (!(p->coord->omega_row = arr_alloc(num_rv+1, int)))
		err_msg("Allocation", "new_prob", "p->coord->omega_row");
	p->coord->omega_row[0] = get_omega_row(sd_global, p->coord->omega_row + 1);

	/*
	 ** Subtract _row_ from each of the random elements' col coordinates
	 ** so that they begin at zero.  Then add one to both row and col coordinates
	 ** to obtain 1-based indices.  Also count the number of random variables
	 ** which occur in R and T, based on the column coordinate of each rv.
	 */
	for (r = 1; r <= num_rv; r++)
		p->coord->omega_row[r] += (1 - row);
	for (c = 1; c <= num_rv; c++)
		if (++p->coord->omega_col[c])
			++p->num->rv_T;
		else
			++p->num->rv_R;

	p->coord->delta_col = find_cols(num_rv, &p->num->rv_cols,
			p->coord->omega_col);
	p->coord->sigma_col = find_cols(p->Tbar->cnt, &p->num->nz_cols,
			p->Tbar->col);
	p->coord->lambda_row = find_rows(num_rv, &p->num->rv_rows,
			p->coord->omega_row);

	return p;
}

/***********************************************************************\
** This function frees the fields of the problem data structure.
 ** Note that it *does* free the arrays and other data associated with
 ** each of the fields.  Everything goes!
 \***********************************************************************/
void free_prob(sdglobal_type* sd_global, prob_type *p)
{

#ifdef TRACE
	printf("Inside free_prob\n");
#endif

	/* Free CPLEX data from both stages */
	free_one_prob(p->master);
	free_one_prob(p->subprob);

	/* Free arrays inside coord structure, then free it. */
	mem_free(p->coord->omega_row);
	mem_free(p->coord->omega_col);
	mem_free(p->coord->delta_col);
	mem_free(p->coord->lambda_row);
	mem_free(p->coord->sigma_col);
	mem_free(p->coord);

	/* Free the sparse vecotrs and matrices */
	mem_free(p->Tbar->row);
	mem_free(p->Tbar->col);
	mem_free(p->Tbar->val);
	mem_free(p->Tbar);
	mem_free(p->Rbar->row);
	mem_free(p->Rbar->val);
	mem_free(p->Rbar);

	/* In regularized QP method, we also need to free A matrix. zl */
	if (sd_global->config.MASTER_TYPE == SDQP)
	{
		mem_free(p->A->row);
		mem_free(p->A->col);
		mem_free(p->A->val);
		mem_free(p->A);
	}

	mem_free(p->num);
	mem_free(p->c);
	mem_free(p);
}

/***********************************************************************\
** This function frees all non-NULL the arrays in a one_problem 
 ** structure, which is used by CPLEX to solve problems.  It is assumed
 ** that the problem has already been unloaded!  Also, if arrays were left 
 ** un-initialized, there will be problems!
 **
 ** Move this into solver.c ???
 \***********************************************************************/
void free_one_prob(one_problem *p)
{

#ifdef TRACE
	printf("Inside free_one_prob\n");
#endif

	if (p)
	{
		if (p->name)
			mem_free(p->name);
		if (p->objx)
			mem_free(p->objx);
		if (p->rhsx)
			mem_free(p->rhsx);
		if (p->senx)
			mem_free(p->senx);
        if (p->typex)
			mem_free(p->typex);
		if (p->rngcol)
			mem_free(p->rngcol);
		if (p->nrowind)
			mem_free(p->nrowind);
		if (p->matbeg)
			mem_free(p->matbeg);
		if (p->matcnt)
			mem_free(p->matcnt);
		if (p->matind)
			mem_free(p->matind);
		if (p->matval)
			mem_free(p->matval);
		if (p->bdl)
			mem_free(p->bdl);
		if (p->bdu)
			mem_free(p->bdu);
		if (p->etype)
			mem_free(p->etype);
		if (p->enzbeg)
			mem_free(p->enzbeg);
		if (p->enzcnt)
			mem_free(p->enzcnt);
		if (p->enzind)
			mem_free(p->enzind);
		if (p->enzval)
			mem_free(p->enzval);
		if (p->dataname)
			mem_free(p->dataname);
		if (p->objname)
			mem_free(p->objname);
		if (p->rhsname)
			mem_free(p->rhsname);
		if (p->rngname)
			mem_free(p->rngname);
		if (p->bndname)
			mem_free(p->bndname);
		if (p->cname)
			mem_free(p->cname);
		if (p->cstore)
			mem_free(p->cstore);
		if (p->rname)
			mem_free(p->rname);
		if (p->rstore)
			mem_free(p->rstore);
		if (p->ename)
			mem_free(p->ename);
		if (p->estore)
			mem_free(p->estore);
		/*
		 printf("zl free 3\n");
		 */
		mem_free(p);
	}
}

void init_param(sdglobal_type* sd_global, prob_type *p)
{
	int multiplier;
	if (sd_global->config.MIN_ITER <= 0)
	{
		/* Min iteration number = 30 x (# of random variabls) + 3 x (# of first stage variables + 3)*/
		if (p->num->rv < 5)
			multiplier = 30;
		else if (p->num->rv < 10)
			multiplier = 20;
		else
			multiplier = 10;

		sd_global->config.MIN_ITER = multiplier * p->num->rv
				+ 3 * (p->num->mast_cols + 3);
	}
    
    /* modified by Yifan 2013.9.14 TEST */
    //sd_global->config.MIN_ITER = 1;

	if (sd_global->config.PI_EVAL_START <= 0)
	{
		sd_global->config.PI_EVAL_START =
				max(1, sd_global->config.MIN_ITER - sd_global->config.SCAN_LEN);
	}

}

/***********************************************************************\
** This function separates a single CPLEX LP into two stages, based on
 ** a (row,column) breakpoint supplied by the caller.  Elements of
 ** the cost vector less than _col_ become costs for the master
 ** while those greater become part of the subproblem.  Elements
 ** in the constraint matrix which are "north west" of _row_ and
 ** _col_ become the constraint matrix of the master program,
 ** while those "south east" of the breakpoint become the constraint
 ** matrix of the subproblem.  There should be no elements "north east"
 ** of the breakpoint.  Elements "south west" of the breakpoint are
 ** considered the T (technology) matrix.  Finally, right hand side
 ** elements above _row_ are for the master, and those below _row_
 ** form the R vector in the subproblem.
 **
 ** cost:                    c  |   g
 **                        -----+-------
 ** constraints:          [  A  |  NULL ] =  b
 **                        -----x------- - - - -
 **                       [  T  |   W   ] =  R
 **
 ** _row_ and _col_ represent the first row and column of the subproblem.
 **
 ** In addition to these elements, the master problem is expanded to
 ** include an extra column (the "eta" column) which will be used in
 ** all of the cut constraints.  It is also allocated with enough
 ** room to hold the maximum number of cuts allowed.
 **
 ** This function should probably be in solver.c, since it is so
 ** dependent upon the solver being used.
 **
 \***********************************************************************/
int decompose(sdglobal_type* sd_global, one_problem *orig, prob_type *p,
		int row, int col)
{
	char *q;
	int c, idx, r, i;
	sd_long m_offset, s_offset; /* modified by Yifan to avoid memory leaks July 26 2011 */
	one_problem *m, *s;
	int status;

#ifdef TRACE
	printf("Inside decompose\n");
#endif

	m = p->master;
	s = p->subprob;

	/* Initialize unused fields in master and subproblem CPLEX data. */
	m->mae = 0;
	s->mae = 0;
	m->maesz = 0;
	s->maesz = 0;
	m->enzsz = 0;
	s->enzsz = 0;
	m->estorsz = 0;
	s->estorsz = 0;
	m->rngcol = NULL;
	s->rngcol = NULL;
	m->nrowind = NULL;
	s->nrowind = NULL;
	m->etype = NULL;
	s->etype = NULL;
	m->ename = NULL;
	s->ename = NULL;
	m->estore = NULL;
	s->estore = NULL;
	m->enzbeg = NULL;
	s->enzbeg = NULL;
	m->enzcnt = NULL;
	s->enzcnt = NULL;
	m->enzind = NULL;
	s->enzind = NULL;
	m->enzval = NULL;
	s->enzval = NULL;
	m->dataname = NULL;
	s->dataname = NULL;
	m->rhsname = NULL;
	s->rhsname = NULL;
	m->rngname = NULL;
	s->rngname = NULL;
	m->bndname = NULL;
	s->bndname = NULL;

	/* Initialize dimensions of each stage based on original data. */
	m->matsz = s->matsz = 0;
	m->mar = m->marsz = row;
	m->mac = m->macsz = col;
	s->mac = s->macsz = orig->mac - col;
	s->mar = s->marsz = orig->mar - row;
	m->cstorsz = (unsigned int) (orig->cname[col] - orig->cname[0]);
	s->cstorsz = orig->cstorsz - m->cstorsz;
	m->rstorsz = (unsigned int) (orig->rname[row] - orig->rname[0]);
	s->rstorsz = orig->rstorsz - m->rstorsz;
	m->objsen = s->objsen = orig->objsen;

	/*
	 printf("zl m: mar = %d, mac = %d, cstorsz = %d, rstorsz = %d, objsen = %d\n",
	 m->mar, m->mac, m->cstorsz, m->rstorsz, m->objsen);
	 printf("zl s: mar = %d, mac = %d, cstorsz = %d, rstorsz = %d, objsen = %d\n",
	 s->mar, s->mac, s->cstorsz, s->rstorsz, s->objsen);
	 */

	p->Tbar->cnt = 0;
	p->Rbar->cnt = 0;
	p->num->rv_R = 0;
	p->num->rv_T = 0;
	p->num->mast_rows = row;
	p->num->mast_cols = col;
	p->num->sub_rows = s->mar;
	p->num->sub_cols = s->mac;
	if (sd_global->config.MASTER_TYPE == SDQP)
		p->A->cnt = 0; /* Master's A matrix. zl */

	/* sd_global->config.MASTER_TYPE indicates whether we are solving the master problem
	 using basic LP method or regularized QP method. If QP, then we only need to     keep m->mac+3 cuts. zl */

	if (sd_global->config.MASTER_TYPE == SDLP)
		p->num->max_cuts = sd_global->config.CUT_MULT * MAX_CUTS(m->mac);
	else
	{
		sd_global->config.CUT_MULT = 1;
		p->num->max_cuts = sd_global->config.CUT_MULT * m->mac + 3;
	}
	/*
	 fprintf(g_FilePointer,
	 "sd_global->config.CUT_MULT = %d, m->mac = %d, p->num->max_cuts = %d. \n\n",
	 sd_global->config.CUT_MULT, m->mac, p->num->max_cuts);
	 */

#ifdef DEBUG
	printf("p->num->max_cuts = %d \n\n", p->num->max_cuts);
#endif

	/*
	 ** Make initial allocations of maximum foreseeable size
	 ** for all special data in the prob_type structure.
	 ** Worst case scenario: entire problem is subproblem.
	 ** Recall that Tbar, Rbar, and c will require 1-norms.
	 */
	if (!(p->Tbar->row = arr_alloc(orig->matsz+1, int)))
		err_msg("Allocation", "decompose", "Tbar->row");
	if (!(p->Tbar->col = arr_alloc(orig->matsz+1, int)))
		err_msg("Allocation", "decompose", "Tbar->col");
	if (!(p->Tbar->val = arr_alloc(orig->matsz+1, double)))
		err_msg("Allocation", "decompose", "Tbar->val");
	if (!(p->Rbar->row = arr_alloc(orig->marsz+1, int)))
		err_msg("Allocation", "decompose", "Rbar->row");
	if (!(p->Rbar->val = arr_alloc(orig->marsz+1, double)))
		err_msg("Allocation", "decompose", "Rbar->val");
	if (!(p->c = arr_alloc(orig->macsz+1, double)))
		err_msg("Allocation", "decompose", "p->c");

	/*
	 ** In the case we are using regularized QP method, we also need to make
	 ** the initial allocation of maximum foreseeable size for A matrix.
	 ** Worst case scenario: entire problem is master problem.  zl
	 */
	if (sd_global->config.MASTER_TYPE == SDQP)
	{
		if (!(p->A->row = arr_alloc(orig->matsz+1, int)))
			err_msg("Allocation", "decompose", "A->row");
		if (!(p->A->col = arr_alloc(orig->matsz+1, int)))
			err_msg("Allocation", "decompose", "A->col");
		if (!(p->A->val = arr_alloc(orig->matsz+1, double)))
			err_msg("Allocaiton", "decompose", "A->val");
	}

	/*
	 ** Make initial allocations of known sizes (based on row
	 ** and col) for master and subproblem CPLEX data.
	 */
	if (!(m->name = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "decompose", "m->name");
	if (!(s->name = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "decompose", "s->name");
	if (!(m->objname = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "decompose", "m->objname");
	if (!(s->objname = arr_alloc(NAME_SIZE, char)))
		err_msg("Allocation", "decompose", "s->objname");
	if (!(m->objx = arr_alloc(m->macsz, double)))
		err_msg("Allocation", "decompose", "m->objx");
	if (!(s->objx = arr_alloc(s->macsz, double)))
		err_msg("Allocation", "decompose", "s->objx");
	if (!(m->bdl = arr_alloc(m->macsz, double)))
		err_msg("Allocation", "decompose", "m->bdl");
	if (!(s->bdl = arr_alloc(s->macsz, double)))
		err_msg("Allocation", "decompose", "s->bdl");
	if (!(m->bdu = arr_alloc(m->macsz, double)))
		err_msg("Allocation", "decompose", "m->bdu");
	if (!(s->bdu = arr_alloc(s->macsz, double)))
		err_msg("Allocation", "decompose", "s->bdu");
	if (!(m->rhsx = arr_alloc(m->marsz, double)))
		err_msg("Allocation", "decompose", "m->rhsx");
	if (!(s->rhsx = arr_alloc(s->marsz, double)))
		err_msg("Allocation", "decompose", "s->rhsx");
	if (!(m->senx = arr_alloc(m->marsz, char)))
		err_msg("Allocation", "decompose", "m->senx");
	if (!(s->senx = arr_alloc(s->marsz, char)))
		err_msg("Allocation", "decompose", "s->senx");
    if (!(m->typex = arr_alloc(m->macsz, char)))                /* Yifan 2013.09.12 */
		err_msg("Allocation", "decompose", "m->typex");
    if (!(s->typex = arr_alloc((s->macsz), char)))                /* Yifan 2013.09.12 */
		err_msg("Allocation", "decompose", "s->typex");
	if (!(m->matbeg = arr_alloc(m->macsz, int)))
		err_msg("Allocation", "decompose", "m->matbeg");
	if (!(s->matbeg = arr_alloc(s->macsz, int)))
		err_msg("Allocation", "decompose", "s->matbeg");
	if (!(m->matcnt = arr_alloc(m->macsz, int)))
		err_msg("Allocation", "decompose", "m->matcnt");
	if (!(s->matcnt = arr_alloc(s->macsz, int)))
		err_msg("Allocation", "decompose", "s->matcnt");
	if (!(m->cname = arr_alloc(m->macsz, string)))
		err_msg("Allocation", "decompose", "m->cname");
	if (!(s->cname = arr_alloc(s->macsz, string)))
		err_msg("Allocation", "decompose", "s->cname");
	if (!(m->cstore = arr_alloc(m->cstorsz, char)))
		err_msg("Allocation", "decompose", "m->cstore");
	if (!(s->cstore = arr_alloc(s->cstorsz, char)))
		err_msg("Allocation", "decompose", "s->cstore");
	if (!(m->rname = arr_alloc(m->marsz, string)))
		err_msg("Allocation", "decompose", "m->rname");
	if (!(s->rname = arr_alloc(s->marsz, string)))
		err_msg("Allocation", "decompose", "s->rname");
	if (!(m->rstore = arr_alloc(m->rstorsz, char)))
		err_msg("Allocation", "decompose", "m->rstore");
	if (!(s->rstore = arr_alloc(s->rstorsz, char)))
		err_msg("Allocation", "decompose", "s->rstore");

	/*
	 ** Make intial allocations of maximum foreseeable sizes
	 ** for the master and subproblem CPLEX data.  Again, in the
	 ** worst case, could be all master, or all subproblem.
	 */
	if (!(m->matval = arr_alloc(orig->matsz, double)))
		err_msg("Allocation", "decompose", "m->matval");
	if (!(s->matval = arr_alloc(orig->matsz, double)))
		err_msg("Allocation", "decompose", "s->matval");
	if (!(m->matind = arr_alloc(orig->matsz, int)))
		err_msg("Allocation", "decompose", "m->matind");
	if (!(s->matind = arr_alloc(orig->matsz, int)))
		err_msg("Allocation", "decompose", "s->matind");

	/*
	 ** Copy the names of each constraint and decision variable for
	 ** both the master and the subproblem.
	 */
	i = 0;
	for (q = orig->cname[0]; q < orig->cname[col]; q++)
		m->cstore[i++] = *q;
	i = 0;
	for (q = orig->cname[col]; q < orig->cname[col] + s->cstorsz; q++)
		s->cstore[i++] = *q;

	/*
	 ** Calculate offsets for name pointers, to be used in transferring
	 ** pointers to names in orignal, over to pointers to names in both
	 ** the subproblem and master.
	 */
	m_offset = m->cstore - orig->cname[0];
	s_offset = s->cstore - orig->cname[col];

	/*
	 ** Loop through all the coefficients in the "western" half of
	 ** the original constraint matrix, and separate A from T
	 ** Also pull out any data associated with master decision
	 ** variables, like cost coefficients and upper/lower bounds.
	 */
	for (c = 0; c < col; c++)
	{
		p->c[c + 1] = orig->objx[c];
		m->cname[c] = orig->cname[c] + m_offset;
		m->objx[c] = orig->objx[c];
        m->typex[c] = orig->typex[c];       /* modified by Yifan 2013.9.12 */
		m->bdl[c] = orig->bdl[c];
		m->bdu[c] = orig->bdu[c];
		m->matbeg[c] = m->matsz;
		m->matcnt[c] = 0;

		for (idx = orig->matbeg[c]; idx < orig->matbeg[c] + orig->matcnt[c];
				idx++)
		{
			if (orig->matind[idx] < row)
			{
				/* The coefficient is part of the master constraint matrix */
				m->matval[m->matsz] = orig->matval[idx];
				m->matind[m->matsz] = orig->matind[idx];
				++m->matcnt[c];
				++m->matsz;

				/* In the regularized QP method, we need to store the master's A
				 matrix. zl */
				if (sd_global->config.MASTER_TYPE == SDQP)
				{
					p->A->val[p->A->cnt + 1] = orig->matval[idx];
					p->A->row[p->A->cnt + 1] = orig->matind[idx] + 1;
					p->A->col[p->A->cnt + 1] = c + 1;
					/*
					 fprintf (g_FilePointer,
					 "A->col[%d] = %d, A->row[%d] = %d, A->val[%d] = %f\n",
					 p->A->cnt+1, p->A->col[p->A->cnt+1],
					 p->A->cnt+1, p->A->row[p->A->cnt+1],
					 p->A->cnt+1, p->A->val[p->A->cnt+1]);
					 */
					++p->A->cnt;
				}
			}
			else
			{
				/* The coefficient is part of the subproblem's T matrix */
				p->Tbar->val[p->Tbar->cnt + 1] = orig->matval[idx];
				p->Tbar->row[p->Tbar->cnt + 1] = orig->matind[idx] - row + 1;
				p->Tbar->col[p->Tbar->cnt + 1] = c + 1;
				++p->Tbar->cnt;
			}
		}
	}
    
    m->typex[col] = 'C';       /* modified by Yifan 2013.9.12 This is for the eta column */


	/* Watch out!!!!  You assume, in find_cols, that Tbar is sorted by col !!!! */

	/*
	 ** Now loop through all the coefficients in the "eastern" half of
	 ** the original constraint matrix, pulling out elements of W.
	 ** Also pull out any data associated with subproblem decision
	 ** variables, like cost coefficients and upper/lower bounds.
	 */
	for (c = col; c < orig->mac; c++)
	{
		s->cname[c - col] = orig->cname[c] + s_offset;
		s->objx[c - col] = orig->objx[c];
        s->typex[c - col] = orig->typex[c];    /* modified by Yifan 2013.9.12 */
		s->bdu[c - col] = orig->bdu[c];
		s->bdl[c - col] = orig->bdl[c];
		s->matbeg[c - col] = s->matsz;
		s->matcnt[c - col] = 0;
		for (idx = orig->matbeg[c]; idx < orig->matbeg[c] + orig->matcnt[c];
				idx++)
		{
			if (orig->matind[idx] < row)
			{
				/* The coefficient is not valid! */
				printf(
						"Error: constraint coefficient exists in invalid quadrant\n");
				return FALSE;
			}
			else
			{
				/* The coefficient is part of the subproblem constraint matrix */
				s->matval[s->matsz] = orig->matval[idx];
				s->matind[s->matsz] = orig->matind[idx] - row;
				++s->matcnt[c - col];
				++s->matsz;
			}
		}
	}

	/*
	 ** Copy the names of each constraint row for
	 ** both the master and the subproblem.
	 */
	i = 0;
	for (q = orig->rname[0]; q < orig->rname[row]; q++)
		m->rstore[i++] = *q;
	i = 0;
	for (q = orig->rname[row]; q < orig->rname[row] + s->rstorsz; q++)
		s->rstore[i++] = *q;

	/*
	 ** Calculate offsets for the rname arrays... when added to a pointer
	 ** in rname for the original, these give a new pointer for rname in
	 ** the master and subproblem.
	 */
	m_offset = m->rstore - orig->rname[0];
	s_offset = s->rstore - orig->rname[row];

	/*
	 ** Now go through the "north" and "south" halves of the original
	 ** problem, to find the right hand side and sense of each constraint
	 ** in the master and subproblem matrices.
	 */
	for (r = 0; r < row; r++)
	{
		m->rhsx[r] = orig->rhsx[r];
		m->senx[r] = orig->senx[r];
		m->rname[r] = orig->rname[r] + m_offset;
	}

	/* Initialize the subproblem rhs and Rbar with the same values */
	for (r = row; r < orig->mar; r++)
	{
		s->rhsx[r - row] = orig->rhsx[r];
		s->senx[r - row] = orig->senx[r];
		s->rname[r - row] = orig->rname[r] + s_offset;
		p->Rbar->val[p->Rbar->cnt + 1] = orig->rhsx[r];
		p->Rbar->row[p->Rbar->cnt + 1] = r - row + 1;
		++p->Rbar->cnt;
	}

	strcpy(m->name, "Master");
	strcpy(s->name, "Subproblem");
	strcpy(m->objname, orig->objname);
	strcpy(s->objname, orig->objname);

	/*
	 ** Our first goal is to shrink down the arrays in the prob, master,
	 ** and subproblem which were previously allocated to the maximum
	 ** foreseeable size, now that we know how large they must be.
	 ** (Recall that elements of the prob_type structure must have room
	 ** for their 1-norms).  But, for the master problem we must reserve
	 ** enough room for all of the cuts, as well as room for the extra
	 ** eta column.  Hence, the arrays will get even bigger.
	 */
	m->matval = (double *) mem_realloc (m->matval, m->matsz*sizeof(double));
	s->matval = (double *) mem_realloc (s->matval, s->matsz*sizeof(double));
	m->matind = (int *) mem_realloc (m->matind, m->matsz*sizeof(int));
	s->matind = (int *) mem_realloc (s->matind, s->matsz*sizeof(int));
	p->c = (double *) mem_realloc (p->c, (m->macsz+1)*sizeof(double));
	p->Tbar->row =
			(int *) mem_realloc (p->Tbar->row, (p->Tbar->cnt+1)*sizeof(int));
	p->Tbar->col =
			(int *) mem_realloc (p->Tbar->col, (p->Tbar->cnt+1)*sizeof(int));
	p->Tbar->val = (double *) mem_realloc (p->Tbar->val,
			(p->Tbar->cnt+1)*sizeof(double));
	p->Rbar->row =
			(int *) mem_realloc (p->Rbar->row, (p->Rbar->cnt+1)*sizeof(int));
	p->Rbar->val = (double *) mem_realloc (p->Rbar->val,
			(p->Rbar->cnt+1)*sizeof(double));

	/* While in regularized QP method, we need to reallocate the space for
	 master's A matrix too.  zl */

	if (SDQP == sd_global->config.MASTER_TYPE)
	{
		p->A->row = (int *) mem_realloc (p->A->row, (p->A->cnt+1)*sizeof(int));
		p->A->col = (int *) mem_realloc (p->A->col, (p->A->cnt+1)*sizeof(int));
		p->A->val = (double *) mem_realloc (p->A->val,
				(p->A->cnt+1)*sizeof(double));
	}

	/* The Cplex default for the number nonzero read limit in Q matrix is 500.
	 * When
	 * 		master->ncols >= 500,
	 * set it to be master->ncols + 1 (including the _eta_ column).
	 * zl, 09/26/05 */
	if ((SDQP == sd_global->config.MASTER_TYPE) && (m->mac >= 500))
	{
		printf("in decompose.c, set qp_nzreadlim to %d\n", m->mac + 1);
		status = set_qp_nzreadlim(m->mac + 1);
		if (0 != status)
		{
			printf("error: status = %d\n", status);
			exit(1);
		}
	}

	/*
	 ** Finally get rid of the original problem
	 */

#ifdef DEBUG
	printf("Printing in decompose\n");
	print_contents(m, "m.out");
	print_contents(s, "s.out");
	printf("Done in decompose\n");
#endif

	return TRUE;
}

/****************************************************************************\
** This function parses the command line paramenters and prompts the user
 ** for startup information required by SD.  It will determine how many
 ** problems will be solved (one or many) based on the command line arguments.
 ** A lack of arguments means that only one problem should be solved, and
 ** its name will be entered.  A number on the command line means that that
 ** many problems will be solved, with filenames "prob0000", "prob0001", etc.
 ** The second command line parameter specifies the starting point for a suite
 ** of test problems.
 \****************************************************************************/
void parse_cmd_line(sdglobal_type* sd_global, int argc, char *argv[],
		char *fname, int *objsen, int *num_probs, int *start, BOOL *read_seeds,
		BOOL *read_iters)
{
  int status;
	if (argc < 2)
	{
		printf("Please enter the name of the problem files (eg. `exags1p'): ");
		status = scanf("%s", fname);
        if (status <= 0) {
          printf("You entered nothing.\n");
        }
		/*
		 printf("Please enter the sense of the objective function (max=-1, min=1)");
		 scanf("%d", objsen);*/
		*objsen = 1;
		*num_probs = 1;
		*read_seeds = TRUE;
		*read_iters = TRUE;
	}
	else if (argc < 3)
	{
		strcpy(fname, "prob   ");
		*num_probs = atoi(argv[1]);
		*objsen = 1;
		*start = 0;
		*read_seeds = TRUE;
		*read_iters = TRUE; /* added by zl. 06/18/02. */
	}
	else if (argc == 5)
	{
		strcpy(fname, argv[1]);
		sscanf(argv[2], "%d", objsen);
		sscanf(argv[3], "%lld", &(sd_global->config.RUN_SEED1));
		sscanf(argv[3], "%lld", &(sd_global->config.RUN_SEED2));
		sscanf(argv[4], "%lld", &(sd_global->config.EVAL_SEED1));

		printf(
				"In sd.c, parse_cmd_line, argc = %d, objsen = %d, RUN_SEED1 = %lld, RUN_SEED2 = %lld, EVAL_SEED1 = %lld\n",
				argc, *objsen, sd_global->config.RUN_SEED1,
				sd_global->config.RUN_SEED2, sd_global->config.EVAL_SEED1);
		*num_probs = 1;
		*read_seeds = FALSE;
		*read_iters = TRUE; /* added by zl. 06/18/02. */
	}
	else if (argc == 7) /* added by zl. 06/18/02. */
	{
		strcpy(fname, argv[1]);
		sscanf(argv[2], "%d", objsen);
		sscanf(argv[3], "%lld", &(sd_global->config.RUN_SEED1));
		sscanf(argv[3], "%lld", &(sd_global->config.RUN_SEED2));
		sscanf(argv[4], "%lld", &(sd_global->config.EVAL_SEED1));
		sscanf(argv[5], "%d", &(sd_global->config.MIN_ITER));
		sscanf(argv[6], "%d", &(sd_global->config.MAX_ITER));
		sd_global->config.START_THIN = sd_global->config.MAX_ITER;
		*num_probs = 1;
		*read_seeds = FALSE;
		*read_iters = FALSE; /* added by zl. 06/18/02. */
		printf(
				"In sd.c, parse_cmd_line, argc = %d, objsen = %d, RUN_SEED1 = %lld, RUN_SEED2 = %lld, EVAL_SEED1 = %lld, MIN = %d, MAX = %d\n",
				argc, *objsen, sd_global->config.RUN_SEED1,
				sd_global->config.RUN_SEED2, sd_global->config.EVAL_SEED1,
				sd_global->config.MIN_ITER, sd_global->config.MAX_ITER);
	}
	else
	{
		strcpy(fname, "prob    ");
		*num_probs = atoi(argv[1]);
		*start = atoi(argv[2]);
		*read_seeds = TRUE;
		*read_iters = TRUE; /* added by zl. 06/18/02. */
		*objsen = 1;
	}
}

/****************************************************************************\
** This function informs the user of an error, and aborts the program.
 ** The first parameter is a string which describes what type of error
 ** has occurred (e.g. Allocation, NULL pointer, etc.)  The second
 ** parameter gives the function in which the error was recognized
 ** (i.e. the function which is calling this one).  The third parameter
 ** specifies which item / calculation / segment of the function caused
 ** the error (e.g. which array was not allocated, which pointer was NULL)
 **
 ** You can change this so that it asks the user if they dare to go on.
 ** Just return if they say yes... (heh heh heh).
 \****************************************************************************/
void err_msg(char *type, char *place, char *item)
{
	printf("\n\n||| %s error in function %s(), item %s.\n", type, place, item);
	exit(1);
}

void get_lower_bound(sdglobal_type* sd_global, one_problem *p, int row, int col)
{
	int ans, cnt;
//	int scr_stat;
	int *indices;
	double *values;
	double lower_bound;
	one_problem *copy;

	if (!(copy = (one_problem *) mem_malloc (sizeof(one_problem))))
		err_msg("Allocation", "get_lower_bound", "copy");
	if (!(sd_global->Bbar = arr_alloc(col+3, double)))
		err_msg("Allocation", "sd", "Bbar");

	/* Before any changes, clone a whole new problem from the master */
	copy->lp = clone_prob(p);

	/* Initialize memory space for indices */
	if (!(indices = arr_alloc(p->mac, int)))
		err_msg("Allocation", "get_lower_bound", "indices");

	/* Initialize memory space for values */
	if (!(values = arr_alloc(p->mac, double)))
		err_msg("Allocation", "get_lower_bound", "values");

	/* First get Beta and Alpha from the mean value problem*/
	calc_alpha_beta(sd_global, p, row, col);

	/* Then use Beta and Alpha to get the lower bound for optimality cut */
	for (cnt = 0; cnt < col; cnt++)
	{
		indices[cnt] = cnt;
		values[cnt] = -sd_global->Bbar[cnt + 1];
	}

	for (cnt = col; cnt < p->mac; cnt++)
	{
		indices[cnt] = cnt;
		values[cnt] = 0.0;
	}

	/*
     write_prob(copy, "before_coef_change.lp");
	 */
	

	/* Change the objective to lower bound calculation */
	change_objective(copy, p->mac, indices, values);

	/* write_prob(copy, "after_coef_change.lp");*/
	/*
	 write_prob(p, "SDprob_after_coef_change.lp");
	 */

	change_solver_primal(p);

	set_intparam(p, PARAM_SCRIND, ON);
	ans = solve_lp(copy);
	if (ans != 0)
	{
		printf("Solution status: %d.\n", ans);
		err_msg("solve_lp", "get_lower_bound",
				"Make sure that the first stage is bounded.");
	}
	set_intparam(p, PARAM_SCRIND, ON);

	lower_bound = get_objective(copy);
	sd_global->Eta0 = lower_bound + sd_global->Abar;

	printf("lower_bound is %f\n", lower_bound);
	printf("Eta0 is %f\n", sd_global->Eta0);

	remove_problem(copy);
	mem_free(copy);
	mem_free(indices);
	mem_free(values);
	mem_free(sd_global->Bbar);
}

/****************************************************************************\
****  added by Yifan Mar 29 2011 * get the dual solution, Tbar, Rbar, alpha, beta  **
 \****************************************************************************/
void calc_alpha_beta(sdglobal_type* sd_global, one_problem *p, int row, int col)
{
	double *y_i; /* Initial dual candidate, from original problem */
	int r, c, idx, j;
	prob_type *prob;

	/* Get the dual multipliers from the mean problem */
	if (!(y_i = arr_alloc(p->mar+1, double)))
		err_msg("Allocation", "main", "y_i");
	get_dual(y_i, p, NULL, p->mar);

	/*
	 ** Make initial allocations of size for prob_type
	 ** structure used in calculating Abar and Bbar
	 */
	if (!(prob = (prob_type *) mem_malloc (sizeof(prob_type))))
		err_msg("Allocation", "main", "prob");
	if (!(prob->Rbar = (sparse_vect *) mem_malloc (sizeof(sparse_vect))))
		err_msg("Allocation", "main", "prob->Rbar");
	if (!(prob->Tbar = (sparse_matrix *) mem_malloc (sizeof(sparse_matrix))))
		err_msg("Allocation", "main", "prob->Tbar");
	if (!(prob->Tbar->row = arr_alloc(p->matsz+1, int)))
		err_msg("Allocation", "sd", "Tbar->row");
	if (!(prob->Tbar->col = arr_alloc(p->matsz+1, int)))
		err_msg("Allocation", "sd", "Tbar->col");
	if (!(prob->Tbar->val = arr_alloc(p->matsz+1, double)))
		err_msg("Allocation", "sd", "Tbar->val");
	if (!(prob->Rbar->row = arr_alloc(p->marsz+1, int)))
		err_msg("Allocation", "sd", "Rbar->row");
	if (!(prob->Rbar->val = arr_alloc(p->marsz+1, double)))
		err_msg("Allocation", "sd", "Rbar->val");

	prob->Tbar->cnt = 0;
	prob->Rbar->cnt = 0;

	/*Procedure used to Extract Rbar*/
	for (r = row; r < p->mar; r++)
	{
		if (p->rhsx[r] > 0.000001)
		{
			prob->Rbar->val[prob->Rbar->cnt + 1] = p->rhsx[r];
			prob->Rbar->row[prob->Rbar->cnt + 1] = r - row + 1;
			++prob->Rbar->cnt;
		}
	}

	/*Procedure used to Extract Tbar*/
	for (c = 0; c < col; c++)
	{
		for (idx = p->matbeg[c]; idx < p->matbeg[c] + p->matcnt[c]; idx++)
		{
			if (p->matind[idx] >= row)
			{
				/* The coefficient is part of the subproblem's T matrix */
				prob->Tbar->val[prob->Tbar->cnt + 1] = p->matval[idx];
				prob->Tbar->row[prob->Tbar->cnt + 1] = p->matind[idx] - row + 1;
				prob->Tbar->col[prob->Tbar->cnt + 1] = c + 1;
				++prob->Tbar->cnt;
			}
		}
	}

	sd_global->Abar = 0.0;
	for (j = 1; j <= prob->Rbar->cnt; j++)
	{
		sd_global->Abar = sd_global->Abar
				+ y_i[row + prob->Rbar->row[j]] * prob->Rbar->val[j];
	}

	for (j = 0; j <= col + 2; j++)
	{
		sd_global->Bbar[j] = 0.0;
	}

	/*Calc Bbar*/
	for (j = 1; j <= prob->Tbar->cnt; ++j)
	{
		sd_global->Bbar[prob->Tbar->col[j]] =
				sd_global->Bbar[prob->Tbar->col[j]]
						+ y_i[row + prob->Tbar->row[j]] * prob->Tbar->val[j];
	}
	sd_global->Bbar[col + 2] = -1;

	/* Free the sparse vecotrs and matrices in the prob_type structure */
	mem_free(y_i);
	mem_free(prob->Tbar->row);
	mem_free(prob->Tbar->col);
	mem_free(prob->Tbar->val);
	mem_free(prob->Tbar);
	mem_free(prob->Rbar->row);
	mem_free(prob->Rbar->val);
	mem_free(prob->Rbar);
	mem_free(prob);
}

void generate_seed(sd_long * seed1, sd_long * seed2)
{
	int idx, cnt;
    sd_long rseed1, rseed2;

	printf("time:%ld\n", time(NULL) % 3600);
    srand((unsigned int) time(NULL));
	for (cnt = 0; cnt < BATCH_SIZE; cnt++)
	{
		for (idx = 0; idx < 4; idx++)
		{
			rseed1 = rand();
			/* printf("rseed1 before:%lx\n",rseed1); */
			rseed1 = rseed1 & 0x000000000000FFFF;
			/* printf("rseed1 after:%lx\n\n",rseed1); */
          if (idx == 0) {
            rseed2 = rseed1;
          }
          else{
            rseed2 = rseed2 << 16;
            rseed2 = rseed1|rseed2 ;
            /* printf("rseed2 after:%lx\n\n",rseed2);*/
          }
		}
		seed1[cnt] = rseed2;
	}
	printf("done\n");
}
