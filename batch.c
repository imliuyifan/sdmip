//
//  batch.c
//  sd1.4
//
//  Created by Yifan Liu on 9/9/12.
//
//

#include <time.h>
#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "log.h"
#include "sdglobal.h"

one_problem * new_batch_problem(one_problem * master, int max_cuts)
{
	one_problem *copy;
	int r, i, j, idx, cnt;
	sd_long col_offset, row_offset, batch_col_offset, batch_row_offset; //modified by Yifan to avoid memory leaks July 26 2011
	char batch_name[NAME_SIZE] =
	{ "   " };
	char *q;

#ifdef TRACE
	printf("Inside new_batch_problem\n");
#endif

	if (!(copy = (one_problem *) mem_malloc (sizeof(one_problem))))
		err_msg("Allocation", "new_batch_problem", "copy");

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

	/* Initialize dimensions of copy based on master. */
	copy->matsz = BATCH_SIZE * master->matsz
			+ BATCH_SIZE * (max_cuts) * (master->mac + 1);
	copy->marsz = BATCH_SIZE * master->mar + BATCH_SIZE * max_cuts;
	copy->mar = BATCH_SIZE * master->mar;
	copy->macsz = BATCH_SIZE * (master->mac + 1);
	copy->mac = BATCH_SIZE * (master->mac + 1);
	copy->cstorsz = BATCH_SIZE
			* (master->cstorsz + BATCH_SUFFIX * master->mac + NAME_SIZE);
	copy->rstorsz = BATCH_SIZE * (master->rstorsz + BATCH_SUFFIX * master->mar)
			+ NAME_SIZE * (BATCH_SIZE * max_cuts);
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
	j = 0;

	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
		batch_name[0] = 'B';
		batch_name[1] = '0' + idx / 10 % 10;
		batch_name[2] = '0' + idx / 1 % 10;

		for (q = master->cname[0]; q < master->cname[0] + master->cstorsz; q++)
		{
			if (*q == '\0')
			{
				copy->cstore[i++] = batch_name[0];
				copy->cstore[i++] = batch_name[1];
				copy->cstore[i++] = batch_name[2];
			}
			copy->cstore[i++] = *q;
		}

		for (q = master->rname[0]; q < master->rname[0] + master->rstorsz; q++)
		{
			if (*q == '\0')
			{
				copy->rstore[j++] = batch_name[0];
				copy->rstore[j++] = batch_name[1];
				copy->rstore[j++] = batch_name[2];
			}
			copy->rstore[j++] = *q;
		}

	}

	strcpy(copy->name, "batch_mean");
	strcpy(copy->objname, master->objname);

	/* Calculate difference in pointers for master/copy row and column names */
	col_offset = copy->cstore - master->cname[0];
	row_offset = copy->rstore - master->rname[0];
	batch_col_offset = master->cstorsz + BATCH_SUFFIX * master->mac;
	batch_row_offset = master->rstorsz + BATCH_SUFFIX * master->mar;

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (i = 0; i < BATCH_SIZE; i++)
	{
		for (j = 0; j < master->mac; j++)
		{
			copy->objx[i * master->mac + j] = (1.0 / BATCH_SIZE)
					* master->objx[j];
			copy->bdu[i * master->mac + j] = master->bdu[j];
			copy->bdl[i * master->mac + j] = master->bdl[j];
			copy->cname[i * master->mac + j] = master->cname[j]
					+ BATCH_SUFFIX * j + col_offset + i * batch_col_offset;
			copy->matbeg[i * master->mac + j] = cnt;
			copy->matcnt[i * master->mac + j] = master->matcnt[j];
			for (idx = master->matbeg[j];
					idx < master->matbeg[j] + master->matcnt[j]; idx++)
			{
				copy->matval[cnt] = master->matval[idx];
				copy->matind[cnt] = master->matind[idx] + i * master->mar;
				cnt++;
			}
			/* This line of code is causing a issue in printing the Batch-Mean */
			//cnt += BATCH_SIZE * max_cuts;
		}
	}

	/* Copy all information concerning rows of master */
	for (i = 0; i < BATCH_SIZE; i++)
	{
		for (r = 0; r < master->mar; r++)
		{
			copy->rhsx[i * master->mar + r] = master->rhsx[r];
			copy->senx[i * master->mar + r] = master->senx[r];
			copy->rname[i * master->mar + r] = master->rname[r]
					+ BATCH_SUFFIX * r + row_offset + i * batch_row_offset;
		}
	}

	/*
	 ** Initialize information for the extra columns in the new batch mean problem.
	 */

	for (i = 0; i < BATCH_SIZE; i++)
	{
		batch_name[0] = 'e';
		batch_name[1] = 't';
		batch_name[2] = 'a';
		batch_name[3] = 'B';
		batch_name[4] = '0' + i / 10 % 10;
		batch_name[5] = '0' + i / 1 % 10;

		/* Be careful! '7' means the size of batch_name here is 7. Yifan 2012-09-09*/
		strcpy(copy->cstore + BATCH_SIZE * batch_col_offset + i * 7,
				batch_name);
		copy->cname[BATCH_SIZE * master->mac + i] = copy->cstore
				+ BATCH_SIZE * batch_col_offset + i * 7;
		copy->objx[BATCH_SIZE * master->mac + i] = (1.0 / BATCH_SIZE);
		copy->bdu[BATCH_SIZE * master->mac + i] = INFBOUND;
		copy->bdl[BATCH_SIZE * master->mac + i] = -INFBOUND;
		copy->matbeg[BATCH_SIZE * master->mac + i] = cnt;
		copy->matcnt[BATCH_SIZE * master->mac + i] = 0;
		/* Be  careful! 0 means no cuts are added until now. Yifan 2012-09-09*/

	}

	/* Load the copy into CPLEX now */
	if (!(setup_problem(copy)))
		err_msg("Problem Setup", "new_master", "copy");

	/*write_prob(copy, "batch.lp");*//* added by Yifan to test batch.lp structure */

	/*
	 ** We're done, and we have room for max_cuts more constraints.
	 */

	return copy;
}

void add_cut_to_batch(sdglobal_type* sd_global, one_cut *cut, prob_type *p,
		cell_type *c, soln_type *s, int batch_id, int cut_position)
{
	int beg_col; /* column where beta coefficients begin */
	int end_col; /* column where beta coefficients finish */
	int *coef_col; /* column number of each beta coefficient */
	double *coef; /* used to store beta coefficient */
	int cnt, i;
	double rhs; /* rhs value in regularized QP method. */

#ifdef TRACE
	printf("Inside add_cut_to_batch\n");
#endif

	cut_position = p->num->mast_rows + cut_position;
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
		err_msg("Allocation", "add_cut_to_batch", "coef_col");
	if (!(coef = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "add_cut_to_batch", "coef");

	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
		coef_col[cnt] = cnt + p->num->mast_cols * batch_id;
	coef_col[p->num->mast_cols] = p->num->mast_cols * BATCH_SIZE + batch_id;

#ifdef DEBUG
	printf("Adding the row:\n");
	print_vect(cut->beta, end_col - beg_col, "c->beta");
#endif

	/* 1. Get the coefficient beta */
	for (i = 0; i <= p->num->mast_cols; i++)
	{
		get_coef(c->master, cut_position, i, &coef[i]);
		/* modified by Yifan 2013.02.15 */
		//cut->beta[i+1] = coef[i];
	}

	/* 2. Get the coefficient in front of eta i=p->num->mast_cols */

	/* 3. Get the RHS of the cut */
	get_rhs(c->master, &rhs, cut_position, cut_position);
    
    /* modified by Yifan 2013.9.18 */
    rhs = cut->alpha;
    
	/* modified by Yifan 2013.02.15 */
	//cut->alpha_incumb = rhs;
	if (!add_row_to_batch(sd_global->batch_problem, beg_col,
			end_col - beg_col + 1, coef_col, coef, GE, rhs, batch_id))
		err_msg("LP solver", "add_cut_to_bach", "ans");
	/* THIS IS A SERIOUS PROBLEM, ROW_NUM IS NOT CORRECT, ONLY ONE CUT IN MASTER!!!*/
	/* Yifan 03/11/2012 Be careful of this row number*/

	mem_free(coef_col);
	mem_free(coef);

	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl
	 
	 print_cut_info(c, p->num, "After adding a cut");
	 */
#ifdef TRACE
	printf("Exiting add_cut_to_batch\n");
#endif

}

void add_fcut_to_batch(sdglobal_type* sd_global, one_cut *cut, prob_type *p,
		cell_type *c, soln_type *s, int batch_id, int iteration_num)
{
	int beg_col; /* column where beta coefficients begin */
	int end_col; /* column where beta coefficients finish */
	int *coef_col; /* column number of each beta coefficient */
	double *coef; /* used to store beta coefficient */
	int cnt, i;
	double rhs; /* rhs value in regularized QP method. */

#ifdef TRACE
	printf("Inside add_cut_to_batch\n");
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
		err_msg("Allocation", "add_cut_to_batch", "coef_col");
	if (!(coef = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "add_cut_to_batch", "coef");

	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
		coef_col[cnt] = cnt + p->num->mast_cols * batch_id;
	coef_col[p->num->mast_cols] = p->num->mast_cols * BATCH_SIZE + batch_id;

#ifdef DEBUG
	printf("Adding the row:\n");
	print_vect(cut->beta, end_col - beg_col, "c->beta");
#endif

	/* 1. Get the coefficient beta */
	for (i = 0; i < p->num->mast_cols; i++)
	{
		coef[i] = cut->beta[i + 1];
	}

	/* 2. Get the coefficient in front of eta i=p->num->mast_cols This value is 0 for fcut */

	/* 3. Get the RHS of the cut */
	rhs = cut->alpha
			- CxX(cut->beta, sd_global->batch_incumb->incumb_x[batch_id],
					p->num->mast_cols);
	rhs += sd_global->config.FEA_TOLER;

	/* modified by Yifan 2013.02.15 */

	if (!add_row_to_batch(sd_global->batch_problem, beg_col,
			end_col - beg_col + 1, coef_col, coef, GE, rhs, batch_id))
		err_msg("LP solver", "add_cut_to_bach", "ans");
	/* THIS IS A SERIOUS PROBLEM, ROW_NUM IS NOT CORRECT, ONLY ONE CUT IN MASTER!!!*/
	/* Yifan 03/11/2012 Be careful of this row number*/

	mem_free(coef_col);
	mem_free(coef);

	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl
	 
	 print_cut_info(c, p->num, "After adding a cut");
	 */
#ifdef TRACE
	printf("Exiting add_cut_to_batch\n");
#endif

}

void update_batch_rhs(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int batch_id)
{

	int cnt;
	double *rhs;
	int *indices;

#ifdef TRACE
	printf("Inside update_batch_rhs\n");
#endif

	if (!(rhs = arr_alloc(p->num->mast_rows, double)))
		err_msg("Allocation", "change_rhs", "rhs");
	if (!(indices = arr_alloc(p->num->mast_rows, int)))
		err_msg("Allocation", "change_rhs", "indices");

	for (cnt = 0; cnt < p->num->mast_rows; cnt++)
	{
		get_rhs(c->master, &rhs[cnt], cnt, cnt);
		indices[cnt] = p->num->mast_rows * batch_id + cnt;
	}

	/* Now we change the rhs of the master problem. */
	change_rhside(sd_global->batch_problem, p->num->mast_rows, indices, rhs);
	mem_free(rhs);
	mem_free(indices);

#ifdef TRACE
	printf("Exiting update_batch_rhs\n");
#endif

}

void update_batch_bounds(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int batch_id)
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
	printf("Inside update_batch_bounds\n");
#endif
	if (!(lbounds = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "change_bounds", "lbounds");
	if (!(lindices = arr_alloc(p->num->mast_cols+1, int)))
		err_msg("Allocation", "change_bounds", "lindices");
	if (!(llu = arr_alloc(p->num->mast_cols+1, char)))
		err_msg("Allocation", "change_bounds", "llu");

	if (!(ubounds = arr_alloc(p->num->mast_cols+1, double)))
		err_msg("Allocation", "change_bounds", "ubounds");
	if (!(uindices = arr_alloc(p->num->mast_cols+1, int)))
		err_msg("Allocation", "change_bounds", "uindices");
	if (!(ulu = arr_alloc(p->num->mast_cols+1, char)))
		err_msg("Allocation", "change_bounds", "ulu");

	get_lbound(c->master, lbounds, 0, p->num->mast_cols);
	get_ubound(c->master, ubounds, 0, p->num->mast_cols);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
	{
		uindices[cnt] = cnt + batch_id * p->num->mast_cols;
		ulu[cnt] = 'U';
	}

	status = change_bound(sd_global->batch_problem, p->num->mast_cols, uindices,
			ulu, ubounds);

	/* Change the Lower Bound */
	for (cnt = 0; cnt < p->num->mast_cols; cnt++)
	{
		lindices[cnt] = cnt + batch_id * p->num->mast_cols;
		llu[cnt] = 'L';
	}

	status = change_bound(sd_global->batch_problem, p->num->mast_cols, lindices,
			llu, lbounds);

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
	printf("Exiting update_batch_bounds\n");
#endif
}

batch_incumb_type * new_batch_incumb(sdglobal_type* sd_global, prob_type *p,
		vector x_k)
{

	/* Remeber to release memory allocated here!!! */
	int idx;
	if (!(sd_global->batch_incumb =
			(batch_incumb_type *) mem_malloc (sizeof(batch_incumb_type))))
		err_msg("Allocation", "new_batch_incumb", "batch_incumb");

	if (!(sd_global->batch_incumb->incumb_x =
			(vector *) mem_calloc (BATCH_SIZE, sizeof(vector))))
		err_msg("Allocation", "incumb_x", "batch_incumb");

	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
		sd_global->batch_incumb->incumb_x[idx] = duplic_arr(x_k,
				p->num->mast_cols);
	}

	/* modified by Yifan 2012.10.05 */
	if (!(sd_global->batch_incumb->R_Master_pi =
			(vector *) mem_calloc (BATCH_SIZE, sizeof(vector))))
		err_msg("Allocation", "new_soln", "Batch_pi");
	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
		sd_global->batch_incumb->R_Master_pi[idx] =
				arr_alloc(p->num->mast_rows+p->num->max_cuts+1,double);
	}
	if (!(sd_global->batch_incumb->R_Master_dj =
			(vector *) mem_calloc (BATCH_SIZE, sizeof(vector))))
		err_msg("Allocation", "new_soln", "Batch_pi");
	for (idx = 0; idx < BATCH_SIZE; idx++)
	{
		sd_global->batch_incumb->R_Master_dj[idx] =
				arr_alloc(p->num->mast_cols+2,double);
	}

	return sd_global->batch_incumb;
}

void save_batch_incumb(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int batch_id)
{
	FILE *bat;
	int idx;
	for (idx = 0; idx <= p->num->mast_cols; idx++)
	{
		sd_global->batch_incumb->incumb_x[batch_id][idx] = s->incumb_x[idx];
	}
	bat = fopen("incumb.out", "a");

	for (idx = 1; idx <= p->num->mast_cols; idx++)
	{
		fprintf(bat, "%f\n", sd_global->batch_incumb->incumb_x[batch_id][idx]);
	}
	fprintf(bat, "\n");
	fclose(bat);

}

void add_batch_equality(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s)
{
	int beg_col; /* column where beta coefficients begin */
	int nzcnt;
	int *coef_col; /* column number of each beta coefficient */
	double *coef; /* used to store beta coefficient */
	int i, j;
	double rhs; /* rhs value in regularized QP method. */

#ifdef TRACE
	printf("Inside add_batch_equality\n");
#endif

	/*
	 ** Initialize an array to specify columns of each coefficient in
	 ** beta.  The one-norm of beta is temporarily used as the coefficeint on
	 ** eta (it is assumed to be replaced in the next step in solve_master()).
	 */
	beg_col = 0;
	nzcnt = 2;
	/*
	 ** Assign the column values for each element in the beta vector.
	 ** The zeroth position will be set to the eta column, which will
	 ** temporarily receive a coefficient equal to the 1-norm of beta.
	 ** (or whatever is stored in beta[0])
	 **
	 ** This ought to be reworked... shouldn't have to alloc every time...
	 */
	if (!(coef_col = arr_alloc(2, int)))
		err_msg("Allocation", "add_cut_to_batch", "coef_col");
	if (!(coef = arr_alloc(2, double)))
		err_msg("Allocation", "add_cut_to_batch", "coef");

	for (i = 0; i < BATCH_SIZE - 1; i++)
	{
		for (j = 0; j < p->num->mast_cols; j++)
		{
			coef_col[0] = i * p->num->mast_cols + j;
			coef_col[1] = i * p->num->mast_cols + j + p->num->mast_cols;
			coef[0] = 1.0;
			coef[1] = -1.0;
			rhs = sd_global->batch_incumb->incumb_x[i + 1][j + 1]
					- sd_global->batch_incumb->incumb_x[i][j + 1];
            /* modified by Yifan 2013.9.18 */
            rhs = 0;
			if (!add_row_to_batch(sd_global->batch_problem, beg_col, nzcnt,
					coef_col, coef, 'E', rhs, i))
				err_msg("LP solver", "add_cut_to_bach", "ans");
		}
	}

	mem_free(coef_col);
	mem_free(coef);
	/* Print all the cuts after adding a cut, for the purpose of cut index
	 checking. zl
	 
	 print_cut_info(c, p->num, "After adding a cut");
	 */
#ifdef TRACE
	printf("Exiting add_batch_equality\n");
#endif
}

BOOL get_beta_x(sdglobal_type* sd_global, soln_type *s, vector Beta,
		one_problem *p, num_type *num, int length)
{
#ifdef TRACE
	printf("Inside get_beta\n");
#endif

	int i, j, k;
	int row, col, mast_rows;
	double alpha, dual;

	mast_rows = BATCH_SIZE * num->mast_rows;

	for (i = 0; i < BATCH_SIZE; i++)
	{

		for (j = 0; j < num->max_cuts; j++)
		{
			dual = s->Batch_pi[i][num->mast_rows + j + 1] * BATCH_SIZE;
			if (dual > sd_global->config.TOLERANCE)
			{
				row = mast_rows + i * num->max_cuts + j;
				for (k = 0; k < num->mast_cols; k++)
				{
					col = i * num->mast_cols + k;
					get_coef(p, row, col, &Beta[k + 1]);
				}
				get_rhs(p, &alpha, row, row);
				/* modified by Yifan 2012.09.28 */
				/* This height is used to estimate the original problem's lower bound */
				s->xc_height[i] += (alpha
						+ CxX(Beta, sd_global->batch_incumb->incumb_x[i],
								num->mast_cols)) * dual;
			}
		}

	}

#ifdef TRACE
	printf("Exiting get_beta\n");
#endif
	return TRUE;
}

