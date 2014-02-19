/***********************************************************************\
**
 ** theta.c
 **
 **
 ** f_k()
 ** get_max_Vi()
 ** append_thetas()
 ** new_theta()
 **
 ** History:
 **   18 Nov 1991 - <Jason Mai> - created.
 **   20 Jan 1992 - <Jason Mai> - prepared for testing / compiling.
 **      Feb 1992 - <Jason Mai> - changed structures.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include <float.h>
#include "theta.h"
#include "utility.h"
#include "log.h"
#include "sdglobal.h"

/***********************************************************************\
** This function makes an approximation of f(X) based on the
 ** current set of cuts from all ancestor cells.  f(X) is the 
 ** weighted sum of a bunch of max functions (one per ancestor cell).  
 ** It is calculated as:
 **
 **     f_k(x) = SUM { V_i * theta_i   |  cell i in current cell }
 ** 
 ** where each Vi is the result of:
 **
 **     Vi = MAX { alpha_k + (beta_k + c) x X  |  cut k in cell i } 
 **
 ** Here, alpha_k and beta_k are the coeffecients for the kth cut, 
 ** and must be multiplied by cut_outcomes / num_outcomes to update.
 **
 ** Each theta_i in the first expression is the cut coefficient for 
 ** cell i (representing of the cumulative effect of being married 
 ** with other cells) and is calculated as:
 **
 **     theta_i =  K_theta_i * P_theta_i / P_cell
 **
 ** The current cell's cuts must also be included in the approximation;
 ** however, they do not have any theta coefficients, so are computed
 ** separately at the end.  The index of the cut from the current cell 
 ** which has the greatest height is returned in the parameter _best_.
 \***********************************************************************/
double f_k(sdglobal_type* sd_global, vector x, vector cost, cut_type *cuts,
		theta_type *theta, num_type *num, cell_type *cur_cell, int *best)
{
	int cell, first;
	double sum, max, cx;

	/*
	 ** But theta is a part of cur_cell.... so are the cuts....
	 */

#ifdef LOOP
	printf("Inside f_k\n");
#endif

	/* Begin with the ancestor cuts, weighted by the theta coefficients */
	first = 0;
	sum = 0.0;
	for (cell = 0; cell < theta->cnt; cell++)
	{
		max = get_max_Vi(sd_global, x, cuts, num, &first, theta->last[cell],
				best, cur_cell);
		sum += max * theta->k[cell] * theta->p[cell] / cur_cell->P;
	}
	if (cur_cell->N)
		sum *= cur_cell->N / (double) (cur_cell->N + cur_cell->k);

	/* Now add the cuts originating in this cell, recognizing that the */
	/* cur_cell->k from the update and from the ancestral weighting cancel out */
	max = get_max_Vi(sd_global, x, cuts, num, &first, cuts->cnt, best,
			cur_cell);
	sum += max / (double) (cur_cell->N + cur_cell->k);

	/* Calculate the fixed cost for this solution to the master problem */
	cx = CxX(cost, x, num->mast_cols);

#ifdef CAL_CHECK
	printf("sum is :%f\n",sum);
	printf("cx is :%f\n",cx);
#endif

	return sum + cx;
}

/***********************************************************************\
** This function loops through a section of cuts specified by the
 ** cut and last_cut parameters, in order to find the one which has
 ** the greatest value at a given X vector.  It returns the
 ** value of the cut-function at this X, multiplied by the number
 ** of observations on which the cut was based.
 \***********************************************************************/
double get_max_Vi(sdglobal_type* sd_global, vector x, cut_type *cuts,
		num_type *num, int *this_cut, int last_cut, int *best, cell_type *c)
{
	double max, Vi;
	/*double	beta;*/
#ifdef LOOP
	printf("Inside get_max_Vi\n");
#endif

	/* The max should not be 0.0....  must be lower !!!!!! */
	/* But if there are no cuts, should have eta = 0 */

	/*
	 max = 0.0;
	 */
	/* Since 0.0 is not a lower bound for negative objective
	 we need negative infinity here Yifan 02/26/2012*/

	/* yifan debug
	 for (idx=0; idx<last_cut; idx++) {
	 printf("*******************************\n");
	 printf("cuts->val[%d]->subfeaflag is %d\n",idx,cuts->val[idx]->subfeaflag);
	 print_vect(cuts->val[idx]->beta, num->mast_cols, "Beta");
	 }*/

	/* Check all the cuts, from this_cut to last_cut, to find the max */
	for (max = -DBL_MAX; *this_cut < last_cut; (*this_cut)++)
	{
		/*added by Yifan to avoid evaluating feasibility cut 02/26/2012*/
		if (cuts->val[*this_cut]->subfeaflag == FALSE)
			continue;

		/* Calculate the part of the cut which depends on X */
		/* Yifan 03/14/2012 Updated for optimality cut height*/
		/*beta = CxX(cuts->val[*this_cut]->beta, x, num->mast_cols);*/

		/*
		 ** Compare the height of this cut with the maximum seen so far.
		 ** Recall that each cut is weighted by the number of samples used
		 ** in its formation.  We divide by the total number of samples later,
		 ** as part of the theta.k[] weight.
		 */
		/* Yifan 03/14/2012 Updated for optimality cut height*/
		/*Vi = (cuts->val[*this_cut]->cut_obs) * (cuts->val[*this_cut]->alpha - beta);*/

		/* Yifan 03/14/2012 Updated for optimality cut height*/
		Vi = c->k * cut_height(sd_global, cuts->val[*this_cut], x, c, num);

		if (Vi > max)
		{
			max = Vi;
			*best = *this_cut;
		}

#ifdef LOOP
		printf("Vi is %f for this_cut = %d\n", Vi, *this_cut);
#endif

	}

#ifdef DEBUG
	printf("Returning a max of %f\n", max);
#endif

	/* Return the value of the highest cut, times the number of cut_observ. */
	return max;
}

/***********************************************************************\
** This function calculates and returns the height of a given cut
 ** at a given X.  It includes the k/(k-1) update, but does not include
 ** the coefficients due to the cell.
 **
 ** Include cell info !
 \***********************************************************************/
double cut_height(sdglobal_type* sd_global, one_cut *cut, vector X,
		cell_type *c, num_type *num)
{
	double height;
	double t_over_k = ((double) cut->cut_obs / (double) c->k);

#ifdef LOOP
	printf("Inside cut_height\n");
#endif

	/* A cut is calculated as Alpha - Beta x X */
	height = cut->alpha;
	height -= CxX(cut->beta, X, num->mast_cols);

	/* Weight cut based on number of observations used to form it */
	height *= t_over_k;

	/* Yifan 03/14/2012 Updated for optimality cut height*/
	if (sd_global->config.LB_TYPE == 1)
	{
		height += (1 - t_over_k) * sd_global->Eta0;
	}
	else
	{
		height += 0.0;
	}

	return height;
}

double c_k(sdglobal_type* sd_global, vector Cost, vector X, cell_type *c,
		num_type *num, int idx)
{
	double ans;

	ans = cut_height(sd_global, c->cuts->val[idx], X, c, num);
	ans += CxX(Cost, X, num->mast_cols);

	return ans;
}

/***********************************************************************\
** This function will combine the theta arrays from two cells
 ** into a single theta array for one (married) cell.  It assumes
 ** that the cuts of the second cell are appended to the END of the
 ** cuts of the first cell.  It also assumes, of course, that the
 ** ordering of the cuts within each cell is not changed.  But,
 ** let me not to the marriage of two cells admit impediments...
 theta_type append_thetas(theta_type cell1, theta_type cell2, int n1, int n2)
 {
 theta_type	th;
 int		cnt, idx;

 th = new_theta(n1 + n2);
 
 for (cnt = 0; cnt < n1; cnt++)
 {

 }

 for (idx = 0; idx < n2; idx++)
 {

 }

 return th;
 }
 \***********************************************************************/

/***********************************************************************\
** This function dynamically allocates a new array of thetas
 ** (coefficients for the cuts of each ancestor cell in the
 ** current cell) based on the number of cells.  It returns
 ** a pointer to the new array.
 \***********************************************************************/
theta_type *new_theta(int num_cells)
{
	theta_type *theta;

#ifdef TRACE
	printf("Inside new_theta\n");
#endif

	if (!(theta = (theta_type *) mem_malloc (sizeof(theta_type))))
		err_msg("Allocation", "new_theta", "theta");

	if (!(theta->last = arr_alloc(num_cells, int)))
		err_msg("Allocation", "new_theta", "theta->last");

	if (!(theta->k = arr_alloc(num_cells, double)))
		err_msg("Allocation", "new_theta", "theta->k");

	if (!(theta->p = arr_alloc(num_cells, double)))
		err_msg("Allocation", "new_theta", "theta->p");

	theta->cnt = num_cells;

	return theta;
}

/***********************************************************************\
\***********************************************************************/
void free_theta(theta_type *theta)
{

#ifdef TRACE
	printf("Inside free_theta\n");
#endif

	mem_free(theta->last);
	mem_free(theta->k);
	mem_free(theta->p);
	mem_free(theta);
}
