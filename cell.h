/**************************************************************************\
**
 ** cell.h
 **
 **
 ** This header file contains the prototypes and type definitions
 ** required for the solution of a problem with one or more cells.
 **
 ** Its data structures are necessary to define one cell, which contains
 ** dual vectors, cuts, and some counts for combining them. 
 **
 **
 **
 ** History:
 **   17 Feb 1992 - <Jason Mai> - created from subprob.h.
 **   16 Nov 1992 - <Jason Mai> - updated documentation.
 **
 **
 \**************************************************************************/

/**************************************************************************\
**   To save time and space, Pi x R and Pi x T are calculated as soon
 ** as possible and stored in structures like sigma and delta.  Toward
 ** this end, pi_R_T_type represents a single calculation of pi X R 
 ** (which is a scalar) and pi X T (which is a vector).
 \**************************************************************************/
#ifndef CELL_H_
#define CELL_H_

#include "sdglobal.h"

typedef struct
{
	double R;
	double *T;
} pi_R_T_type;

/**************************************************************************\
**   The lambda structure stores some of the dual variable values from every 
 ** distinct dual vector obtained during the program.  Each vector contains 
 ** only those dual variables whose corresponding rows in the subproblem
 ** constraint matrix contain random elements.  _val_ is an array of 
 ** these dual vectors (thus it is 2-D).  _row_ gives the corresponding 
 ** row number for a given dual variable in _val_.  _cnt_ represents 
 ** the number of dual vectors currently stored in lambda.
 \**************************************************************************/
typedef struct
{
	int cnt;
	int *row;
	double **val;
} lambda_type;

/**************************************************************************\
**   The sigma matrix contains the values of Pi x Rbar and Pi x Tbar
 ** for all values of pi obtained so far (note it does not depend on
 ** observations of omega).  _col_ gives the column number of 
 ** each non-zero element in pi X Tbar.  _val_ is an array of values
 ** for pi X Rbar and pi X Tbar, one entry for each pi.  Note that values
 ** which are always zero (because Rbar or Tbar is zero there) are not 
 ** stored.  The _lamb_ array is the same size as the _val_ array, and for 
 ** each element in _val_ the corresponding element in _lamb_ references 
 ** the dual vector in lambda that was used to calculate that entry in sigma. 
 \**************************************************************************/
typedef struct
{
	int cnt;
	int *col;
	pi_R_T_type *val;
	int *lamb;
	int *ck; //record the iteration # of the sigma
} sigma_type;




/**************************************************************************\
**   The current cell's long array of cuts is composed of the cuts of all 
 ** the ancestors' cells, as well as its own cuts.  The theta structure 
 ** describes all the attributes of an ancestor's cuts.  Each ancestor
 ** is represented by one element in the _last_, _k_, and _p_ arrays.
 ** The _last_ array specifies the position in this cell's cuts->val array 
 ** of the last cut belonging to an ancestor.  The _k_ array represents the 
 ** coefficient on all cuts of an ancestor due to the cumulative multiplication 
 ** by k/N when cells are married.  The _p_ field represents the original 
 ** probability of the ancestor.  The _cnt_ field simply gives the number of 
 ** ancestors (equivalently, the length of each array in theta).
 **
 **   Should theta be one long array, as before, to ease the marriage rountine?
 ** harder to reference; easier to merge...
 **
 \**************************************************************************/
typedef struct
{
	int cnt;
	int *last;
	double *k;
	double *p;
} theta_type;

/**************************************************************************\
**   In order to marry two cells together, the dual vectors and cuts
 ** from each cell must be combined, using weighting factors based
 ** on the number of samples taken in each cell and the probability 
 ** associated with each cell.  So, each cell maintains the dual
 ** vector data (_lambda_ and _sigma_), the cuts (_cuts_), and counts
 ** for weighting the cuts (_N_, _theta_, _P_).  Each cell
 ** also knows its id number as _id_num_.  _k_ represents the number 
 ** of samples taken thus far in the solution of the current cell.  
 ** 
 **   In order to generate observations from a married cell, we must
 ** know the first generation ancestors (those which once existed as 
 ** single cells) which have combined to make this cell.  The id numbers
 ** of these ancestors are stored in _members_, which is an array
 ** of length _num_members_.  When cells are running in parallel, each
 ** one should have its own copy of the master and subproblem data, so 
 ** there are no collisions or locks.  Pointers to these problems are stored
 ** in the _master_ and _subprob_ fields.
 **
 **   Should the theta structure be updated at the end of a cell's
 ** solution, or at the beginning of the married cell's solution?
 ** i.e. should num_samples (k) be stored in the cell, or later in theta?
 \**************************************************************************/
typedef struct
{
	int id_num;
	int num_members;
	int *members;
	one_problem *master;
	one_problem *subprob;
	/* Yifan 06/25/2012 batch mean */
	cut_type *cuts;
	lambda_type *lambda;
	sigma_type *sigma;
	theta_type *theta;
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	cut_type *feasible_cuts_pool;
	cut_type *feasible_cuts_added;
	lambda_type *feasible_lambda;
	sigma_type *feasible_sigma;
	theta_type *feasible_theta;
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	double quad_scalar; /* the quadratic scalar, 'sigma'.  zl */
	int LP_cnt; /* the LPs solved. zl. 06/30/02 */
	int LP_test; /* LPs solved for testing. zl. 07/02/02 */
	int N;
	double P;
	int k;
	BOOL opt_mode; /* Yifan 03/23/2012 used in deciding the mode for solve_QP master */
	BOOL incumb_infea; /* Yifan 03/23/2012 if incumbent infeasibility is encountered in form_fea_cut */
	int fea_count; /* Yifan 03/23/2012 count how many iterations the SD goes into the feasibility mode*/
} cell_type;

typedef struct
{
    int u_idx;
    double u;
} u_type; /* modified by Yifan 2013.9.25 */

typedef struct
{
    u_type uvalue[total_binary];
}u_group_type; /* modified by Yifan 2013.9.25 */

/***********************************************************************\
 **	 Batch type vs cell_type is similar to cut_type vs one_cut.
 **  20 cells can be saved in the batch_type and _cnt_
 **  will record that number.
 \***********************************************************************/
typedef struct
{
	int cnt;
	cell_type **val;
} batch_type;

/* This structure is used in reduce_cut for sorting cuts*/
typedef struct
{
	int obs;
	int cut_index;
} sort_type;

void solve_cell(sdglobal_type* sd_global, cell_type *cell, prob_type *prob,
		vector x_k, char *fname);
cell_type *new_cell(sdglobal_type* sd_global, prob_type *p, int id_num);
void free_cell(cell_type *c, num_type *num);
void write_cell(prob_type *p, cell_type *c, char *fname);

#endif

