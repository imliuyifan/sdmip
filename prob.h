/***********************************************************************\
**
 ** prob.h
 **
 **  In the first place, this file contains all the definitions necessary 
 **  for the routines in solver.c to solve LP's using a given optimizer.
 **  (Currently, the CPLEX optimizer is being used.  If the LP solver
 **   is changed, all the definitions in this file, and all the
 **   functions in solver.c, need to be re-written)
 **
 **  Second, this file contains fundamental typedefs and definitions 
 **  required for every source file in the SD executable.  It must be 
 **  included in all of them.
 **
 ** History:
 **   ?? ??? 91 - <Jason Mai> - created.
 **
 \**********************************************************************/

#ifndef PROB_H_
#define PROB_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sdconstants.h"
#include "sdglobal.h"
#include <string.h>

/*
 ** This definition ought to be in defs.inc, which all modules 
 ** must include in order to use CPLEX, but it isn't.  
 */
typedef void cpxlp;

/* The SOLVER_TYPE is to define the choice of the two sd method: basic LP
 or Regularized QP.  zl */

/**********************************************************************\
** The coord structure provides the row / column coordinates of several
 ** sparse vectors and matrices used in the SD algorithm.  Each field is
 ** an array of integers, such that array[i] gives the coordinate of the
 ** ith element of another array.  _omega_row_ and _omega_col_ contain
 ** the row and column number, respectively, of each random variable in 
 ** the problem.  _delta_col_ contains a list of all the columns in the
 ** T matrix which contain random variables -- and thus have non-zero
 ** values for PixTomega.  _lambda_row_, on the other hand, contains 
 ** a list of all rows of the subproblem which contain random variables. 
 \**********************************************************************/
typedef struct
{
	int *omega_row;
	int *omega_col;
	int *sigma_col;
	int *delta_col;
	int *lambda_row;
} coord_type;

/**********************************************************************\
** The num structure specifies the dimensions of the SD problem.
 ** _mast_rows_ and _mast_cols_ give the number of rows and columns,
 ** respectively, in the master problem constraint matrix (not including
 ** the eta column nor the cuts, which we will add later).  Likewise,
 ** _sub_rows_ and _sub_cols_ give the number of rows and columns in
 ** the subproblem constraint matrix.  _rv_rows_ and _rv_cols_ give the
 ** number of rows and columns in the T matrix which contain random
 ** elements (if a row has two rv's, it is only counted once.  Also,
 ** random elements of R aren't counted in _rv_cols_, but are in _rv_rows_).
 ** _nz_cols_ gives the number of columns in T which have at least one non-zero 
 ** element.  The _iter_ field specifies the maximum number of iterations allowed
 ** for solving the problem -- this determines how large some structures will
 ** be.  Similarly, the _max_cuts_ field specifies the maximum number of 
 ** cuts which may be in the master problem at any given time.  The _cipher_ 
 ** field specifies the number of integers required to encode a single 
 ** observation of omega. _rv_ simply gives the total number of random 
 ** variables in the problem.  The _rv_R and _rv_T_ fields specify the number 
 ** of random variables in the R vector and the T matrix, respectively.
 \**********************************************************************/
typedef struct
{
	int mast_rows;
	int mast_cols;
	int sub_rows;
	int sub_cols;
	int rv_rows;
	int rv_cols;
	int nz_cols;
	int max_cuts;
	int cipher;
	int iter;
	int rv;
	int rv_R;
	int rv_T;
	int batch_id;
} num_type;

/**********************************************************************\
** In a sparse vector, the non-zero elements (only) are stored 
 ** consecutively in the _val_ field.  Each non-zero element also
 ** has its corresponding row number stored in the _row_ field.
 ** Finally, _cnt_ provides the number of non-zero elements stored
 ** in the vector.
 \**********************************************************************/
typedef struct
{
	int cnt;
	int *row;
	double *val;
} sparse_vect;

/**********************************************************************\
** In a sparse matrix, the non-zero elements (only) are stored
 ** in consecutive positions in the _val_ array.  Each one of these
 ** non-zero elements has its row and column coordinates in the matrix
 ** stored in the _row_ and _col_ fields, respectively.  The _cnt_ field
 ** gives the total number of non-zero elements present in the matrix.
 \**********************************************************************/
typedef struct
{
	int cnt;
	int *row;
	int *col;
	double *val;
} sparse_matrix;

/**********************************************************************\
** A string is, obviously, an array of characters, hopefully null-terminated.
 \**********************************************************************/
typedef char *string;

/**********************************************************************\
** The deterministic part of the problem is stored in a structure of
 ** type prob_type.  This includes all the data for the master problem,
 ** all the data for the subproblem, the matrix Tbar and the vector Rbar
 ** used to calculate the right hand side of the subproblem, and the
 ** master problem cost coefficients, c, used to calculate f(x).
 ** It also contains _num_, which specifies various dimensions of the 
 ** problem (like the number of columns containing random variables, etc.),
 ** and _tau_, which specifies how often the incumbent cut is re-evaluated.
 ** The _coord_ field provides the row / column coordinates for important
 ** structures of SD which only need to be calculated once per problem.
 \**********************************************************************/
typedef struct
{
	one_problem *master;
	one_problem *subprob;
	sparse_vect *Rbar;
	sparse_matrix *Tbar;
	sparse_matrix *A; /* to store the A matrix in regularized QP method. zl */
	coord_type *coord;
	num_type *num;
	vector c;
	int tau;
	int current_batch_id; /* to record the current batch's number */
	sd_long eval_seed;
} prob_type;

prob_type *new_prob(sdglobal_type* sd_global, one_problem *original, int num_rv,
		int num_cipher, int row, int col);
void free_prob(sdglobal_type* sd_global, prob_type *p);
void init_param(sdglobal_type* sd_global, prob_type *p);
void solve_SD(sdglobal_type* sd_global, one_problem *original, vector x_k,
		int num_rv, int num_cipher, int row, int col, char *fname, int batch_id);
int decompose(sdglobal_type* sd_global, one_problem *orig, prob_type *p,
		int row, int col);
void free_one_prob(one_problem *p);
void generate_seed(sd_long * seed1, sd_long * seed2);
void err_msg(char *type, char *place, char *item);
void parse_cmd_line(sdglobal_type* sd_global, int argc, char *argv[],
		char *fname, int *objsen, int *num_probs, int *start, BOOL *read_seeds,
		BOOL *read_iters);
void get_lower_bound(sdglobal_type* sd_global, one_problem *SDptr, int row,
		int col);
void calc_alpha_beta(sdglobal_type* sd_global, one_problem *SDprobptr, int row,
		int col);

#endif
