/**************************************************************************\
**
 ** soln.h
 **
 **
 ** This header file contains the typedefs and prototypes used by 
 ** solve_cell() and all its subroutines for solving a single SD problem.
 **
 ** These structures are required for each serial SD solution, but they may 
 ** be freed once a solution is obtained (whereas the cell data must be 
 ** stored for future marriages, and the prob data is always needed).
 **
 **
 ** History:
 **   17 Feb 1992 - <Jason Mai> - created from subprob.h.
 **   16 Nov 1992 - <Jason Mai> - updated documentation.
 **
 \**************************************************************************/

#ifndef SOLN_H_
#define SOLN_H_

#include <time.h>
#include "cell.h"
#include "sdglobal.h"

/**************************************************************************\
** Omega stores the set of observations which have been made so far.
 **
 **   Each observation consists of a vector of realizations of random
 ** variables with discrete distributions.  Since every distribution is
 ** discrete, an observation is just stored as a vector of indices into a 
 ** distribution array containing the possible values.  _idx_ is an array
 ** of such vectors.  Each rv occurs in the R vector or T matrix which
 ** (along with the candidate X) make up the rhs of the subproblem.  
 **
 **   The _row_ and _col_ arrays give the coordinates in R and T of each rv 
 ** realization in a vector.  If _col_ is zero for a given entry, then the
 ** realization comes from R; otherwise, it comes from T.  The field _RT_ 
 ** represents a vector of actual realizations (as opposed to indices) of 
 ** omega for both the Romega and Tomega structures (only one observation's 
 ** worth of omega). 
 ** 
 **   The _weight_ field specifies the number of times a particular outcome 
 ** has been observed (it starts at 1 when the outcome is first generated, 
 ** and increments every time the same outcome is observed again).  _cnt_ 
 ** just specifies the number of distinct outcomes which have been observed 
 ** and stored in the omega structure.  
 **
 **   If the problem gets too large, some unexciting omegas may be dropped.  
 ** So, _filter_ (same length as _weight_) describes which vectors in the 
 ** _idx_ array are actually filled with observations.  _next_ references 
 ** the place to start in the _filter_ array when trying to find the next 
 ** available position in _idx_.  _most_ represents the number of
 ** elements needed to store all elements in the _idx_ array, from
 ** the first to the last.  (It is the greatest index at which an 
 ** observation is stored in the _idx_ array, plus one).   Finally,
 ** _last_ indicates the index of omega from which we started dropping
 ** omegas last time -- so we only "need" to go from omega.most down
 ** to omega.last when dropping them this time (we'll miss some).
 \**************************************************************************/
typedef struct
{
	int cnt;
	int next;
	int most;
	int last;
	int k; /* modified by Yifan 2012.07.02 iteration number*/
	int *row;
	int *col;
	sd_small **idx;
    sd_long *fidx;
	int *weight;
	int *filter;
	int **batch_idx;
	BOOL *used_opt_cut; /* Yifan 03/20/2012 Test for omega issues*/
	double *RT;
} omega_type;

/**************************************************************************\
**   The delta matrix contains the values of lambda_pi X Romega and
 ** lambda_pi X Tomega for all values of pi and all observations of omega.  
 ** _col_ gives the column number of the each non-zero element in the 
 ** multiplication of lambda_pi X Tomega (the same elements are non-zero 
 ** each time).  _val_ is an array of vectors of (lambda_pi X Romega, 
 ** lambda_pi X Tomega) pairs of calculations.  A row in _val_ corresponds
 ** to a distinct dual vector, and a column in _val_ corresponds to a 
 ** distinct observation of omega.  Thus, every pi-omega combination is
 ** represented here, and the size of the delta matrix can be determined
 ** from lambda->cnt and omega->cnt.
 **
 **   Note that when elements of omega get dropped, vacant columns appear 
 ** in delta.  This is ok, but be sure to loop carefully!
 \**************************************************************************/
typedef struct
{
	int *col;
	pi_R_T_type **val;
} delta_type;

/**************************************************************************\
**   When calculating istar for a cut, it is useful to have two separate
 ** references into the sigma and delta structures, since each dual vector
 ** is stored in two places -- part in sigma and part in delta.  The final
 ** entry in cut->istar[] will just be the _sigma_ field of this structure.
 \**************************************************************************/
typedef struct
{
	int delta;
	int sigma;
} i_type;

/**************************************************************************\
** The time_type struct is to record the time spent on different procedures
 ** of the SD algorithm: namely, the total time spent on SD algorithm; the
 ** time spent on one complete iteration; the solution time on solving master
 ** LPs (QPs in Regularized SD); solution time on solving subproblem LPs; 
 ** time spent the full optimality test; and time on the argmax procedures. 
 ** The time spent on the last five of the above procedures are recorded for
 ** both the current iteration and the accumutively of the algorithm. 
 ** added by zl, 06/29/04. 
 \**************************************************************************/
typedef struct
{
	double total_time;
	double iteration_time;
	double soln_master_iter;
	double soln_subprob_iter;
	double full_test_iter;
	double argmax_iter;
	double iteration_accum;
	double soln_master_accum;
	double soln_subprob_accum;
	double full_test_accum;
	double argmax_accum;
} time_type;

/**************************************************************************\
**   Whereas the prob and cell structures contain relatively permanent data,
 ** the soln structure contains temporary data necessary for solving a 
 ** single SD problem.  It keeps track of all realizations observed 
 ** (in _omega_) and calculations based on the realizations (_delta_).
 ** To avoid allocating and freeing vectors every iteration, the soln
 ** structure also contains a number of solution vectors (_Pi_, _candid_x_,
 ** and _incumb_x_) allocated once for the duration of the solution. 
 ** The _candid_est_, _incumb_est_, and _incumb_stdev_ fields store info
 ** about each of the solution vectors, for use in controlling the algorithm.
 ** _last_update_ stores the iteration of the last time the incumbent cut
 ** was updated.  _gamma_ stores the expected (in the previous iteration)
 ** improvement of the objective function, and is used in deciding whether
 ** or not to update the incumbent solution. 
 ** Master_pi added by JH 4/8/98 for pre_test_2.
 ** Master_dj added by JH 5/19/98 for pre_test_2. 
 ** norm_d_k_1 added by zl 6/10/02 for updating the scaling factor in QP.
 ** norm_d_k added by zl 6/10/02 for updating the scaling factor in QP. 
 \**************************************************************************/
typedef struct
{
	vector Pi; /* Sub poblem dual multipliers*/
	double subobj_est; /*added by Yifan 02/01/12*/
	vector Master_pi;
	vector Master_dj;
	vector *Batch_pi;
	vector *Batch_dj;
	vector Batch_pi_mean; /* modified by Yifan 2012.10.05 */
	vector Batch_dj_mean; /* modified by Yifan 2012.10.05 */
	vector Batch_pi_stdev; /* modified by Yifan 2012.10.05 */
	vector Batch_dj_stdev; /* modified by Yifan 2012.10.05 */
	vector R_Master_pi_mean; /* modified by Yifan 2012.10.05 */
	vector R_Master_dj_mean; /* modified by Yifan 2012.10.05 */
	vector R_Master_pi_stdev; /* modified by Yifan 2012.10.05 */
	vector R_Master_dj_stdev; /* modified by Yifan 2012.10.05 */
	/* modified by Yifan 2012.10.05 */
	double Obj_lb_mean;
	double Obj_lb_stdev;
	double Obj_lb_L;
	double Obj_lb_U;
	double Obj_comp_mean;
	double Obj_comp_stdev;
	double Obj_mean_mean;
	double Obj_mean_stdev;
	double rep_mean;
	double rep_stdev;
    double incumb_eta;  /* modified by Yifan 2013.9.16 */
	vector candid_x;
	double candid_est;
    vector mip_x;
	vector incumb_x;
	int incumb_k;
	vector incumb_d; /*Initially store d, later become x since we add incumbent to it directly*/
	vector incumb_avg;
	double c_xc;
	double alpha;
	double *beta;
	double *xc_height; /* Used to store cut height estimated by compromise solutions*/
	double incumb_est;
	double opt_value;
	double norm_d_k_1;
	double norm_d_k;
	double incumb_stdev;
	int incumb_cut;
	int last_update;
	double gamma;
	time_type *run_time; /* added by zl, 06/29/04. */
	BOOL optimality_flag;/* added by zl, 06/30/04. */
	BOOL smpl_ever_flag; /* added by zl, 08/20/04. */
	BOOL smpl_test_flag; /* added by zl, 08/17/04. */
	BOOL *dual_statble_flag; /*added by Yifan 09/27/2011*/
	double full_test_error;/* added by zl, 08/17/04. average error when failed in full test. */
	double *pi_ratio; /*added by Yifan 09/27/2011*/
	double max_ratio; /*added by Yifan 09/27/2011*/
	double min_ratio; /*added by Yifan 09/27/2011*/
	int passed; /* added by zl, 06/30/04. # of replications passed in full test. */
	double sub_lb_checker; /* addedy by zl, 07/01/04. check the valid lower bound of the subproblem obj function value. */
	omega_type *omega;
	delta_type *delta;
	BOOL incumbent_change;
	delta_type *feasible_delta;
} soln_type;

int *find_rows(int num_rv, int *num_rv_rows, int *omega_row);
int *find_cols(int num_rv, int *num_rv_cols, int *omega_col);
soln_type *new_soln(sdglobal_type* sd_global, prob_type *p, vector x_k);
void free_soln(prob_type *p, cell_type *c, soln_type *s);
void print_soln(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, char *fname);
int print_detailed_soln(sdglobal_type* sd_global, soln_type *s, prob_type *p,
		char *fname, int status);
int test_average_soln(sdglobal_type* sd_global, soln_type *s, prob_type *p);
void evaluate_inc(sdglobal_type* sd_global, cell_type *cell, prob_type *prob,
		soln_type *s, vector x_k, char *fname, double *conf_interval,
		int status);
void replace_incumbent(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, double *phi);
void resolve_infeasibility(sdglobal_type* sd_global, cell_type *c, prob_type *p,
		soln_type *s, int omeg_idx, BOOL new_omega, double *phi);
void calc_quad_scalar(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, double candid_est);
double calc_q_vec(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, double *q_vec);
void record_master_dual_and_obj(sdglobal_type* sd_global, prob_type *prob,
		soln_type *soln, FILE *master_dual, FILE *master_obj);
void process_batch_prob(sdglobal_type* sd_global, prob_type *prob,
		soln_type *soln, double *obj, FILE *batch_dual, FILE *batch_d,
		FILE *batch_obj);
void calc_lowerbound(sdglobal_type* sd_global, prob_type *prob, soln_type *soln);

#endif
