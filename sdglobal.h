/*
 * sdglobal.h
 *
 *  Created on: Jun 12, 2013
 *      Author: lian
 */

#ifndef SDGLOBAL_H_
#define SDGLOBAL_H_
#include "sdconstants.h"
#include "time.h"
#include <stdlib.h>

typedef enum
{
	FALSE, TRUE
} BOOL;

/***** MAKE THIS unsigned int SO THAT BITWISE OPERATIONS WORK RIGHT ******/
typedef int sd_small; /* Type used for elements of omega */

typedef struct
{
	sd_small upper;
	sd_small lower;
} range;

/***************************************************************************\
**  Struct cell_map contains all stochastic information relative to
 ** an individual cell.  This includes the distribution according to which
 ** the cell generates observations of omega (_map_)
 \***************************************************************************/

typedef struct
{
	range *map; /* prob. range structures ordered by omega */
	sd_small cell_ID; /* cell number */
	float prob; /* conditional prob. of occurence of random element in cell*/
	float fam_prob; /* conditional probability within marriage */
} cell_map;

/**********************************************************************\
** The one_key structure contains the information required to bitwise
 ** encode or decode an integer into a segment of another integer.
 ** The source integer must lie within a known interval, and if an array
 ** of integers are to be encoded (this is the expected case) then
 ** each one should have its own key (so there'll be an array of keys,
 ** corresponding to the source array of integers).  Thus, _element_
 ** specifies which item in the ciphered array contains the given source
 ** integer.  _shift_ describes how many bits to the left the given source
 ** integer has been shifted in that cipher element.  And _mask_
 ** provides a bit pattern ready for AND-ing to strip off other coded
 ** integers from the bits corresponding to the desired integer.
 \**********************************************************************/
typedef struct
{
	int element;
	int shift;
	int mask;
} one_key;

/**********************************************************************\
**  struct identity is used to maintain a record of the order of
 ** columns and rows within the constraint matrix.  the names of the
 ** columns and rows are stored in the arrays col_name and row_name
 ** by index number.
 **
 **  This structure will be loaded by load_core during the reading of
 ** the core file.  This structure will be used by load omega for
 ** addressing the stochastic elements of the constraint matrix.  The
 ** structure will also be used for the decomposition of the constraint
 ** matrix following the running of the main problem in cplex.
 \**********************************************************************/

typedef struct
{
	char **row_name; /* array of row names ordered by index */
	char **col_name; /* array of column names ordered by index */
	sd_small num_rows, num_cols;
} identity;

/**********************************************************************\
** The config structure contains constants used by SD during the
 ** solution process.  They are read from a special configuration file,
 ** "config.sd", but they revert to default values if this file
 ** doesn't specify them.  See load_config() in input.c for info.
 \**********************************************************************/
typedef struct
{
	double CONFID_HI;
	double CONFID_LO;
	double PERCENT_PASS;
	double EPSILON;
	double PRE_EPSILON; /* for pretests. zl */
	double SMOOTH_PARM;
	int SMOOTH_I;
	double TOLERANCE; /* for zero identity test. zl */
	double FEA_TOLER; /* modified by Yifan 2013.05.05 */
	double THIN_TOLER;
	double R;
	double R2; /* for updating the scaling_factor. zl */
	double R3; /* for updating the scaling_factor. zl */
	double MIN_QUAD_SCALAR; /* Yifan 03/27 */
	double MAX_QUAD_SCALAR; /* added by zl 06/20/02.  */
	double SUBPROB_LB; /* the lb of the subprob obj values. zl, 07/01/04. */
	int SUB_LB_CHECK; /* for subproblem LB check. zl, 09/20/05. */
	double ITER_FACT;
	int EVAL_FLAG;
	double EVAL_ERROR;
	double MEAN_DEV;
	int CUT_MULT;
	int M;
	int TAU;
	int MIN_ITER;
	int MAX_ITER;
    int OVERRIDE;
	int START_THIN;
	int THIN_CYCLE;
	int PRINT_CYCLE;
	int EVAL_RUN_FLAG;
	int DROP_TIME;
	sd_long RUN_SEED;
	sd_long RUN_SEED1;
	sd_long RUN_SEED2;
	sd_long RUN_SEED3;
	sd_long RUN_SEED4;
	sd_long RUN_SEED5;
	sd_long RUN_SEED6;
	sd_long RUN_SEED7;
	sd_long RUN_SEED8;
	sd_long RUN_SEED9;
	sd_long RUN_SEED10;
	sd_long RUN_SEED11;
	sd_long RUN_SEED12;
	sd_long RUN_SEED13;
	sd_long RUN_SEED14;
	sd_long RUN_SEED15;
	sd_long RUN_SEED16;
	sd_long RUN_SEED17;
	sd_long RUN_SEED18;
	sd_long RUN_SEED19;
	sd_long RUN_SEED20;
	sd_long RUN_SEED21;
	sd_long RUN_SEED22;
	sd_long RUN_SEED23;
	sd_long RUN_SEED24;
	sd_long RUN_SEED25;
	sd_long RUN_SEED26;
	sd_long RUN_SEED27;
	sd_long RUN_SEED28;
	sd_long RUN_SEED29;
	sd_long RUN_SEED30;
	sd_long EVAL_SEED1;
	int MASTER_TYPE;
	int LB_TYPE;
	int TEST_TYPE;
	int PI_EVAL_START;
	int PI_CYCLE;
	int SCAN_LEN;
	int DETAILED_SOLN;
	int MULTIPLE_REP;
	int AUTO_SEED;
    double BOX_SIZE;
} config_type;

/**************************************************************************\
**  Struct omegastuff contains all pertinent information related to the
 ** stochastic elements in the SD problem.
 \**************************************************************************/
typedef struct
{
	sd_small num_omega; /* number of stochastic elements stored in structure */
	sd_small num_cells; /* number of cells being used in program */
	sd_small num_cipher; /* number of ints needed to encode an observation */
	sd_small *indices; /* array temporarily containing one observation's indices */
	sd_small *row; /* row number array */
	sd_small *col; /* column number array */
	sd_small *num_vals;
	double **omega_vals; /* indexed array of discrete elements in omega */
	double **omega_probs;/* indexed array of probabilities associated with discrete elements */
    double *mean;         /*added by Yifan to record the mean of each rv */
    char *file_name;
	cell_map *mapping; /* array of structures mapping cell ranges to omega. */
	one_key *key; /* array of keys for encoding / decoding obsevations */
} omegastuff;

/**********************************************************************\
** struct one_problem is used to store the LP problem information in
 ** a format usable by cplex.  The information can be loaded through the
 ** use of the function loadprob() which returns a pointer *lp to the
 ** data it receives.
 **
 ** struct cpxlp *
 ** loadprob(probname, mac, mar, mae, objsen, objx, rhsx, senx, matbeg,
 **	    matcnt, matind, matval, bdl, bdu, rngcol, nrowind, etype,
 **	    enzbeg, enzcnt, enzind, enzval, dataname, objname, rhsname,
 **	    rngname, bndname, cname, cstore, rname, rstore, ename,
 **	    estore, macsz, marsz, matsz, maesz, enzsz, cstorsz, rstorsz,
 **	    estorsz);
 **
 \**********************************************************************/
typedef struct
{
	char *name; /* problem name (up to 15 characters) */
	int mac; /* number of columns in the constraint matrix */
	int mar; /* number of rows in the constraint matrix */
	int mae; /* number of free-row, bound, range and rhs vectors */
	int objsen; /* indicates min (1) or max (-1) problem */
	double *objx; /* array containing the objective row coefficients */
	double *rhsx; /* array containing rhs term for each constraint */
	char *senx; /* array containing the sense of each constraint */
    char *typex;  /* array containing the type of each decision variable */
	int *rngcol; /* don't worry about it */
	int *nrowind;/* don't worry about it */
	int *matbeg; /* array of indices to column groups in matval */
	int *matcnt; /* array containing number of elements in groups */
	int *matind; /* array of row numbers of coefficients in matval */
	double *matval; /* array of coefficients in constraint matrix */
	double *bdl; /* array containing the lower bound on each coeff. */
	double *bdu; /* array containing the upper bound on each coeff. */
	int *etype; /* array containing type of each extra rim vector */
	int *enzbeg; /* rim vector stuff... don't worry about it */
	int *enzcnt; /* rim vector stuff... don't worry about it */
	int *enzind; /* rim vector stuff... don't worry about it */
	double *enzval; /* rim vector stuff... don't worry about it */
	char *dataname;/* don't worry about it */
	char *objname;/* don't worry about it */
	char *rhsname;/* don't worry about it */
	char *rngname;/* don't worry about it */
	char *bndname;/* don't worry about it */
	char **cname; /* array of pointers to character strings in cstore */
	char *cstore; /* end to end col. names separated by '\0' */
	char **rname; /* array of pointers to row strings in rstore */
	char *rstore; /* end to end row names separated by '\0' */
	char *estore; /* don't worry about it */
	char **ename; /* don't worry about it */
	int macsz; /* length of objx, matbeg, matcnt, bdl, bdu, cname. */
	int marsz; /* length of rhsx, senx, rname, rngcol. */
	int matsz; /* length of matind, matval */
	int maesz; /* length of nrowind, etype, ename, enzbeg, enzcnt. */
	int enzsz; /* length of enzind and enzval */
	int nzcnt; /* length of non-zero elements in coefficients matrix*/
	unsigned cstorsz; /* length of cstore */
	unsigned rstorsz; /* length of rstore */
	unsigned estorsz; /* length of estore */
	BOOL feaflag; /*added by Yifan to record the feasibility of the sub-problem 08/11/2011 */
	void *lp; /* pointer to the problem returned by CPLEX */
} one_problem;

/**********************************************************************\
** A vector is just an array of values, of whatever size is allocated.
 ** It is assumed that all the elements are stored, consecutively.
 \**********************************************************************/
typedef double *vector;

typedef struct
{
	vector *incumb_x;
	vector *R_Master_pi; /* modified by Yifan 2012.10.05 */
	vector *R_Master_dj; /* modified by Yifan 2012.10.05 */
} batch_incumb_type;

/**************************************************************************\
**   Each cut consists of a scalar value for the right hand side, stored
 ** in _alpha_, a vector of coefficients for the master program's primal
 ** variables, stored in _beta_, and an array of indices to the maximal pi
 ** for each observation of omega, stored in _istar_ (these are references
 ** into sigma).  In order to weight the cuts properly, _cut_obs_ gives
 ** the number of samples on which the given cut was based.  In contrast,
 ** _omega_cnt_ gives the number of *distinct* observations on which the
 ** cut was based (this is also the length of istar).
 **
 **   Cuts which are "loose" when the master is solved will be dropped
 ** periodically, based on the _slack_cnt_ field, which counts the number
 ** of consecutive iterations the cut's constraint was slack.  Finally,
 ** in order to interface with the LP solver, each cut should know what
 ** its row number is in the master constraint matrix, so this is stored
 ** in _row_num_.
 **
 **   Note that all cuts in a given cell will have the same "num_obs" (from
 ** Sen's earlier implementation paper) so this is stored in theta.
 \**************************************************************************/
typedef struct
{
	int *istar;
	int slack_cnt;
	int cell_num;
	int omega_cnt;
	int cut_obs;
	int row_num;
	double alpha;
	double alpha_incumb; /* added by Yifan 04/04/2012 */
	double *beta;
	BOOL subfeaflag; /*added by Yifan*/
	BOOL is_incumbent; /*added by Yifan 02/02/12*/
	double *subobj_omega; /*added by Yifan 02/02/12/ */
	int *subobj_freq; /*added by Yifan 02/02/12/ */
} one_cut, *cut_ptr;

/**************************************************************************\
**   A collection of the single cuts described above is stored here.
 ** The _val_ array holds pointers to cut structures, while the _cnt_
 ** field tells how many cuts are currently stored in _val_.
 \**************************************************************************/
typedef struct
{
	int cnt;
	one_cut **val;
} cut_type, *cuts_ptr;

typedef struct
{
	int b_size;
	cut_type **batch;
} batch_cut_type, *bcuts_ptr;

typedef struct
{
	config_type config;
	omegastuff omegas;
	one_problem * batch_problem;
	int MALLOC;
	batch_incumb_type *batch_incumb;
	vector Obj_lb;
	vector quad_v;
	batch_cut_type *bcuts;
	batch_cut_type *bfcuts;
	batch_cut_type *bfcuts_pool;
	double ck[BATCH_SIZE];
	double Eta0;
	FILE *fptrALLOC;
	FILE *fptrFREE;
    FILE *fptrOMEGA;
	sd_long MEM_USED;
	clock_t LAST_CLOCK;
	double Abar;
	vector Bbar;
    int average_flag;
    int obj_flag;
    double obj_mean;
    double obj_stdev;
} sdglobal_type;

#endif /* SDGLOBAL_H_ */
