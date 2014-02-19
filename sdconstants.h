/*
 * sdconstants.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef SDCONSTANTS_H_
#define SDCONSTANTS_H_

//utility.c
#define MAX_BITS	sizeof(int) * 8
#define KEY_PAGE	500 / MAX_BITS

//input.c
#define WORDSIZE 16
#define BLOCKSIZE 20
#define STRSIZE 180

//optimal.c
#define  HOOPS 3

//prob.h
/*
 ** "#define" or "#undef" these flags to control debugging print statements.
 */
#define RUN    /* provides streamlined output when running SD, including     */
/* solution & objective estimate every PRINT_CYCLE iterations */
#undef WRITE  /* causes SD output to be written to files, including the     */
/* solution, estimate, structure sizes, and Pi-histograms     */
#undef PIS     /* prints messages about pi-thinning in lambda, sigma, delta  */
#undef OMEG    /* prints messages concerning the thinning of omegas          */
#undef ENCODE  /* prints omega information to check the bitwise encoding     */
#undef OPT     /* prints debug messages in the code for checking optimality  */
#undef DEBUG   /* prints lots of stuff, like structure contents and messages */
#undef TRACE   /* prints a message at the beginning every major function     */
#undef SAVE    /* causes program to save every master & subproblem solved    */
#undef LOOP    /* prints a message at the beginning of minor functions,      */
/* which often occur (many, many times) inside loops          */
#undef TIM     /* controls statements in Tim's code (load_core/stoch/time)   */
#undef LOAD_CORE_CPX  /* added by zl. to check the load_core_cpx routine.    */
#undef LOAD_CORE      /* added by zl. to check the load_core routine.        */
#undef CAL_CHECK      /* prints data/calculations on checking QP. zl */
#undef TEST_SDAPI     /* prints debug messages in API's Yifan 11/02/2011 */
#undef OPT_TEST    /* prints messages about full optimality test Yifan 03/18/2012*/
#undef TEST_SOLN  /* prints candidate and incumbent solutions before and after improvement checking Yifan 03/25/12*/
#undef RECOURSE_OBJ
#undef OMEGA_FILE
#undef REC_OMEGA

/*
 ** Some constants used in various places
 */
#define NAME_SIZE	16	/* length of CPLEX problem names */
#define RHS_COL		-1	/* index of the right hand side in CPLEX */
#define GE		'G'	/* code for a >= constraint in CPLEX */
#define LE		'L'

#define NUM_INTS	25	/* ~ # of integers that can fit on one line */
#define NUM_DBLS	6	/* ~ # of doubles that can fit on one line */

#define DROPPED		-1	/* Code for a dropped Pi in a cut's istar */
#define UNUSED		0	/* Code for an available location in array */
#define USED		1	/* Code for an occupied location in array */
#define SDLP              0       /* Code for MASTER_TYPE as Basic LP zl */
#define SDQP              1       /* Code for MASTER_TYPE as Regularized QP zl */

#define BATCH_SUFFIX 3
#define BUFFER_SIZE 248
/* If you wish to run fewer than 30 replications, the normality property of CLT may not be questionable. */
#define BATCH_SIZE 30

#define ITER_DAT	"iter.dat"
#define CNT_DAT		"cnt.mat"
#define L_HIST_DAT	"lamb_histo.mat"
#define S_HIST_DAT	"sig_histo.mat"

//quad.h
#define ZERO_IDENTITY 1.0e-7

//solver.h

/* External Solver Options*/
#define CPLEX
#undef GUROBI

/************************************************************************\
 **           Parameters definition for different solvers.  2012.04.27 Yifan
 \************************************************************************/
/* If CPLEX is the solver for SD, then define the following parameters 2012.04.27 Yifan */
#ifdef CPLEX
#include "cplex.h"
typedef CPXENVptr ENVptr;

// The following is used in several places including:
// input.c, master.c, optimal.c
#define INFBOUND         CPX_INFBOUND

// The followings are used in quad.c
#define PARAM_QPMETHOD   CPX_PARAM_QPMETHOD
#define PARAM_LPMETHOD   CPX_PARAM_LPMETHOD
#define ALG_BARRIER      CPX_ALG_BARRIER
#define ALG_AUTOMATIC    CPX_ALG_AUTOMATIC
#define ALG_CONCURRENT   CPX_ALG_CONCURRENT
#define PROB_QP          CPXPROB_QP
#define PARAM_BARCROSSALG CPX_PARAM_BARCROSSALG
#define PARAM_BARALG     CPX_PARAM_BARALG

// The followings are used in utility.c
#define AT_LOWER         CPX_AT_LOWER
#define BASIC            CPX_BASIC
#define AT_UPPER         CPX_AT_UPPER
#define FREE_SUPER       CPX_FREE_SUPER

// The following is used in input.c
#define NEGATIVE_SURPLUS CPXERR_NEGATIVE_SURPLUS

// The followings are used in sd.c
#define PARAM_SCRIND CPX_PARAM_SCRIND
#define ON CPX_ON

#define BASIC_SOLN       CPX_BASIC_SOLN
ENVptr env;
#endif

/* If GUROBI is the solver for SD, then define the following parameters 2012.04.27 Yifan */
#ifdef GUROBI
#include "gurobi_c.h"
typedef GRBenv ENVptr;
// The following is used in several places including:
// input.c, master.c, optimal.c
#define INFBOUND         GRB_INFINITY

// The followings are used in quad.c
#define PARAM_QPMETHOD   GRB_INT_PAR_METHOD
#define PARAM_LPMETHOD   GRB_INT_PAR_METHOD
#define ALG_BARRIER      GRB_INT_PAR_CROSSOVER
#define ALG_AUTOMATIC    GRB_METHOD_AUTO
#define PARAM_BARCROSSALG GRB_INT_PAR_CROSSOVER
#define PARAM_CROSSOVERBASIS  GRB_INT_PAR_CROSSOVERBASIS
#define PARAM_BARALG     GRB_INT_PAR_BARORDER
#define PROB_QP 1

// The followings are used in utility.c
#define AT_LOWER         GRB_NONBASIC_LOWER
#define BASIC            GRB_BASIC
#define AT_UPPER         GRB_NONBASIC_UPPER
#define FREE_SUPER       GRB_SUPERBASIC

// The following is used in input.c
#define NEGATIVE_SURPLUS GRB_ERROR_SIZE_LIMIT_EXCEEDED

// The followings are used in sd.c (Gurobi only has presolve parameters)
#define PARAM_SCRIND GRB_INT_PAR_OUTPUTFLAG
#define ON 1

#define BASIC_SOLN       CPX_BASIC_SOLN  //not used

GRBmodel *model_i;
ENVptr *env;
#endif


/************************************************************************\
                        OS Platform Options 
\************************************************************************/
#define SD_unix
#undef SD_win

#ifdef SD_unix
#define sd_long long long
#endif

#ifdef SD_win
#define sd_long long long
#endif

/************************************************************************\
**			Macro Definitions
 \************************************************************************/

/*
 ** This macro swaps the contents of _a_ and _b_, using the temporary
 ** variable _t_.  All parameters must be of the same type.
 */
#define 	swap(a,b,t)     	(t) = (a), (a) = (b), (b) = (t)

/*
 ** Returns the absolute value of a number (not just for int's).
 */
#define		DBL_ABS(x)		((x) > 0.0 ? (x) : -(x))

/*
 ** Returns the square value of a number (not just for int's).
 */
#define		SQR(x)			((x) * (x))

/*
 ** Determines the maximum number of cuts allowed, given the number
 ** of primal variables in the master problem.
 */
#define		MAX_CUTS(c)		(3 * (c) + 3)

/*
 ** Copies the first n+1 elements of array _b_ into array _a_
 ** (the 1-norm stored at 0 is also copied).
 */
#define copy_arr(a,b,n)  {int i;for(i=0;i<=(n);i++) *((a)+i)=*((b)+i);}

/*
 ** Compute the minimum of two numeric values. added by zl. 06/18/02
 */
#define min(X, Y) ((X) < (Y) ? (X) : (Y))
#define max(X, Y) ((X) > (Y) ? (X) : (Y))

/*
 ** Returns a pointer to an array of type _type_ and size _n_.
 ** Uses homemade allocation routines, so logging / accounting may be done.
 */
#define		arr_alloc(n,type)	(type *)mem_calloc((n),sizeof(type))

/************************************************************************\
**		Macro Definitions for dynamic memory allocations
 \************************************************************************/

/*
 ** All dynamic allocation in the SD program uses the macros defined here.
 ** So, if you want to record a log of all allocations and frees performed,
 ** or if you want keep a running total of memory used, you can redefine
 ** these macros and insert the appropriate code.
 */

/*
 ** Returns a pointer to an array of type _type_ and size _n_
 ** Also logs information about the memory allocated.
 **
 #define mem_calloc(n,size) calloc(n,size)
 **
 */
#define mem_calloc(n,size) log_alloc("calloc : " #n " : " #size, \
                                      calloc((n),(size)), ((n) * size))

/*
 ** Returns a pointer to allocated memory of size _n_
 ** Also logs information about the memory allocated.
 **
 #define mem_malloc(n)		malloc(n)
 **
 */
#define mem_malloc(n) log_alloc("malloc : " #n,malloc((n)), (n))

/*
 ** Returns a pointer to reallocated memory of new size _n_ from location _ptr_
 ** Also logs information about the memory allocated.
 **
 #define mem_realloc(ptr, n) realloc(ptr,n)
 **
 */
#define mem_realloc(ptr, n) log_realloc("realloc : " #ptr " : " #n,(ptr), \
                                         realloc((ptr),(n)), (n))

/*
 ** Frees the memory associated with the pointer _ptr_
 ** Also compares information about the memory freed with that allocated
 **
 #define mem_free(ptr) log_free(#ptr,(ptr))
 **
 */
#define mem_free(ptr) free(ptr)

#endif /* SDCONSTANTS_H_ */
