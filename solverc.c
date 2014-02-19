/***********************************************************************\
**
 ** solver.c
 **
 ** This file contains all the routines necessary to solve LP's.
 ** Currently, it uses the CPLEX optimizer.  If the LP solver is
 ** to be changed, all the functions in this file, and all the
 ** typedef's in solver.h, should be re-written.  All other SD
 ** functions use exclusively the functions specified here to
 ** access the LP solver.
 **
 ** You might consider dropping this file altogether... or just
 ** include it as a list of descriptions of CPLEX functions which
 ** another programmer can use to write the appropriate functions.
 **
 ** open_Solver()  	       // initialize CPLEX environment. zl
 ** close_Solver()               // release the CPLEX environment. zl
 ** setup_problem()  
 ** print_problem()
 ** solve_problem()
 ** remove_problem()
 ** change_coef()
 ** change_column()
 ** get_dual()
 ** get_dual_slacks()    
 ** get_lb()           // get the lower bounds of variables. zl 
 ** add_row()
 **
 **
 \**********************************************************************/

/* Most CPLEX 7.0 functions return a status value of zero on success, and a non-zero integer if an error occurs. While on the other hand, all the other part of
 this software treat a returned status value of zero as an indication of failure.So one must be careful on this when updating this software. modified. zl
 */

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "log.h"
#include "sdglobal.h"


/***********************************************************************\
** Release the CPLEX environment.
 \***********************************************************************/
void close_Solver(void)
{

	/* Free up the CPLEX environment, if necessary. */

	int status;

#ifdef TRACE
	printf("Inside close_Solver\n");
#endif

	if (env != NULL)
	{
		status = CPXcloseCPLEX(&env);

		if (status)
		{
			char errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

#ifdef TRACE
	printf("Exiting close_Solver\n");
#endif

}

/***********************************************************************\
** Initialize the CPLEX environment. zl
 \***********************************************************************/
void open_Solver(void)
{
	int status;

#ifdef TRACE
	printf("Inside open_Solver\n");

#endif

	/* Initialize the CPLEX environment */

	env = CPXopenCPLEX(&status);

	/* If an error occurs, the status value indicates the reason for
	 failure. */

	if (env == NULL)
	{
		char errmsg[1024];
		fprintf(stderr, "Could not open CPLEX environment. \n");
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "%s", errmsg);
		goto TERMINATE;
	}

	/* Turn on output to the screen */

	status = set_intparam(NULL, PARAM_SCRIND, ON);
	if (status)
	{
		fprintf(stderr, "Failed to turn on screen indicator, error %d.\n",
				status);
		goto TERMINATE;
	}

	/* Require memory growth. */
	/* These parameters have been removed from CPLEX 10.0  Yifan 06/10/2011 */
	/*status = CPXsetintparam (env, CPX_PARAM_COLGROWTH, 1);
	 if (status) {
	 fprintf (stderr, 
	 "Failed to require variable (column) memory growth, error %d.\n", status);
	 goto TERMINATE;
	 }
	 
	 status = CPXsetintparam (env, CPX_PARAM_ROWGROWTH, sd_global->config.MAX_ITER);
	 if (status) {
	 fprintf (stderr, 
	 "Failure to require constraint (row) memory growth, error %d.\n", status);
	 goto TERMINATE;
	 }
	 
	 status = CPXsetintparam (env, CPX_PARAM_NZGROWTH, 100*sd_global->config.MAX_ITER);
	 if (status) {
	 fprintf (stderr, 
	 "Failure to require non-zero element memory growth, error %d.\n", status);
	 goto TERMINATE;
	 }
	 */

	TERMINATE:

	/* Free up the CPLEX environment, if necessary. */
	if (status != 0)
		close_Solver();
}

/***********************************************************************\ 
 ** Load the LP problem into CPLEX 
 \***********************************************************************/
BOOL setup_problem(one_problem *current)
{
	int status;

#ifdef TRACE
	printf("Inside setup_problem\n");
#endif

	/* A returned pointer of NULL may indicate a memory allocation error. */

	/* modified. zl

	 if ((current->lp = (void *) loadprob(
	 current->name, 
	 current->mac, 
	 current->mar,
	 0, 
	 current->objsen,
	 current->objx,
	 current->rhsx, 
	 current->senx, 
	 current->matbeg,
	 current->matcnt, 
	 current->matind, 
	 current->matval,
	 current->bdl, 
	 current->bdu, 
	 NULL,
	 NULL, 
	 NULL, 
	 NULL, 
	 NULL, 
	 NULL, 
	 NULL,
	 NULL, 
	 current->objname, 
	 NULL,
	 NULL, 
	 NULL, 
	 current->cname,
	 current->cstore, 
	 current->rname, 
	 current->rstore,
	 NULL, 
	 NULL, 
	 current->macsz,
	 current->marsz, 
	 current->matsz, 
	 0, 
	 0,
	 current->cstorsz, 
	 current->rstorsz, 
	 0)))
	 return 1;
	 else
	 {
	 err_msg("Loading", "setup_problem", "current->lp");
	 return 0;
	 }
	 */

	current->lp = CPXcreateprob(env, &status, current->name);

	if ((status = CPXcopylpwnames(env, current->lp, current->mac, current->mar,
			current->objsen, current->objx, current->rhsx, current->senx,
			current->matbeg, current->matcnt, current->matind, current->matval,
			current->bdl, current->bdu, NULL, current->cname, current->rname)))
	{
		err_msg("Loading", "setup_problem", "current->lp");
		return 0;
	}

	else
		return 1;

}

/***********************************************************************\
 ** Tell the external Solver to create a new problem and read data into
 ** that problem from an MPS file.
 \***********************************************************************/
void *read_problem(one_problem *p, char *filename, char *filetype)
{
  /* changed from _void * lp_ to _CPXLPptr lp_ Yifan 11/08/2011 */
  CPXLPptr lp;
  int status;
  
  printf("the filename is : %s\n", filename);
  lp = CPXcreateprob(env, &status, p->name);
  
  if (lp == NULL)
  {
    printf(" read_problem: original wasn't created! \n");
    return lp;
  }
  
  /* the pointer in the second argument is changed from _p->lp_ to _lp_ Yifan 11/08/2011 */
  status = CPXreadcopyprob(env, lp, filename, filetype);
  
  if (status != 0)
  {
    printf(" read_problem: read MPS file failed \n");
    CPXfreeprob(env, &lp);
    lp = NULL;
  }
  
  return lp;
}

/***********************************************************************\
** Tell CPLEX to write out an MPS file for the LP problem,
 ** using the filename provided.  It returns whatever CPLEX returns.
 \***********************************************************************/
BOOL print_problem(one_problem *p, char *filename)
{
	int status;
#ifdef TRACE
	printf("Inside print_problem\n");
#endif

	status = CPXwriteprob(env, p->lp, filename, "LP");

#ifdef TRACE
	printf("Exiting print_problem\n");
#endif
	/*  return mpswrite(p->lp, filename);  modified. zl */
	return status;

	/* return CPXwriteprob(env, p->lp, filename, "LP"); zl */

}

/***********************************************************************\
** This function solves the given LP problem by calling on the
 ** optimizer.  If the optimizer was unable to solve the problem,
 ** it returns FALSE; otherwise TRUE.  Note that no answers are
 ** provided, but must be called for with a separate function.
 \***********************************************************************/
#undef PRINT_X
BOOL solve_problem(sdglobal_type* sd_global, one_problem *p)
{
	int status;
	BOOL ans;
	int qpnzlim;
	int ctr = 0;

#ifdef PRINT_X
	int i, stat1;
	vector x;
#endif
	int param = 1;
#ifdef TRACE
	printf("Inside solve_problem\n");
#endif

#ifdef SAVE
	print_problem(p, "solve.mps");
#endif

#ifdef PRINT_X
	x = (double *) malloc (p->mac+1 * sizeof(double));
	if (!(x))
	{
		printf("ERROR: out of memory, solve_prob\n");
		exit(1);
	}
#endif

	/*status = CPXsetdblparam (env, CPX_PARAM_TILIM, 30000.0);*/
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);

	/* Modified for the regularized QP method (LP=0, QP=1). zl */
	//printf("p->name is: %s \n",p->name);
	if ((strcmp(p->name, "Subproblem") == 0))
	{
		/*turn_off_presolve: if(p->feaflag == FALSE && pre_solve == TRUE) {pre_stat = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF); pre_solve = FALSE;}*/
		change_solver_primal(p);
		CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
		ans = !CPXlpopt(env, p->lp);
		CPXsetintparam(env, CPX_PARAM_PREIND, CPX_ON);
	}
	else if (sd_global->config.MASTER_TYPE == 0)
	{
		ans = !CPXprimopt(env, p->lp);
	}
	else
	{
		resolve_master: ans = !CPXbaropt(env, p->lp);
		/* Set it back to default */
		CPXsetintparam(env, CPX_PARAM_SCAIND, 0);
        change_barrier_algorithm(p, 0); /* Change Barrier Algorithm to its default setting*/
	}

	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	CPXsolution(env, p->lp, &status, NULL, NULL, NULL, NULL, NULL);
	/*Notice!! Won't work with CPXprimopt   added by Yifan 09/13/2011*/

	/* Error messages for checking purpose. added by zl. */
	/* Change from cplex70 to cplex81. zl 04/13/05. */
	/* if (status != CPX_STAT_OPTIMAL) */
	if (status != CPX_STAT_OPTIMAL)
	{
		if (strcmp(p->name, "Subproblem") == 0)
		{ //printf("\nSubproblem non-optimal. ");
		  //printf("Please check complete recourse condition. \n");
		  //printf("CPLEX solution status = %d in solve_problem for %s.\n", 
		  //       status, p->name);
			if (status == CPX_STAT_INFEASIBLE)
			{
				p->feaflag = FALSE; /*added by Yifan to generate feasibility cut 08/11/2011*/
				/*if (pre_solve == TRUE) {
				 goto turn_off_presolve;
				 }*/
				printf("*****   SUB PROBLEM INFEASIBLE   *****\n");
				return 1; //modified by Yifan to test infeasibility of sub-prob
			}
			else
			{
				printf("*****   SUB PROBLEM NO SOLUTION   *****\n");
				return 1;
			}

#ifdef PRINT_X  
			get_primal(x, p, p->mac);
			printf("X :: ");
			for (i=0; i<p->mac; i++)
			printf("%s = %f ", p->cname[i], x[i+1]);
			printf("\n");
			print_contents(p, "solve_err.out");
			print_problem(p, "solve_error.mps");
#endif

		}
		else
		{
			/* Change from cplex70 to cplex81. zl 04/13/05. */
			/* if (status != CPX_STAT_INFEASIBLE) { */
			if (status != CPX_STAT_INFEASIBLE)
			{
				/*printf("CPLEX solution status = %d in solve_problem for %s.\n", status, p->name);*/
				if (status == 6)
				{
					/* Included for aggresive scaling by GJH 04/19/11*/
					/* Yifan 03/25/2012 Solution Status 6 encountered. Change to Equilibration Scaling*/
					/* printf("changing scaling strategy\n"); */
					ctr++;
					param = -param;
					CPXsetintparam(env, CPX_PARAM_SCAIND, param);
					CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
					/* Set 0 if you want to change algorithm immediately after status 6 */
					/* Set 2 if you want to change scaling strategy twice */
					if (ctr > 2)
					{
						/*err_msg("Change to another SEED", "solve_problem", "ill_conditioned problem");*/
						change_barrier_algorithm(p, 2); /* Change Barrier Algorithm to Infeasibility-constant start*/
					}
					goto resolve_master;
				}

				if (38 == status)
				{
					qpnzlim = get_qp_nzreadlim();
					printf("CPX_PARAM_QPNZREADLIM = %d\n", qpnzlim);
					print_problem(p, "mast_err.mps");
					exit(1);
				}
			}
			else
			{
				printf("\nMaster problem infeasible. \n");
				print_problem(p, "Infeasible_Master");
				printf("CPLEX solution status = %d in solve_problem for %s.\n",
						status, p->name);
				printf("Please check the problem data files.\n");
#ifdef PRINT_X  
				get_primal(x, p, p->mac);
				printf("X :: ");
				for (i=0; i<p->mac; i++)
				printf("%s = %f ", p->cname[i], x[i+1]);
				printf("\n");
				print_contents(p, "solve_err.out");
				print_problem(p, "solve_error.mps");
#endif
				exit(1);
			}
		}
	}

#ifdef PRINT_X
	if (x)
	free(x);
#endif

#ifdef TRACE
	printf("Exiting solve_problem\n");
#endif

	return ans;
}


/***********************************************************************\
 ** This function unloads the given LP problem from the
 ** optimizer's memory.  The data and arrays in the one_problem
 ** structure remain intact, and are returned to their original
 ** values.  They must be freed by the user.
 \***********************************************************************/
void remove_problem(one_problem *p)
{
  
  int status;
#ifdef TRACE
  printf("Inside remove_problem\n");
#endif
  
  /* modified.  zl
   unloadprob(&p->lp);
   */
  
  /* Free up the LP problem, if necessary. */
  
  if (p->lp != NULL)
  {
    status = CPXfreeprob(env, (CPXLPptr *) &(p->lp));
    if (status)
    {
      fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
    }
  }
  
#ifdef TRACE
  printf("Exiting remove_problem\n");
#endif
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the number of columns in the argument problem.
 \**********************************************************************************/
int get_numcols(one_problem *p)
{
  return CPXgetnumcols(env, p->lp);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the number of rows in the argument problem.
 \**********************************************************************************/
int get_numrows(one_problem *p)
{
  return CPXgetnumrows(env, p->lp);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the number of nonzeroes in the argument problem.
 \**********************************************************************************/
int get_numnz(one_problem *p)
{
  return CPXgetnumnz(env, p->lp);
}


/***********************************************************************\
** This function retreives and returns the value of the objective function 
 ** for a given problem.  It assumes the problem has already been solved!
 \***********************************************************************/
double get_objective(one_problem *p)
{
	double ans;
	int status;

#ifdef TRACE
	printf("Inside get_objective\n");
#endif

	/* Get the objective function value only */

	CPXsolution(env, p->lp, &status, &ans, NULL, NULL, NULL, NULL);

	/* modified.  zl
	 solution(p->lp, &status, &ans, NULL, NULL, NULL, NULL);
	 */

#ifdef TRACE
	printf("Exiting get_objective\n");
#endif

	return ans;
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the coefficients of the selected columns in
 ** the argument problem, in the arrays cmatbeg, cmatind, cmatval.
 \**********************************************************************************/
int get_obj(one_problem *p, double *obj, int begin, int end)
{
  return CPXgetobj(env, p->lp, obj, begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the coefficients of the selected columns in
 ** the argument problem, in the arrays cmatbeg, cmatind, cmatval.
 \**********************************************************************************/
int get_cols(one_problem *p, int *pnzcnt, int *cmatbeg, int *cmatind,
             double *cmatval, int cmatspace, int *psurplus, int begin, int end)
{
  return CPXgetcols(env, p->lp, pnzcnt, cmatbeg, cmatind, cmatval, cmatspace,
                    psurplus, begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the coefficients of the selected rows in
 ** the argument problem, in the arrays rmatbeg, rmatind, rmatval.
 \**********************************************************************************/
int get_rows(one_problem *p, int *pnzcnt, int *rmatbeg, int *rmatind,
             double *rmatval, int rmatspace, int *psurplus, int begin, int end)
{
  return CPXgetrows(env, p->lp, pnzcnt, rmatbeg, rmatind, rmatval, rmatspace,
                    psurplus, begin, end);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the right hand side values of the selected rows
 ** in the argument problem, in the array rhsx.
 \**********************************************************************************/
int get_rhs(one_problem *p, double *rhsx, int begin, int end)
{
  return CPXgetrhs(env, p->lp, rhsx, begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the "sense" of the selected rows in the
 ** argument problem, in the array senx.
 \**********************************************************************************/
int get_sense(one_problem *p, char *senx, int begin, int end)
{
  return CPXgetsense(env, p->lp, senx, begin, end);
}

/**********************************************************************************\
 ** 2013.09.12: This function returns the "type" of the columns (decision variables)
 ** in the master problem.
 \**********************************************************************************/
int get_type(one_problem *p, char *typex, int begin, int end)
{
    return CPXgetctype (env, p->lp, typex, begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the basis (status of either columns or rows)
 ** in the argument problem, in the arrays cstat and rstat.
 \**********************************************************************************/
int get_basis(one_problem *p, int *cstat, int *rstat)
{
  return CPXgetbase(env, p->lp, cstat, rstat);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns variable lower bounds for the argument problem.
 ** Added as a "clean wrapper" for CPXgetlb(), slightly different from get_lb().
 \**********************************************************************************/
int get_lbound(one_problem *p, double *lb, int begin, int end)
{
  return CPXgetlb(env, p->lp, lb, begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns variable upper bounds for the argument problem.
 ** Added as a "clean wrapper" for CPXgetub(), slightly different from get_ub().
 \**********************************************************************************/
int get_ubound(one_problem *p, double *ub, int begin, int end)
{
  return CPXgetub(env, p->lp, ub, begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function returns the variable values from the argument problem.
 ** Added as a "clean wrapper" for CPXgetx(), slightly different from get_primal().
 \**********************************************************************************/
int get_x(one_problem *p, double * x, int begin, int end)
{
  return CPXgetx(env, p->lp, x, begin, end);
}




/***********************************************************************\
** This function obtains the vector of optimal primal variables for
 ** a given problem.  Note that, like most other vectors,
 ** the solution vector reserves the 0th location for its 1-norm, and
 ** the _length_ parameter is assumed not to include it.
 ** It returns 1 if the query was successful; 0 otherwise.
 \***********************************************************************/
BOOL get_primal(vector X, one_problem *p, int length)
{
	BOOL failed;

#ifdef TRACE
	printf("Inside get_primal\n");
#endif

	/* modified.  zl
	 failed = getx(p->lp, X+1, 0, length-1);
	 */

	failed = CPXgetx(env, p->lp, X + 1, 0, length - 1);
	X[0] = one_norm(X + 1, length);

#ifdef DEBUG
	print_vect(X, length, "Primal solution");
#endif

#ifdef TRACE
	printf("Exiting get_primal\n");
#endif

	return !failed;
}

/***********************************************************************\
** This function obtains the vector of optimal dual variables for
 ** a given problem.  The vector is allocated and filled here, and
 ** must be freed by the user.  Note that, like most other vectors,
 ** the dual vector reserves the 0th location for its 1-norm, and it is
 ** assumed that _length_ does not include the 0th position.
 ** It returns TRUE if the query was successful; FALSE otherwise.
 \***********************************************************************/
BOOL get_dual(vector Pi, one_problem *p, num_type *num, int length)
{
	BOOL failed;

#ifdef TRACE
	printf("Inside get_dual\n");
#endif

	/* modified.  zl
	 failed = getpi(p->lp, Pi+1, 0, length-1);
	 */
	if (strcmp(p->name, "batch_mean") == 0)
	{
		/*This part is for original master constraints*/
		failed = CPXgetpi(env, p->lp, Pi + 1, num->mast_rows * num->batch_id,
				num->mast_rows * num->batch_id + length - 1);
		/*This part is for final optimality cuts*/
		failed = CPXgetpi(
				env,
				p->lp,
				Pi + 1 + length,
				BATCH_SIZE * num->mast_rows + num->batch_id * num->max_cuts,
				BATCH_SIZE * num->mast_rows + num->batch_id * num->max_cuts
						+ num->max_cuts - 1);
	}
	else
	{
		failed = CPXgetpi(env, p->lp, Pi + 1, 0, length - 1);
	}

	Pi[0] = one_norm(Pi + 1, length);

#ifdef TRACE
	printf("Exiting get_dual\n");
#endif

	return !failed;
}

/***********************************************************************\
** This function obtains the vector of optimal dual slacks  for
 ** a given problem.  The vector is allocated and filled here, and
 ** must be freed by the user.  Note that, like most other vectors,
 ** the dual vector reserves the 0th location for its 1-norm, and it is
 ** assumed that _length_ does not include the 0th position.
 ** It returns TRUE if the query was successful; FALSE otherwise.
 \***********************************************************************/
BOOL get_dual_slacks(vector Dj, one_problem *p, num_type *num, int length)
{
	BOOL failed;

#ifdef TRACE
	printf("Inside get_dual_slacks\n");
#endif

	/* modified.  zl
	 failed = getdj(p->lp, Dj+1, 0, length-1);
	 */

	if (strcmp(p->name, "batch_mean") == 0)
	{
		failed = CPXgetdj(env, p->lp, Dj + 1, num->mast_cols * num->batch_id,
				num->mast_cols * num->batch_id + length - 1);
	}
	else
	{
		failed = CPXgetdj(env, p->lp, Dj + 1, 0, length - 1);
	}
	Dj[0] = one_norm(Dj + 1, length);

#ifdef TRACE
	printf("Exiting get_dual_slacks\n");
#endif

	return !failed;
}

/**********************************************************************************\
 ** 2011.10.30: This function gets the character string name of the objective in
 ** the argument problem.
 \**********************************************************************************/
int get_objname(one_problem *p, char *buf, int bufspace, int *psurplus)
{
  return CPXgetobjname(env, p->lp, buf, bufspace, psurplus);
}

/**********************************************************************************\
 ** 2011.10.30: This function gets the character string names of the selected rows in
 ** the argument problem.
 \**********************************************************************************/
int get_rowname(one_problem *p, char **name, char *namestore, int storespace,
                int *psurplus, int begin, int end)
{
  return CPXgetrowname(env, p->lp, name, namestore, storespace, psurplus,
                       begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function gets the character string names of the selected columns
 ** in the argument problem.
 \**********************************************************************************/
int get_colname(one_problem *p, char **name, char *namestore, int storespace,
                int *psurplus, int begin, int end)
{
  return CPXgetcolname(env, p->lp, name, namestore, storespace, psurplus,
                       begin, end);
}

/**********************************************************************************\
 ** 2011.10.30: This function changes the type (LP, QP, etc) of the argument problem.
 \**********************************************************************************/
int change_probtype(one_problem *p, int type)
{
  return CPXchgprobtype(env, p->lp, type);
}

/**********************************************************************************\
 ** 2011.10.30: This function changes one or more coefficients of the objective in
 ** the argument problem.
 \**********************************************************************************/
int change_objective(one_problem *p, int cnt, int *indices, double *values)
{
  return CPXchgobj(env, p->lp, cnt, indices, values);
}

/*********************************************************************************\
 ** 2011.10.30: This function changes one or more right hand sides of constraints/rows
 ** in the argument problem.
 \*********************************************************************************/
int change_rhside(one_problem *p, int cnt, int *indices, double *values)
{
  return CPXchgrhs(env, p->lp, cnt, indices, values);
}

/**********************************************************************************\
 ** 2011.10.30: This function changes one or more bounds on variables/columns in the
 ** argument problem.
 \**********************************************************************************/
int change_bound(one_problem *p, int cnt, int *indices, char *lu, double *bd)
{
  return CPXchgbds(env, p->lp, cnt, indices, lu, bd);
}

/***********************************************************************\
** This function will change the constraint, objective, or right hand
 ** side coefficients of a given problem.  The coefficients are specified
 ** in the form of a sparse matrix, whose rows and columns correspond
 ** to the rows and columns in the problem (objective row is -1; rhs
 ** column is -1).  It returns FALSE if one or more of the coefficients
 ** could not be changed; TRUE otherwise.
 \***********************************************************************/
BOOL change_coef(one_problem *p, sparse_matrix *coef)
{
	int cnt;

#ifdef TRACE
	printf("Inside change_coef\n");
#endif

	/*  modified.  zl
	 for (cnt = 0; cnt < coef->cnt; cnt++)
	 if (chgcoef(p->lp, coef->row[cnt], coef->col[cnt], coef->val[cnt]))
	 return FALSE;
	 */

	for (cnt = 0; cnt < coef->cnt; cnt++)
		if (CPXchgcoef(env, p->lp, coef->row[cnt], coef->col[cnt],
				coef->val[cnt]))
			return FALSE;

#ifdef TRACE
	printf("Exiting change_coef\n");
#endif

	return TRUE;
}

/***********************************************************************\
** This function will change the coefficients of one column of the
 ** constraint matrix (or the right hand side if column is specified as -1)
 ** The non-zero coefficients of the new column are specified in the
 ** _coef_ vector, and are assumed to be contiguous from _start_ to _stop_.
 ** It returns FALSE if one or more of the coefficients could not be
 ** changed; TRUE otherwise.
 \***********************************************************************/
BOOL change_col(one_problem *p, int column, vector coef, int start, int stop)
{
	int row;

#ifdef TRACE
	printf("Inside change_col\n");
#endif

#ifdef DEBUG
	printf("Changing the following coefficients in column %d.\n", column);
	for (row = start; row < stop; row++)
	printf("%d:%f, ", row, coef[row-start]);
	printf("\n");
#endif 

	/* modified.  zl
	 for (row = start; row < stop; row++)
	 if (chgcoef(p->lp, row, column, coef[row-start]))
	 return FALSE;
	 */

	for (row = start; row < stop; row++)
	{
		if (CPXchgcoef(env, p->lp, row, column, coef[row - start]))
			return FALSE;
	}

#ifdef TRACE
	printf("Exiting change_col\n");
#endif

	return TRUE;
}

/***********************************************************************\
** Like change_col(), this function will change the coefficients of
 ** one row of the constraint matrix (or the ojective function if the
 ** row is specified as -1). The non-zero coefficients of the new row
 **  are specified in the _coef_ vector, and are assumed to be contiguous
 ** from _start_ to _stop_. It returns FALSE if one or more of the
 ** coefficients could not be changed; TRUE otherwise.
 \***********************************************************************/
BOOL change_row(one_problem *p, int row, vector coef, int start, int stop)
{
	int col;

#ifdef TRACE
	printf("Inside change_row\n");
#endif

	/* modified.  zl
	 for (col = start; col < stop; col++)
	 if (chgcoef(p->lp, row, col, coef[col-start]))
	 return FALSE;
	 */

	for (col = start; col < stop; col++)
		if (CPXchgcoef(env, p->lp, row, col, coef[col - start]))
			return FALSE;

#ifdef TRACE
	printf("Exiting change_row\n");
#endif

	return TRUE;
}

/***********************************************************************\
** This function inserts a new row (constraint) into a given problem.
 ** It requires an array of all the non-zero coefficients in the row
 ** and their column location.  It returns TRUE if the row was
 ** successfully added; FALSE otherwise.
 ** Perhaps the new row MUST have a name... 
 **
 ** Note: any scalar value which passed to CPLEX is first copied into
 ** an array of length 1.  (assuming CPLEX was expecting an array)
 **
 ** This function has been very hacked, in order to comply with CPLEX.
 \***********************************************************************/
BOOL add_row(one_problem *p, int start, int stop, int *coef_col, double *coef,
		char sense, double yrhs)
{
	int ans;
	static int cumul_num = 0;

	char **r_names;
	double xrhs[1];
	char xsense[1] =
	{ 'G' };
	int rmatbeg[1] =
	{ 0 };

#ifdef DEBUG
	int idx;
#endif

#ifdef TRACE
	printf("Inside add_row\n");
#endif

	/* Give the cut a name */
	r_names = arr_alloc(2, string);
	r_names[0] = arr_alloc(NAME_SIZE, char);

	strcpy(r_names[0], "Cut    ");
	r_names[0][3] = '0' + cumul_num / 10000 % 10;
	r_names[0][4] = '0' + cumul_num / 1000 % 10;
	r_names[0][5] = '0' + cumul_num / 100 % 10;
	r_names[0][6] = '0' + cumul_num / 10 % 10;
	r_names[0][7] = '0' + cumul_num / 1 % 10;
	cumul_num++;

#ifdef DEBUG
	printf("p=%d, nzcnt=%d, rhs=%f, sense=%c, beg[0]=%d, cname=%s, rname=%s.\n",
			p->lp, stop-start+1, yrhs, sense, start, NULL, r_names[0]);
	printf("rmatind:");
	for (idx = 0; idx < stop-start+1; idx++)
	printf(" %d", coef_col[idx]);
	printf("\n");
	printf("rmatval:");
	for (idx = 0; idx < stop-start+1; idx++)
	printf(" %f", coef[idx]);
	printf("\n");
#endif

#ifdef SAVE
	print_problem(p, "addrow.mps");
	print_contents(p, "addrow.out");
	printf("Problem successfully printed\n");
#endif

	/*
	 ** Add zero new columns and one new row to the constraint matrix,
	 ** according to the _coef_ vector.
	 */
	xrhs[0] = yrhs;
	xsense[0] = sense;
	rmatbeg[0] = start;

	/* modified.  zl
	 ans = addrows (p->lp, 0, 1, stop-start+1, xrhs, xsense, rmatbeg,
	 coef_col, coef, NULL, r_names);
	 */

	ans = CPXaddrows(env, p->lp, 0, 1, stop - start + 1, xrhs, xsense, rmatbeg,
			coef_col, coef, NULL, r_names);

	mem_free(r_names[0]);
	mem_free(r_names);

#ifdef SAVE
	print_problem(p, "afteraddrow.mps");
	print_contents(p, "afteraddrow.out");
	printf("Problem successfully printed\n");
#endif

#ifdef TRACE
	printf("Exiting add_row: %d\n", ans);
#endif

	return (!ans);
}

BOOL remove_row(one_problem *p, int row_num)
{
  
#ifdef TRACE
  printf("Inside remove_row\n");
#endif
  
  /* modified.  zl
   return !delrows(p->lp, row_num, row_num);
   */
  
  return !CPXdelrows(env, p->lp, row_num, row_num);
  
}

/***********************************************************************\
 \***********************************************************************/
BOOL add_row_to_master(one_problem *p, int start, int stop, int *coef_col,
		double *coef, char sense, double yrhs)
{
	int ans;
	static int cumul_num = 0;

	char **r_names;
	double xrhs[1];
	char xsense[1] =
	{ 'G' };
	int rmatbeg[1] =
	{ 0 };

#ifdef DEBUG
	int idx;
#endif

#ifdef TRACE
	printf("Inside add_row_to_master\n");
#endif

	/* Give the feasibilty cut a row name */
	r_names = arr_alloc(2, string);
	r_names[0] = arr_alloc(NAME_SIZE, char);

	strcpy(r_names[0], "FeaCut    ");
	r_names[0][3] = '0' + cumul_num / 10000 % 10;
	r_names[0][4] = '0' + cumul_num / 1000 % 10;
	r_names[0][5] = '0' + cumul_num / 100 % 10;
	r_names[0][6] = '0' + cumul_num / 10 % 10;
	r_names[0][7] = '0' + cumul_num / 1 % 10;
	cumul_num++;

#ifdef DEBUG
	printf("p=%d, nzcnt=%d, rhs=%f, sense=%c, beg[0]=%d, cname=%s, rname=%s.\n",
			p->lp, stop-start+1, yrhs, sense, start, NULL, r_names[0]);
	printf("rmatind:");
	for (idx = 0; idx < stop-start+1; idx++)
	printf(" %d", coef_col[idx]);
	printf("\n");
	printf("rmatval:");
	for (idx = 0; idx < stop-start+1; idx++)
	printf(" %f", coef[idx]);
	printf("\n");
#endif

#ifdef SAVE
	print_problem(p, "AddFeaRow.mps");
	print_contents(p, "AddFeaRow.out");
	printf("Problem successfully printed\n");
#endif

	/*
	 ** Add zero new columns and one new row to the constraint matrix,
	 ** according to the _coef_ vector.
	 */
	xrhs[0] = yrhs;
	xsense[0] = sense;
	rmatbeg[0] = start;

	ans = CPXaddrows(env, p->lp, 0, 1, stop - start + 1, xrhs, xsense, rmatbeg,
			coef_col, coef, NULL, r_names);

	mem_free(r_names[0]);
	mem_free(r_names);

#ifdef SAVE
	if (cumul_num==2)
	{
		print_problem(p, "AfterAddFeaRow.mps");
		print_contents(p, "AfterAddFeaRow.out");
		printf("Problem successfully printed\n");
	}

#endif

#ifdef TRACE
	printf("Exiting add_row_to_master: %d\n", ans);
#endif

	return (!ans);
}

/***********************************************************************\
 \***********************************************************************/
BOOL add_row_to_batch(one_problem *p, int start, int nzcnt, int *coef_col,
                      double *coef, char sense, double yrhs, int batch_id)
{
  int ans;
  static int cumul_num = 0;
  
  char **r_names;
  double xrhs[1];
  char xsense[1] =
  { 'G' };
  int rmatbeg[1] =
  { 0 };
  
#ifdef DEBUG
  int idx;
#endif
  
#ifdef TRACE
  printf("Inside add_row_to_batch\n");
#endif
  
  /* Give the feasibilty cut a row name */
  r_names = arr_alloc(2, string);
  r_names[0] = arr_alloc(NAME_SIZE, char);
  
  strcpy(r_names[0], "BCUT     ");
  r_names[0][4] = '0' + cumul_num / 10000 % 10;
  r_names[0][5] = '0' + cumul_num / 1000 % 10;
  r_names[0][6] = '0' + cumul_num / 100 % 10;
  r_names[0][7] = '0' + cumul_num / 10 % 10;
  r_names[0][8] = '0' + cumul_num / 1 % 10;
  r_names[0][9] = 'B';
  r_names[0][10] = '0' + batch_id / 10 % 10;
  r_names[0][11] = '0' + batch_id / 1 % 10;
  
  cumul_num++;
  
#ifdef DEBUG
  printf("p=%d, nzcnt=%d, rhs=%f, sense=%c, beg[0]=%d, cname=%s, rname=%s.\n",
         p->lp, nzcnt, yrhs, sense, start, NULL, r_names[0]);
  printf("rmatind:");
  for (idx = 0; idx < nzcnt; idx++)
	printf(" %d", coef_col[idx]);
  printf("\n");
  printf("rmatval:");
  for (idx = 0; idx < nzcnt; idx++)
	printf(" %f", coef[idx]);
  printf("\n");
#endif
  
#ifdef SAVE
  print_problem(p, "AddFeaRow.mps");
  print_contents(p, "AddFeaRow.out");
  printf("Problem successfully printed\n");
#endif
  
  /*
   ** Add zero new columns and one new row to the constraint matrix,
   ** according to the _coef_ vector.
   */
  xrhs[0] = yrhs;
  xsense[0] = sense;
  rmatbeg[0] = start;
  
  ans = CPXaddrows(env, p->lp, 0, 1, nzcnt, xrhs, xsense, rmatbeg, coef_col,
                   coef, NULL, r_names);
  
  mem_free(r_names[0]);
  mem_free(r_names);
  
#ifdef SAVE
  if (cumul_num==2)
  {
    print_problem(p, "AfterAddFeaRow.mps");
    print_contents(p, "AfterAddFeaRow.out");
    printf("Problem successfully printed\n");
  }
  
#endif
  
#ifdef TRACE
  printf("Exiting add_row_to_master: %d\n", ans);
#endif
  
  return (!ans);
}

/***********************************************************************\
  Write out the problem in the specified format by the file name. zl
 \***********************************************************************/
void write_prob(one_problem *p, char *file_name)
{
	int status;

#ifdef TRACE
	printf("inside write_prob\n");
#endif

	status = CPXwriteprob(env, p->lp, file_name, NULL);

	if (status)
	{
		fprintf(stderr, "Failed to write the problem -- %s.\n", file_name);
		exit(1);
	}

#ifdef  TRACE
	printf("Exiting write_prob\n");
#endif

}


/**********************************************************************************\
 ** 2011.10.30: This function sets an integer parameter affecting Solver behavior.
 \**********************************************************************************/
int set_intparam(one_problem *p, int whichparam, int newvalue)
{
  return CPXsetintparam(env, whichparam, newvalue);
}

/***********************************************************************\
 Get the the number of Q matrix nonzeros that can be read.
 \***********************************************************************/
int get_qp_nzreadlim(void)
{
  int param = 0;
  
  CPXgetintparam(env, CPX_PARAM_QPNZREADLIM, &param);
  printf("CPX_PARAM_QPREADLIM = %d", param);
  
  return param;
}

/***********************************************************************\
  Set the the number of Q matrix nonzeros that can be read. 
 It returns zero on success. zl, 09/25/05. 
 \***********************************************************************/
int set_qp_nzreadlim(int nzreadlim)
{
	int status = 0;
	int param = 0;

	status = CPXsetintparam(env, CPX_PARAM_QPNZREADLIM, nzreadlim);
	status = CPXgetintparam(env, CPX_PARAM_QPNZREADLIM, &param);
	printf("CPX_PARAM_QPREADLIM = %d", param);

	return status;
}

/**********************************************************************************\
 ** 2011.10.30: This function sets the QP coefficients of the argument problem from
 ** the vector of diagonal elements, for a separable QP.
 \**********************************************************************************/
int copy_qp_separable(one_problem *p, double *qsepvec)
{
  /* NOTE: CPLEX evaluates the corresponding objective with a factor of 0.5
   in front of the quadratic objective term.*/
  return CPXcopyqpsep(env, p->lp, qsepvec);
}

/**********************************************************************************\
 ** 2013.09.12: This function sets the MIX variable type of the master problem.
 \**********************************************************************************/
int copy_type(one_problem *p, char *typex)
{
    return CPXcopyctype (env, p->lp, typex);
}

/***********************************************************************\
** This function obtains the vector of lower bounds (lb) for a given problem.  
 ** Note that, like most other vectors, the lb vector reserves the 0th location 
 ** for its 1-norm, and it is assumed that _length_ does not include the 0th 
 ** position. It returns TRUE if the query was successful; FALSE otherwise.
 \***********************************************************************/
BOOL get_lb(vector lb, one_problem *p, int length)
{
	BOOL failed;

#ifdef TRACE
	printf("Inside get_lb\n");
#endif

	failed = CPXgetlb(env, p->lp, lb + 1, 0, length - 1);
	lb[0] = one_norm(lb + 1, length);

#ifdef TRACE
	printf("Exiting get_lb\n");
#endif

	return !failed;
}

/***********************************************************************\
 ** This function obtains the vector of upper bounds (ub) for a given problem.  
 ** Note that, like most other vectors, the ub vector reserves the 0th location 
 ** for its 1-norm, and it is assumed that _length_ does not include the 0th 
 ** position. It returns TRUE if the query was successful; FALSE otherwise.
 \***********************************************************************/
BOOL get_ub(vector ub, one_problem *p, int length)
{
	BOOL failed;

#ifdef TRACE
	printf("Inside get_ub\n");
#endif

	failed = CPXgetub(env, p->lp, ub + 1, 0, length - 1);
	ub[0] = one_norm(ub + 1, length);

#ifdef TRACE
	printf("Exiting get_ub\n");
#endif

	return !failed;
}

/**********************************************************************************\
 ** 20112.04.26: This function solves an lp problem with simplex algorithm
 \**********************************************************************************/
int solve_lp(one_problem *p)
{
  int status;
  CPXlpopt(env, p->lp);
  status = CPXsolution(env, p->lp, NULL, NULL, NULL, NULL, NULL, NULL);
  return status;
}

/**********************************************************************************\
 ** 20112.04.26: This function solves an lp problem with simplex algorithm
 \**********************************************************************************/
void *clone_prob(one_problem *p)
{
  void *lp;
  int status = 0;
  lp = CPXcloneprob(env, p->lp, &status);
  if (status != 0 || lp == NULL)
  {
    printf(" clone_prob: clone problem failed \n");
    exit(0);
  }
  return lp;
}

/****************************************************************************\
 This function change the CPLEX LP optimization method to barrier.
 \****************************************************************************/
void change_solver_barrier(one_problem *p)
{
  int status = 0;
  
#ifdef TRACE
  printf("Inside change_solver_barrier.\n");
#endif
  //Change qpmethod to barrier modified by Yifan 06/15/2011
  status = set_intparam(NULL, PARAM_QPMETHOD, ALG_CONCURRENT); /* 2011.10.30 */
  if (status)
  {
    fprintf(stderr, "Failed to set the optimization method, error %d.\n",
            status);
    exit(0);
  }
  
  status = 0;
  
  /*
   -1   No crossover
   0    Automatic: let CPLEX choose; default
   1    Primal crossover
   2    Dual crossover
   */
  
  status = set_intparam(NULL, PARAM_BARCROSSALG, 0); /* 2011.10.30 */
  
  
  if (status)
  {
    fprintf(stderr, "Failed to set the barcrossover off, error %d.\n",
            status);
    exit(0);
  }
}

/****************************************************************************\
 This function change the algorithm of barrier from standard to Infeasibility-constant start.
 When you encouter the "status 6", we should consider a different barrier algorithm.
 CPLEX Paremeter: PARAM_BARALG
 0: Default setting
 1: Infeasibility-estimate start
 2: Infeasibility-constant start
 3: Standard barrier
 \****************************************************************************/
void change_barrier_algorithm(one_problem *p, int k)
{
  int status = 0;
  
#ifdef TRACE
  printf("Inside change_barrier_algorithm.\n");
#endif
  
  status = set_intparam(NULL, PARAM_BARALG, k); /* 2011.10.30 */
  
  if (status)
  {
    fprintf(stderr, "Failed to set the barrier algorithm, error %d.\n",
            status);
    exit(0);
  }
}


/****************************************************************************\
 This function change the CPLEX LP optimization method to dual simplex. 
 \****************************************************************************/
void change_solver_primal(one_problem *p)
{
	int status = 0;

#ifdef TRACE
	printf("Inside change_solver_primal.\n");
#endif

	status = set_intparam(NULL, PARAM_LPMETHOD, ALG_AUTOMATIC); /* 2011.10.30 */

	if (status)
	{
		fprintf(stderr, "Failed to set the optimization method, error %d.\n",
				status);
		exit(0);
	}
}

/**********************************************************************************\
 ** 2012.09.10: This function returns the coefficients of the selected position in
 ** the argument problem, in the pointer to a double coefficient.
 \**********************************************************************************/
int get_coef(one_problem *p, int row, int col, double *coef)
{
  return CPXgetcoef(env, p->lp, row, col, coef);
}


