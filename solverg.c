/***********************************************************************\
 **
 ** solver.c
 **
 ** This file contains all the routines necessary to solve LP's.
 ** Currently, it uses the CPLEX optimizer.  If the LP solver is
 ** to be changed, all the functions in this file, and all the
 ** typedef's in prob.h, should be re-written.  All other SD
 ** functions use exclusively the functions specified here to
 ** access the LP solver.
 **
 ** 2011.10.30 - Modified by Frontline Systems to complete the
 ** abstraction/isolation of use of an external LP/QP Solver:
 ** Ensure that ALL calls to the Solver, incuding problem query
 ** calls and file read/write calls, are made through this source
 ** file, and all typedefs and constants related to the external
 ** LP/QP Solver are defined in prob.h.
 **
 ** open_Solver()  	       // initialize Solver environment.
 ** close_Solver()          // release the Solver environment.
 ** setup_problem()  
 ** read_problem()          // 2011.10.30 - added
 ** print_problem()
 ** solve_problem()
 ** remove_problem()
 ** get_numcols()           // 2011.10.30 - added
 ** get_numrows()           // 2011.10.30 - added
 ** get_numnz()             // 2011.10.30 - added
 ** get_objective()         // 2011.10.30 - documented
 ** get_cols()              // 2011.10.30 - added
 ** get_rows()              // 2011.10.30 - added
 ** get_rhs()               // 2011.10.30 - added
 ** get_sense()             // 2011.10.30 - added
 ** get_basis()             // 2011.10.30 - added
 ** get_lbound()            // 2011.10.30 - added
 ** get_ubound()            // 2011.10.30 - added
 ** get_x()                 // 2011.10.30 - added
 ** get_primal()            // 2011.10.30 - documented
 ** get_dual()
 ** get_dual_slacks()    
 ** get_objname()           // 2011.10.30 - added
 ** get_rowname()           // 2011.10.30 - added
 ** get_colname()           // 2011.10.30 - added
 ** change_probtype()       // 2011.10.30 - added
 ** change_objective()      // 2011.10.30 - added
 ** change_rhside()         // 2011.10.30 - added
 ** change_bound()          // 2011.10.30 - added
 ** change_coef()
 ** change_col()            // 2011.10.30 - name fixed
 ** change_row()            // 2011.10.30 - documented
 ** add_row()
 ** remove_row()            // 2011.10.30 - documented
 ** write_prob()            // 2011.10.30 - documented
 ** set_intparam()          // 2011.10.30 - added
 ** get_qp_nzreadlim()      // 2011.10.30 - documented
 ** set_qp_nzreadlim()      // 2011.10.30 - documented
 ** copy_qp_separable()     // 2011.10.30 - added
 ** get_lb()           // get the lower bounds of variables. zl 
 ** get_ub()                // 2011.10.30 - documented
 **
 \**********************************************************************/


#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "solver.h"
#include "utility.h"
#include "log.h"
#include "sdglobal.h"


/***********************************************************************\
 ** Release the Solver environment.
 \***********************************************************************/
void close_Solver()
{
  
  /* Free up the Gurobi environment, if necessary. */
  GRBfreeenv(env);
}


/***********************************************************************\
 ** Initialize the Solver environment. zl
 \***********************************************************************/

void open_Solver()
{ int error = 0;      
  error = GRBloadenv(&env, "lp.log");
  if (error || env == NULL) {
    fprintf(stderr, "Error: could not create environment\n");
    exit(1);
  }
}

/***********************************************************************\ 
 ** Load an already-created LP problem into the external Solver. 
 \***********************************************************************/

BOOL setup_problem(one_problem *current)
{
  int status;

  if ((status = GRBloadmodel(env, (GRBmodel**)(&current->lp), NULL, current->mac, current->mar, current->objsen, 0, current->objx, current->senx, current->rhsx, current->matbeg, current->matcnt, current->matind, current->matval, current->bdl, current->bdu, NULL,current->cname, current->rname)))
  {
    err_msg("Loading", "setup_problem", "current->lp");
    return FALSE;
  }

  return TRUE;
}


/***********************************************************************\ 
 ** Tell the external Solver to create a new problem and read data into
 ** that problem from an MPS file.
 \***********************************************************************/
void *read_problem(one_problem *p, char *filename, char *filetype)
{
   int status;
   
   if((status = GRBreadmodel(env, filename, &model_i)))
   {
    printf(" read_problem: read MPS file failed ; status: %d \n", status);
    exit(1);
   }
  
   return model_i;
}


/***********************************************************************\
 ** Tell the external Solver to write out an MPS file for the LP problem,
 ** using the filename provided.  It returns whatever the Solver returns.
 \***********************************************************************/
BOOL print_problem(one_problem *p, char *filename)
{

  int status;
  
  GRBupdatemodel((GRBmodel*)p->lp);
  
  if ((status = GRBwrite((GRBmodel*)p->lp, filename)))
  {
   err_msg("Loading", "print_problem", "p->lp");
   return FALSE;
  }
  
  return TRUE;
}  



/***********************************************************************\
 ** This function solves the given LP problem by calling on the
 ** optimizer.  If the optimizer was unable to solve the problem,
 ** it returns FALSE; otherwise TRUE.  Note that no answers are
 ** provided, but must be called for with a separate function.
 \***********************************************************************/
BOOL solve_problem(sdglobal_type* sd_global, one_problem *p)
{
  int status;
  BOOL	ans;

  /* Turn off output log 2012.05.02 Yifan */
  status = set_intparam (p, PARAM_SCRIND, 0);
  
  GRBupdatemodel((GRBmodel*)p->lp);
  
  if (GRBoptimize((GRBmodel*)p->lp))
    ans = FALSE;
  else
    ans = TRUE;

  GRBgetintattr((GRBmodel*)p->lp,GRB_INT_ATTR_STATUS, &status);
  
  if (status != GRB_OPTIMAL)
  {
    if (strcmp(p->name, "Subproblem") == 0)
    {
      if (status==GRB_INFEASIBLE) {
        printf("Sub problem INFEASIBLE.\n");
        p->feaflag=FALSE ;  /*added by Yifan to generate feasibility cut 08/11/2011*/
      }
      return status;  //modified by Yifan to test infeasibility of sub-prob
    } 
    else {
      if (status == GRB_SUBOPTIMAL) {
        printf("s");
        /* Change to concurrent solver for solving first stage robustly */
        /* set_intparam (p, PARAM_QPMETHOD, 3); not working??? Ok just print s*/
      }
      else{
        printf("\nMaster problem unsolvable. \n");
        printf("Gurobi solution status = %d in solve_problem for %s.\n", status, p->name);
        printf("Please check the problem data files.\n");
        exit(1);
      }
    }
  }
  return status;
}


/***********************************************************************\
 ** This function unloads the given LP problem from the
 ** optimizer's memory.  The data and arrays in the one_problem
 ** structure remain intact, and are returned to their original
 ** values.  They must be freed by the user.
 \***********************************************************************/
void remove_problem(one_problem *p)
{
  /* Free up the LP problem, if necessary. */
  
  GRBfreemodel(p->lp);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the number of columns in the argument problem.
 \**********************************************************************************/
int get_numcols(one_problem *p)
{
  int ans;
  
  GRBgetintattr((GRBmodel*)p->lp, GRB_INT_ATTR_NUMVARS, &ans);
  return ans;
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the number of rows in the argument problem.
 \**********************************************************************************/
int get_numrows(one_problem *p)
{
  int ans;
  
  GRBgetintattr((GRBmodel*)p->lp, GRB_INT_ATTR_NUMCONSTRS, &ans);
  
  return ans;
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the number of nonzeroes in the argument problem.
 \**********************************************************************************/
int get_numnz(one_problem *p)
{
  int ans;
  
  GRBgetintattr((GRBmodel*)p->lp, GRB_INT_ATTR_NUMNZS, &ans);
  
  return ans;
}


/***********************************************************************\
 ** This function retreives and returns the value of the objective function 
 ** for a given problem.  It assumes the problem has already been solved!
 \***********************************************************************/
double get_objective(one_problem *p)
{
  double ans;
  
  GRBgetdblattr((GRBmodel*)p->lp, GRB_DBL_ATTR_OBJVAL, &ans);
  
  return ans;
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the coefficients of the selected columns in
 ** the argument problem, in the arrays cmatbeg, cmatind, cmatval.
 \**********************************************************************************/
int get_obj(one_problem *p, double *obj, int begin, int end)
{
   return GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_OBJ, begin, end-begin+1, obj);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the coefficients of the selected columns in
 ** the argument problem, in the arrays cmatbeg, cmatind, cmatval.
 \**********************************************************************************/
int get_cols(one_problem *p, int *pnzcnt, int *cmatbeg, int *cmatind, double *cmatval,
             int cmatspace, int *psurplus, int begin, int end)
{
   return GRBgetvars((GRBmodel*)p->lp, pnzcnt,cmatbeg, cmatind, cmatval, begin, end-begin+1);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the coefficients of the selected rows in
 ** the argument problem, in the arrays rmatbeg, rmatind, rmatval.
 \**********************************************************************************/
int get_rows(one_problem *p, int *pnzcnt, int *rmatbeg, int *rmatind, double *rmatval,
             int rmatspace, int *psurplus, int begin, int end)
{
   return GRBgetconstrs((GRBmodel*)p->lp, pnzcnt, rmatbeg, rmatind, rmatval, begin, end-begin+1);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the right hand side values of the selected rows
 ** in the argument problem, in the array rhsx.
 \**********************************************************************************/
int get_rhs(one_problem *p, double *rhsx, int begin, int end)
{
   return GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_RHS, begin, end-begin+1, rhsx);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the "sense" of the selected rows in the 
 ** argument problem, in the array senx.
 \**********************************************************************************/
int get_sense(one_problem *p, char *senx, int begin, int end)
{
   return GRBgetcharattrarray((GRBmodel*)p->lp, GRB_CHAR_ATTR_SENSE, begin, end-begin+1, senx);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the basis (status of either columns or rows)
 ** in the argument problem, in the arrays cstat and rstat.
 \**********************************************************************************/
int get_basis(one_problem *p, int *cstat, int *rstat)
{
  int status;
  status = GRBgetintattr((GRBmodel*)p->lp, GRB_INT_ATTR_VBASIS, cstat);
  status = GRBgetintattr((GRBmodel*)p->lp, GRB_INT_ATTR_CBASIS, rstat);
  return status;
}


/**********************************************************************************\
 ** 2011.10.30: This function returns variable lower bounds for the argument problem.
 ** Added as a "clean wrapper" for CPXgetlb(), slightly different from get_lb().
 \**********************************************************************************/
int get_lbound(one_problem *p, double *lb, int begin, int end)
{
   return GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_LB, begin, end-begin+1, lb);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns variable upper bounds for the argument problem.
 ** Added as a "clean wrapper" for CPXgetub(), slightly different from get_ub().
 \**********************************************************************************/
int get_ubound(one_problem *p, double *ub, int begin, int end)
{
   return GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_UB, begin, end-begin+1, ub);
}


/**********************************************************************************\
 ** 2011.10.30: This function returns the variable values from the argument problem.
 ** Added as a "clean wrapper" for CPXgetx(), slightly different from get_primal().
 \**********************************************************************************/
int get_x(one_problem *p, double * x, int begin, int end)
{
  return GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_X, begin, end-begin+1, x);
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
  BOOL	failed;
  
  failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_X, 0, length, X+1) == 0 ? FALSE: TRUE);
  X[0] = one_norm(X+1, length);

  if (failed)
    return FALSE;
  return TRUE;
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
  
  /* The _len_ argument should be "length" instead of "length-1" */
  if (strcmp(p->name, "batch_mean") == 0)
  {
    /*This part is for original master constraints*/
    failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_PI, num->mast_rows * num->batch_id, length, Pi+1) == 0 ? FALSE : TRUE);
    /*This part is for final optimality cuts*/
    failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_PI, BATCH_SIZE * num->mast_rows + num->batch_id * num->max_cuts, num->max_cuts, Pi + 1 + length) == 0 ? FALSE : TRUE);
  }
  else
  {
    failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_PI, 0, length, Pi+1) == 0 ? FALSE : TRUE);
  }
  
  Pi[0] = one_norm(Pi+1, length);
  
  if (failed)
    return FALSE;
  return TRUE;
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
  BOOL	failed;

  if (strcmp(p->name, "batch_mean") == 0)
  {
    failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_RC, num->mast_cols * num->batch_id, length, Dj+1) == 0 ? FALSE : TRUE);
  }
  else
  {
    failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_RC, 0, length, Dj+1) == 0 ? FALSE : TRUE);
  }
  Dj[0] = one_norm(Dj+1, length);
  
  if (failed)
	  return FALSE;
  return TRUE;
}


/**********************************************************************************\
 ** 2011.10.30: This function gets the character string name of the objective in 
 ** the argument problem.
 \**********************************************************************************/
int get_objname(one_problem *p, char *buf, int bufspace, int *psurplus)
{ int status = 0;
  char str[] = "Obj";
  /* status = GRBgetstrattrarray((GRBmodel*)p->lp, GRB_STR_ATTR_CONSTRNAME,0, 1,&buf);*/
  /* No way to get objective name in Gurobi, so all objective name will be "Obj"*/
  strcpy(buf, str);
  return status;
}


/**********************************************************************************\
 ** 2011.10.30: This function gets the character string names of the selected rows in 
 ** the argument problem.
 \**********************************************************************************/

#define WORDSIZE 16

int get_rowname(one_problem *p, char **name, char *namestore, int storespace,
                int *psurplus, int begin, int end)
{
  int status;
  status = GRBgetstrattrarray((GRBmodel*)p->lp, GRB_STR_ATTR_CONSTRNAME, begin, end-begin+1, name);
  return status;
}


/**********************************************************************************\
 ** 2011.10.30: This function gets the character string names of the selected columns
 ** in the argument problem.
 \**********************************************************************************/
int get_colname(one_problem *p, char **name, char *namestore, int storespace,
                int *psurplus, int begin, int end)
{
  return GRBgetstrattrarray((GRBmodel*)p->lp, GRB_STR_ATTR_VARNAME, begin, end-begin+1, p->cname);
}


/**********************************************************************************\
 ** 2011.10.30: This function changes the type (LP, QP, etc) of the argument problem.
 \**********************************************************************************/
int change_probtype(one_problem *p, int type)
{
	//not used in gurobi....
	return 0;
  //return CPXchgprobtype (env, p->lp, type);
}


/**********************************************************************************\
 ** 2011.10.30: This function changes one or more coefficients of the objective in 
 ** the argument problem.
 \**********************************************************************************/
int change_objective(one_problem *p, int cnt, int *indices, double *values)
{
  int status;
  status = GRBsetdblattrlist((GRBmodel*)p->lp, GRB_DBL_ATTR_OBJ, cnt, indices, values);
  GRBupdatemodel((GRBmodel*)p->lp);
  return status;
}


/*********************************************************************************\
 ** 2011.10.30: This function changes one or more right hand sides of constraints/rows
 ** in the argument problem.
 \*********************************************************************************/
int change_rhside(one_problem *p, int cnt, int *indices, double *values)
{ 
  int status;
  status = GRBsetdblattrlist((GRBmodel*)p->lp, GRB_DBL_ATTR_RHS, cnt, indices, values);
  GRBupdatemodel((GRBmodel*)p->lp);
  return status;
}


/**********************************************************************************\
 ** 2011.10.30: This function changes one or more bounds on variables/columns in the
 ** argument problem.
 \**********************************************************************************/
int change_bound(one_problem *p, int cnt, int *indices, char *lu, double *bd)
{
	
	int i;
	for (i=0; i < cnt; i++)
	{
		if (lu[i] == 'L')
		{
			GRBsetdblattrelement((GRBmodel*)p->lp, GRB_DBL_ATTR_LB, indices[i], bd[i]);
		}
		else if (lu[i] == 'U')
		{
			GRBsetdblattrelement((GRBmodel*)p->lp, GRB_DBL_ATTR_UB, indices[i], bd[i]);
		}
		else
		{
      GRBsetdblattrelement((GRBmodel*)p->lp, GRB_DBL_ATTR_LB, indices[i], bd[i]);
      GRBsetdblattrelement((GRBmodel*)p->lp, GRB_DBL_ATTR_UB, indices[i], bd[i]);
		}
	}
  
  GRBupdatemodel((GRBmodel*)p->lp);
  
	return 0;
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
  int status;

  status = GRBchgcoeffs((GRBmodel*)p->lp, coef->cnt, coef->row, coef->col, coef->val);

  return (status == 0 ? TRUE : FALSE);
}


/***********************************************************************\
 ** This function will change the coefficients of one column of the
 ** constraint matrix (or the right hand side if column is specified as -1)
 ** The non-zero coefficients of the new column are specified in the
 ** _coef_ vector, and are assumed to be contiguous from _start_ to _stop_.
 ** It returns FALSE if one or more of the coefficients could not be
 ** changed; TRUE otherwise.
 \***********************************************************************/
BOOL change_col(one_problem *p, int column, vector coef,
                int start, int stop)
{
  int		row;
	
  for (row = start; row < stop; row++)
  {
    if (column == -1)
    {
      if (GRBsetdblattrelement((GRBmodel*)p->lp,GRB_DBL_ATTR_RHS, row, coef[row-start]))
        return FALSE;
    }
    else
    {
       if (GRBchgcoeffs((GRBmodel*)p->lp, 1, &row, &column, &coef[row-start]))
         return FALSE;
    }
  } 
  
  GRBupdatemodel((GRBmodel*)p->lp);

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
BOOL change_row(one_problem *p, int row, vector coef,
                int start, int stop)
{
	int		col;

  for (col = start; col < stop; col++)
  {
    if (row == -1)
    {
      if (GRBsetdblattrelement((GRBmodel*)p->lp, GRB_DBL_ATTR_OBJ, col, coef[col-start]))
        return FALSE;
    }
    else
    {
      if (GRBchgcoeffs((GRBmodel*)p->lp, 1, &row, &col, &coef[col-start]))
        return FALSE;
    }
  }
  
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
BOOL add_row(one_problem *p, int start, int stop, int *coef_col,
             double *coef, char sense, double yrhs)
{
  int  status;
  
  static int cumul_num = 0;
  
  char		**r_names;
  double	xrhs[1];
  char		xsense[1] = {'G'};
  int		rmatbeg[2] = {0};
  
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
  rmatbeg[1] = stop-start+1;

  status = GRBaddconstr((GRBmodel*)p->lp, stop-start+1, coef_col, coef, sense, yrhs, r_names[0]);
  
  GRBupdatemodel((GRBmodel*)p->lp);
    
  mem_free(r_names[0]); mem_free(r_names);

  return (status == 0 ? TRUE : FALSE);
}


BOOL remove_row(one_problem *p, int row_num)
{
  int status;

  status = GRBdelconstrs((GRBmodel*)p->lp, 1, &row_num);

  GRBupdatemodel((GRBmodel*)p->lp);

  return (status == 0 ? TRUE : FALSE);  
}

/***********************************************************************\
 \***********************************************************************/
BOOL add_row_to_master(one_problem *p, int start, int stop, int *coef_col,
                       double *coef, char sense, double yrhs)
{
  static int cumul_num = 0;
  int status;
  char		**r_names;
  double	xrhs[1];
  char		xsense[1] = {'G'};
  int		rmatbeg[1] = {0};
  
#ifdef DEBUG
  int idx;
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
  
  status = GRBaddconstr((GRBmodel*)p->lp, stop-start+1, coef_col, coef, sense, yrhs, r_names[0]);
  
  GRBupdatemodel((GRBmodel*)p->lp);
  
  mem_free(r_names[0]); mem_free(r_names);
  
  return (!status);
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
  
  ans = GRBaddconstr((GRBmodel*)p->lp, nzcnt, coef_col, coef, sense, yrhs, r_names[0]);
  
  GRBupdatemodel((GRBmodel*)p->lp);
  
  mem_free(r_names[0]); mem_free(r_names);
  
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
void write_prob (one_problem *p, char *file_name)
{
  int status;
  
  GRBupdatemodel((GRBmodel*)p->lp);
  
  status = GRBwrite( (GRBmodel*)p->lp, file_name);
  
  if ( status ) {
    fprintf (stderr, "Failed to write the problem -- %s.\n", file_name);
    exit (1);
  }
}


/**********************************************************************************\
 ** 2011.10.30: This function sets an integer parameter affecting Solver behavior.
 \**********************************************************************************/
int set_intparam(one_problem *p, const char *whichparam, int newvalue)
{
  return GRBsetintparam(GRBgetenv(p->lp), whichparam, newvalue);
}

/***********************************************************************\
 Get the the number of Q matrix nonzeros that can be read. 
 \***********************************************************************/
int get_qp_nzreadlim()
{
	/*
   int	status = 0;
   int param = 0;
   
   status = CPXgetintparam(env, CPX_PARAM_QPNZREADLIM, &param); 
   printf("CPX_PARAM_QPREADLIM = %d", param);
   
   return param;
   */
	return 0; // not used
}

/***********************************************************************\
 Set the the number of Q matrix nonzeros that can be read. 
 It returns zero on success. zl, 09/25/05. 
 \***********************************************************************/
int set_qp_nzreadlim(int nzreadlim)
{
	/*
   int	status = 0;
   int param = 0;
   
   status = CPXsetintparam(env, CPX_PARAM_QPNZREADLIM, nzreadlim);
   status = CPXgetintparam(env, CPX_PARAM_QPNZREADLIM, &param); 
   printf("CPX_PARAM_QPREADLIM = %d", param);
   
   return status;
   */
	return 0;
}


/**********************************************************************************\
 ** 2011.10.30: This function sets the QP coefficients of the argument problem from
 ** the vector of diagonal elements, for a separable QP.
 \**********************************************************************************/
int copy_qp_separable(one_problem *p, double *qsepvec)
{
    int i, ret = 0;
  GRBdelq((GRBmodel*)p->lp);
	for (i = 0; i< p->mac; i++)
	{
    /* NOTE: CPLEX evaluates the corresponding objective with a factor of 0.5 
     in front of the quadratic objective term.*/
    qsepvec[i] *= 0.5; 
    
		ret += GRBaddqpterms((GRBmodel*)p->lp, 1, &i, &i, &qsepvec[i]);
	}
  GRBupdatemodel((GRBmodel*)p->lp);
	return ret;
}


/***********************************************************************\
 ** This function obtains the vector of lower bounds (lb) for a given problem.  
 ** Note that, like most other vectors, the lb vector reserves the 0th location 
 ** for its 1-norm, and it is assumed that _length_ does not include the 0th 
 ** position. It returns TRUE if the query was successful; FALSE otherwise.
 \***********************************************************************/
BOOL get_lb(vector lb, one_problem *p, int length)
{
  BOOL	failed;

  failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_LB, 0, length-1, lb+1) == 0 ? FALSE: TRUE);
  lb[0] = one_norm(lb+1, length);

  if (failed)
    return FALSE;
  return TRUE;
}

/***********************************************************************\
 ** This function obtains the vector of upper bounds (ub) for a given problem.  
 ** Note that, like most other vectors, the ub vector reserves the 0th location 
 ** for its 1-norm, and it is assumed that _length_ does not include the 0th 
 ** position. It returns TRUE if the query was successful; FALSE otherwise.
 \***********************************************************************/
BOOL get_ub(vector ub, one_problem *p, int length)
{
  BOOL	failed;

  failed = (GRBgetdblattrarray((GRBmodel*)p->lp, GRB_DBL_ATTR_UB, 0, length-1, ub+1) == 0 ? FALSE: TRUE);
  ub[0] = one_norm(ub+1, length);
   
  if (failed)
    return FALSE;
  return TRUE;
}


/**********************************************************************************\
 ** 20112.04.26: This function solves an lp problem with simplex algorithm
 \**********************************************************************************/
int solve_lp(one_problem *p)
{
  GRBupdatemodel((GRBmodel*)p->lp);
  /* Turn off output log 2012.05.02 Yifan */
  set_intparam (p, PARAM_SCRIND, 0);
  
  return GRBoptimize(p->lp);
}


/**********************************************************************************\
 ** 20112.04.26: This function solves an lp problem with simplex algorithm
 \**********************************************************************************/
void *clone_prob(one_problem *p)
{ 
  return GRBcopymodel(p->lp);
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
  /*
   -1 = automatic, 
   0  = primal simplex, 
   1  = dual simplex, 
   2  = barrier, 
   3  = concurrent, 
   4  = deterministic concurrent
   */
  status = set_intparam (p, PARAM_QPMETHOD, 2); /* 2011.10.30 */
  if (status) {
    fprintf (stderr, 
             "Failed to set the optimization method, error %d.\n", status);
    exit(0); 
  }
  
  status = 0;
  
  /*status = set_intparam(p, GRB_INT_PAR_CROSSOVERBASIS, 0);*/
  /*
   Automatic (-1)                         GUROBI decides.
   Disable, use interior solution (0)     No crossover performed.
   Dual first,finish with primal (1)      Pushes dual variables first then primal, finishes with primal.
   Dual first, finsih with dual (2)       Pushes dual variables first then primal, finishes with dual.
   Primal first, finish with primal (3)   Pushes primal variables first then dual, finishes with primal.
   Primal first, finsih with dual (4)     Pushes primal variables first then dual, finishes with dual.
   */
  status = set_intparam (p, PARAM_BARCROSSALG, -1); /* 2011.10.30 */
  
  /* 0: choose initial basis quickly but unstable
     1: take longer time but return a stable basis
   */
  status = set_intparam (p, PARAM_CROSSOVERBASIS, 1); /* 2011.10.30 */
  
  if (status) {
    fprintf (stderr, 
             "Failed to set the barcrossover off, error %d.\n", status);
    exit(0); 
  }
}

/****************************************************************************\
 This function chooses the barrier sparse matrix fill-reducing algorithm.
 When you encouter the "status 12", we should consider a different barrier algorithm.
 Gurobi Paremeter: GRB_INT_PAR_BARORDER
 -1: Automatic (Default)
  0: Approximate Minimum Degree ordering
  1: Nested Dissection ordering.
 \****************************************************************************/
void change_barrier_algorithm(one_problem *p, int k)
{
  int status = 0;
  
#ifdef TRACE
  printf("Inside change_barrier_algorithm.\n");
#endif
  
  status = set_intparam (p, PARAM_BARALG, -1); /* 2011.10.30 */
  
  if (status) {
    fprintf (stderr, 
             "Failed to set the barrier algorithm, error %d.\n", status);
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
  /*
   -1: automatic
    0: primal simplex
    1: dual simplex
    2: barrier
    3: concurrent
    4: deterministic concurrent
   
   */

  
  status = set_intparam (p, PARAM_LPMETHOD, -1); /* 2011.10.30 */
  
  if (status) {
    fprintf (stderr, 
             "Failed to set the optimization method, error %d.\n", status);
    exit(0); 
  }
}

/**********************************************************************************\
 ** 2012.09.10: This function returns the coefficients of the selected position in
 ** the argument problem, in the pointer to a double coefficient.
 \**********************************************************************************/
int get_coef(one_problem *p, int row, int col, double *coef)
{
  return GRBgetcoeff((GRBmodel*)p->lp, row, col, coef);
}


