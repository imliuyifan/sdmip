//
//  solver.h
//  sd1.3
//
//  Created by Yifan Liu on 4/26/12.
//  Copyright (c) 2012 The Data Driven Decision Lab. All rights reserved.
//
#ifndef SOLVER_H_
#define SOLVER_H_

#include "sdconstants.h"
#include "sdglobal.h"
#ifdef CPLEX
int set_intparam(one_problem *p, int whichparam, int newvalue);
#else
int set_intparam(one_problem *p, const char *whichparam, int newvalue);
#endif
int change_probtype(one_problem *p, int type);
int copy_qp_separable(one_problem *p, double *qsepvec);
int copy_type(one_problem *p, char *typex);
int change_rhside(one_problem *p, int cnt, int *indices, double *values);
int change_bound(one_problem *p, int cnt, int *indices, char *lu, double *bd);
int get_basis(one_problem *p, int *cstat, int *rstat);
int get_x(one_problem *p, double * x, int begin, int end);
int get_numrows(one_problem *p);
int get_numcols(one_problem *p);
int get_numnz(one_problem *p);
int get_rows(one_problem *p, int *pnzcnt, int *rmatbeg, int *rmatind,
		double *rmatval, int rmatspace, int *psurplus, int begin, int end);
int get_cols(one_problem *p, int *pnzcnt, int *cmatbeg, int *cmatind,
		double *cmatval, int cmatspace, int *psurplus, int begin, int end);
int get_coef(one_problem *p, int row, int col, double *coef);
void *read_problem(one_problem *p, char *filename, char *filetype);
int get_obj(one_problem *p, double *obj, int begin, int end);
int get_rhs(one_problem *p, double *rhsx, int begin, int end);
int get_sense(one_problem *p, char *senx, int begin, int end);
int get_type(one_problem *p, char *typex, int begin, int end);
int get_lbound(one_problem *p, double *lb, int begin, int end);
int get_ubound(one_problem *p, double *ub, int begin, int end);
int get_objname(one_problem *p, char *buf, int bufspace, int *psurplus);
int get_rowname(one_problem *p, char **name, char *namestore, int storespace,
		int *psurplus, int begin, int end);
int get_colname(one_problem *p, char **name, char *namestore, int storespace,
		int *psurplus, int begin, int end);
int change_objective(one_problem *p, int cnt, int *indices, double *values);
int solve_lp(one_problem *p);
void *clone_prob(one_problem *p);
void change_solver_barrier(one_problem *p);
void change_barrier_algorithm(one_problem *p, int k);
void change_solver_primal(one_problem *p);
BOOL setup_problem(one_problem *current);
BOOL print_problem(one_problem *p, char *filename);
BOOL solve_problem(sdglobal_type* sd_global, one_problem *p);
double get_objective(one_problem *p);
BOOL get_primal(vector X, one_problem *p, int length);
BOOL get_dual(vector Pi, one_problem *p, num_type *num, int length);
BOOL get_dual_slacks(vector Dj, one_problem *p, num_type *num, int length);
void remove_problem(one_problem *p);
BOOL change_col(one_problem *p, int column, vector coef, int start, int stop);
BOOL change_row(one_problem *p, int row, vector coef, int start, int stop);
BOOL add_row(one_problem *p, int start, int stop, int *coef_col, double *coef,
		char sense, double rhs);
BOOL add_row_to_master(one_problem *p, int start, int stop, int *coef_col,
		double *coef, char sense, double yrhs);
BOOL add_row_to_batch(one_problem *p, int start, int nzcnt, int *coef_col,
		double *coef, char sense, double yrhs, int batch_id);
BOOL remove_row(one_problem *p, int row_num);
void write_prob(one_problem *p, char *file_name);
int get_qp_nzreadlim(void);
int set_qp_nzreadlim(int nzreadlim);
void close_Solver(void);
void open_Solver(void);
BOOL change_coef(one_problem *p, sparse_matrix *coef);
BOOL get_lb(vector lb, one_problem *p, int length);
BOOL get_ub(vector ub, one_problem *p, int length);

#endif

