/*
 * master.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef MASTER_H_
#define MASTER_H_
#include "sdglobal.h"

BOOL solve_master(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
BOOL solve_QP_master(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
one_problem *new_master(one_problem *master, cut_type *cuts, int extra_cuts,
		vector x_k);
one_problem *orig_new_master(one_problem *master, cut_type *cuts,
		int extra_cuts);
void change_eta_col(one_problem *p, cut_type *cuts, int k, soln_type *soln,
		num_type *num);
void free_master(one_problem *copy);
void update_rhs(sdglobal_type* sd_global, prob_type *prob, cell_type *cell, soln_type *soln);

#endif /* MASTER_H_ */
