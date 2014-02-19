/*
 * subprob.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef SUBPROB_H_
#define SUBPROB_H_
#include "sdglobal.h"

int solve_subprob(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, vector Xvect, int omeg_idx);
one_problem *new_subprob(one_problem *subprob);
vector compute_rhs(sdglobal_type* sd_global, num_type *num, sparse_vect *Rbar, sparse_matrix *Tbar,
		vector X, omega_type *omega, int omeg_idx);
void free_subprob(one_problem *subprob);

#endif /* SUBPROB_H_ */
