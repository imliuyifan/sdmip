/*
 * improve.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef IMPROVE_H_
#define IMPROVE_H_
#include "sdglobal.h"

BOOL check_improvement(sdglobal_type* sd_global, prob_type *prob,
		cell_type *cell, soln_type *soln, double *phi);
BOOL new_check_improvement(sdglobal_type* sd_global, prob_type *prob,
		cell_type *cell, soln_type *soln, double *phi);
double max_cut_height(sdglobal_type* sd_global, cut_type *cuts, vector X, cell_type *c, num_type *num);
void new_incumbent(prob_type *p, cell_type *c, soln_type *s, double est);

#endif /* IMPROVE_H_ */
