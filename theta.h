/*
 * theta.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef THETA_H_
#define THETA_H_

#include "prob.h"
#include "cell.h"
#include "sdglobal.h"

double c_k(sdglobal_type* sd_global, vector Cost, vector X, cell_type *c,
		num_type *num, int idx);
double cut_height(sdglobal_type* sd_global, one_cut *cut, vector X,
		cell_type *c, num_type *num);
double f_k(sdglobal_type* sd_global, vector x, vector cost, cut_type *cuts,
		theta_type *theta, num_type *num, cell_type *cur_cell, int *best);
double get_max_Vi(sdglobal_type* sd_global, vector x, cut_type *cuts,
		num_type *num, int *this_cut, int last_cut, int *best, cell_type *c);
theta_type *new_theta(int num_cells);
void free_theta(theta_type *theta);

#endif /* THETA_H_ */
