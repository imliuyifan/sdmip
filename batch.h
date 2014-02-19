/*
 * batch.h
 *
 *  Created on: Jun 11, 2013
 *      Author: lian
 */

#ifndef BATCH_H_
#define BATCH_H_

#include "sdglobal.h"

batch_incumb_type * new_batch_incumb(sdglobal_type* sd_global, prob_type *p, vector x_k);
BOOL get_beta_x(sdglobal_type* sd_global, soln_type *s, vector Beta,
		one_problem *p, num_type *num, int length);
one_problem * new_batch_problem(one_problem * master, int max_cuts);
void add_batch_equality(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
void add_cut_to_batch(sdglobal_type* sd_global, one_cut *cut, prob_type *p,
		cell_type *c, soln_type *s, int batch_id, int cut_position);
void add_fcut_to_batch(sdglobal_type* sd_global, one_cut *cut, prob_type *p,
		cell_type *c, soln_type *s, int batch_id, int iteration_num);
void save_batch_incumb(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s, int batch_id);
void update_batch_bounds(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int batch_id);
void update_batch_rhs(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int batch_id);

#endif /* BATCH_H_ */
