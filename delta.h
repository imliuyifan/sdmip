/*
 * delta.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef DELTA_H_
#define DELTA_H_
#include "sdglobal.h"

delta_type *new_delta(int num_iter, coord_type *coord);
void calc_delta_col(sdglobal_type* sd_global, delta_type *delta, lambda_type *lambda, omega_type *omega,
		num_type *num, int obs);
void calc_delta_row(sdglobal_type* sd_global, delta_type *delta,
		lambda_type *lambda, omega_type *omega, num_type *num, int pi_idx);
void drop_delta_col(delta_type *delta, lambda_type *lambda, int col);
void drop_delta_row(delta_type *delta, lambda_type *lambda, omega_type *omega,
		int row);
void free_delta(delta_type *delta, omega_type *omega, int num_lamb);
void print_delta(delta_type *delta, num_type *num, int idx, int obs);

#endif /* DELTA_H_ */
