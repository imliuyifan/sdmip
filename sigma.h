/*
 * sigma.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef SIGMA_H_
#define SIGMA_H_
#include "sdglobal.h"

int calc_sigma(sdglobal_type* sd_global, cell_type *c, sigma_type *sigma,
		num_type *num, vector pi_k, sparse_vect *Rbar, sparse_matrix *Tbar,
		int lamb_idx, BOOL new_lamb, BOOL *new_sigma);
sigma_type *new_sigma(int num_iter, int num_nz_cols, int num_pi,
		coord_type *coord);
void free_sigma(sigma_type *sigma);
void print_sigma(sigma_type *sigma, num_type *num, int idx);
void write_sigma(FILE *fptr, sigma_type *sigma, num_type *num);

#endif /* SIGMA_H_ */
