/*
 * memory.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef MEMORY_H_
#define MEMORY_H_
#include "sdglobal.h"

BOOL duplic_delta_col(sdglobal_type* sd_global, delta_type *delta,
		sigma_type *sigma, num_type *num, int *marked, int omeg_idx, int obs);
BOOL duplicate_omega(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int obs, int *dupl);
void drop_lambda(lambda_type *lambda, delta_type *delta, omega_type *omega,
		sigma_type *sigma, int idx);
void drop_omega(omega_type *omega, delta_type *delta, lambda_type *lambda,
		cut_type *cuts, int drop, int keep);
void drop_sigma(sigma_type *sigma, cut_type *cuts, int idx);
void mark_used_pis(int *marked, cut_type *cuts, int obs);
void thin_data(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
void thin_omega(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
void thin_pi(prob_type *p, cell_type *c, soln_type *s);
void write_histo_file(int *histo, int arr_cnt, int iter, char *fname);
void write_iter(sdglobal_type* sd_global, FILE *fout, prob_type *p,
		cell_type *c, soln_type *s);
void write_statistics(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);

#endif /* MEMORY_H_ */
