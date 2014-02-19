/*
 * optimal.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef OPTIMAL_H_
#define OPTIMAL_H_
#include "sdglobal.h"

BOOL full_test(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
BOOL optimal(sdglobal_type* sd_global, prob_type *p, cell_type *c, soln_type *s,
		double *phi);
BOOL pre_test_1(sdglobal_type* sd_global, soln_type *s);
BOOL pre_test_2(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
cut_type *choose_cuts(prob_type *p, cell_type *c, soln_type *s);
double cal_temp_lb(prob_type *p, cell_type *c, soln_type *s, cut_type *T);
double calculate_ULm(prob_type *p, cell_type *c, cut_type *T, soln_type *s);
double solve_temp_master(sdglobal_type* sd_global, prob_type *p, cut_type *T,
		cell_type *c);
int randfun(int greatest, sd_long *seed);
void empirical_distrib(omega_type *omega, int *cdf);
void reform_cuts(sdglobal_type* sd_global, sigma_type *sigma, delta_type *delta,
		omega_type *omega, num_type *num, cut_type *T, int *observ, int k);
void sample_omega(sdglobal_type* sd_global, int *cdf, int *observ, int k);

#endif /* OPTIMAL_H_ */
