/*
 * cuts.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef CUTS_H_
#define CUTS_H_
#include "sdglobal.h"

void SD_cut(sdglobal_type* sd_global, sigma_type *sigma, delta_type *delta,
		omega_type *omega, num_type *num, one_cut *cut, vector Xvect,
		double *pi_ratio, double max_ratio, double min_ratio, int num_samples,
		BOOL *dual_statble_flag);
i_type compute_istar(int obs, one_cut *cut, sigma_type *sigma,
		delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
		double *argmax, BOOL pi_eval, int ictr);
i_type compute_new_istar(int obs, one_cut *cut, sigma_type *sigma,
		delta_type *delta, vector Xvect, num_type *num, vector Pi_Tbar_X,
		double *argmax, int ictr);
void free_cut(one_cut *cut);
one_cut *new_fea_cut(int num_x, int num_istar, int num_samples);
int FEA_cut(sdglobal_type* sd_global, cell_type *cell, soln_type *soln,
		sigma_type *sigma, delta_type *delta, omega_type *omega, num_type *num,
		int num_samples, BOOL *dual_statble_flag, BOOL new_omega,
		BOOL new_sigma);
int add_to_cutpool(sdglobal_type* sd_global, double *alpha, double *beta,
		cell_type *c, soln_type *s, int mast_cols);
int FEA_cut_check_add(sdglobal_type* sd_global, cell_type *cell,
		prob_type *prob, soln_type *soln, vector x_k);
void decrease_feacut_rownum(cell_type *c);
void update_dual_size(cell_type *c, soln_type *s, prob_type *p);
cut_type *new_cuts(int num_cuts, int num_x, int num_betas);
one_cut *new_cut(int num_x, int num_istar, int num_samples);
void free_cuts(cut_type *cuts);
void print_cut(cut_type *cuts, num_type *num, int idx);
void print_cut_info(cell_type *c, num_type *num, char *phrase);
BOOL stochastic_updates(sdglobal_type* sd_global, cell_type *c,
		lambda_type *lambda, sigma_type *sigma, delta_type *delta,
		omega_type *omega, num_type *num, sparse_vect *Rbar,
		sparse_matrix *Tbar, one_problem *subprob, vector Pi, int omeg_idx,
		BOOL new_omega);
void form_new_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int omeg_idx, BOOL new_omega);
void remove_all_cuts(sdglobal_type* sd_global, prob_type *p, cell_type *c,
                     soln_type *s);
void form_ip_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
                 soln_type *s, int omeg_idx, BOOL new_omega, int mip_idx);
void form_fea_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int omeg_idx, vector x_k, BOOL new_omega);
BOOL form_incumb_cut(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s, int omeg_idx);
void reduce_cuts(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
void thin_cuts(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		soln_type *s);
void drop_cut(int cnt, prob_type *p, cell_type *c, soln_type *s);
int add_cut(sdglobal_type* sd_global, one_cut *cut, prob_type *p, cell_type *c,
		soln_type *s);
void add_cut_to_master(sdglobal_type* sd_global, one_cut *cut, prob_type *p,
		cell_type *c, soln_type *s, int idx);
batch_cut_type *new_bcuts(prob_type *p, int bsize, batch_cut_type *batch_cuts);
void free_bcuts(batch_cut_type *bcuts);

#endif /* CUTS_H_ */
