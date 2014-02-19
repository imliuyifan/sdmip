/*
 * omega.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef OMEGA_H_
#define OMEGA_H_
#include "sdglobal.h"

BOOL equal_obs(int *a, int *b, int len);
BOOL valid_omega_idx(omega_type *omega, int idx);
int generate_observ(sdglobal_type* sd_global, omega_type *omega, num_type *num,
		BOOL *new_omeg, sd_long *RUN_SEED);
int get_observ(sdglobal_type* sd_global, omega_type *omega, num_type *num, BOOL *new_omeg);
int next_omega_idx(omega_type *omega);
omega_type *new_omega(int num_iter, int num_rv, coord_type *coord);
void free_omega(omega_type *omega);
void get_R_T_omega(sdglobal_type* sd_global, omega_type *omega, int obs_idx);
void init_R_T_omega(sparse_vect *Romega, sparse_matrix *Tomega,
		omega_type *omega, num_type *num);
void print_omega(omega_type *omega, num_type *num, int idx);

#endif /* OMEGA_H_ */
