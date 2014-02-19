/*
 * utility.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef UTILITY_H_
#define UTILITY_H_
#include "sdglobal.h"

BOOL decode(int *cipher, one_key *key, int *plain, int len);
BOOL encode(int *plain, one_key *key, int *cipher, int len);
BOOL equal_arr(double *a, double *b, int len, double tolerance);
double *duplic_arr(double *a, int len);
double *reduce_vect(double *f_vect, int *row, int num_elem);
double calc_var(sdglobal_type* sd_global, double *x, double *mean_value,
		double *stdev_value, int batch_size);
double compute_Mu(one_problem *p, int sub_cols);
double CxX(vector c, vector x, int len);
double MuxR(one_problem *p, int sub_cols, vector dj);
double one_norm(double *a, int len);
double PIxR(vector pi_k, sparse_vect *R);
int form_key(one_key *key, int *ranges, int num_ranges);
int get_num_bits(int num);
vector expand_vect(double *s_vect, int *row, int num_elem, int length);
vector PIxT(vector pi_k, sparse_matrix *T, int length);
vector TxX_plus(sparse_matrix *T, vector X, vector ans);
vector TxX(sparse_matrix *T, vector X, vector ans);
void calc_mean_stdev(vector *x, vector mean_value, vector stdev_value,
		int num_element, int batch_size);
void filename_number(char *string, int beg, int max, int num);
void fprint_vect(FILE *fptr, vector X, int size, char *string);
void print_num(sdglobal_type* sd_global, num_type *num);
void print_sparse_matrix(sparse_matrix *V, char *string);
void print_sparse_vect(sparse_vect *V, char *string);
void print_vect(vector X, int size, char *string);
double sup_norm(vector a, int size);
void rdiff(sdglobal_type *sd_global, double *a, double *b, int size, vector rd);

#endif /* UTILITY_H_ */
