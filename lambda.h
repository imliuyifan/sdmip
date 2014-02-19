/*
 * lambda.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef LAMBDA_H_
#define LAMBDA_H_
#include "sdglobal.h"

int calc_lambda(sdglobal_type* sd_global, lambda_type *lambda, num_type *num,
		vector Pi, BOOL *new_lamb);
lambda_type *new_lambda(int num_iter, int num_lambda, int num_rv_rows,
		coord_type *coord);
void free_lambda(lambda_type *lambda);
void print_lambda(lambda_type *lambda, num_type *num, int idx);
void write_lambda(FILE *fptr, lambda_type *lambda, num_type *num);

#endif /* LAMBDA_H_ */
