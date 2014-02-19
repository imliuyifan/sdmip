/*
 * testout.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef TESTOUT_H_
#define TESTOUT_H_
#include "sdglobal.h"

void cplex_err_msg(sdglobal_type* sd_global, char *string, prob_type *p, cell_type *c, soln_type *s);
void print_contents(one_problem *, char *);

#endif /* TESTOUT_H_ */
