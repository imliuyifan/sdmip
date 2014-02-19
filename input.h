/*
 * input.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef INPUT_H_
#define INPUT_H_
#include "sdglobal.h"

int copy_names(char **x_name, char ***xname, char **xstore, unsigned *xstorsz,
		int max);
/*--------------------------------------------------------------------------*/
/*--------------------------------input.c-----------------------------------*/
/*--------------------------------------------------------------------------*/

int get_file(FILE **, char *, char *, char *);

/**********************************************************************\
**   The function load_core calls the function get_line() which returns
 ** the data obtained from the input file, separated into fields, data-
 ** type, and number of fields.  The function then determines the nature
 ** of the data and stores the data in the struct original of type
 ** one_problem.
 \**********************************************************************/

int load_core(one_problem **, identity **, char *, int);
/*  Added by zl. */
int load_core_cpx(one_problem **, identity **, char *, int);
void allocate_arrays(one_problem **, identity **);
void getnamespace(one_problem *, int *, int *, int, int); /* 2011.10.30 */

/**********************************************************************\
**   The function load_stoch calls the function get_line() which returns
 ** the data obtained from the input file, separated into fields, data-
 ** type, and number of fields.  The function then determines the nature
 ** of the data and stores the data in the struct omega of type omega.
 **
 ** for DISCRETE values, the fields are as follows:
 ** 	field1 = ' '
 ** 	field2 = column name
 ** 	field3 = row name
 ** 	field4 = value
 ** 	field6 = prob of value
 ** Options other than discrete are: normal, uniform, exponential,
 ** gamma and geometric.  All distributions share the same field values
 ** for fields 1-3 but use field 4 and 6 differently.  The uses are as
 ** follows:
 ** 	NORMAL		field4 = mean			field6 = var
 **	UNIFORM		field4 = lowerbound		field6 = upperbound
 **	EXPONENTIAL	field4 = mean
 **	GEOMETRIC	field4 = mean
 **
 \**********************************************************************/

int load_stoch(sdglobal_type* sd_global, one_problem *, identity *, char *);

/**********************************************************************\
**   The function load_time receives pointers to smalls and a pointer
 ** to a identity structure.  The function sets the values of the
 ** smalls to the indices of the row and column where the main
 ** problem is to be split.
 \**********************************************************************/

int load_time(sd_small *, sd_small *, identity *, char *);
int load_config(sdglobal_type* sd_global, BOOL read_seeds, BOOL read_iters);

/**********************************************************************\
** These functions free the memory required for reading the files
 ** and providing observations of omega.
 \**********************************************************************/
void free_ident(identity *);
void free_omegas(sdglobal_type* sd_global);

/**********************************************************************\
**  The function str_to_float converts a string passed to it as an
 ** argument and returns a double.  The string may contain an exponent
 ** part as well as a real part.
 \**********************************************************************/

double str_to_float(char *);

/************************************************************************\
**  The function get_mean() receives an array of probabilities and an
 ** array of corresponding values for a discrete distribution along
 ** with an integer representing the number of values contained in each
 ** array.  The function then calculates and returns the mean of the
 ** distribution.
 \************************************************************************/

double get_mean(double *, double *, int);

/*************************************************************************\
**  The function uniform receives an integer representing the number of
 ** discrete elements to be generated, an empty array of doubles for
 ** storage of the discrete elements and two doubles representing the
 ** upper and lower bounds on the distribution.  The function then
 ** calculates the discrete elements on the distribution corresponding
 ** to the probabilities and loads the storage array.
 \*************************************************************************/

void uniform(int, double *, double, double);

#endif /* INPUT_H_ */
