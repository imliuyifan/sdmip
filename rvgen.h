/*
 * rvgen.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef RVGEN_H_
#define RVGEN_H_

/*-----------------------------------------------------------------------*/
/*-------------------------rvgen.c---------------------------------------*/
/*-----------------------------------------------------------------------*/

/***********************************************************************\
**  The following inverse normal variate generator was published by
 ** Micheal J. Wichura, University of Chicago in Applied Statistics, as
 ** Algorithm AS 241.  The C function normal() was converted from the
 ** Fortran function PPND7 and produces normal random variates for the
 ** lower tail of a normal distribution accurate to approx. 7 significant
 ** figures.
 \***********************************************************************/

int normal(double *, int, double *, float, float);

/**********************************************************************\
**  The function exponential() generates random variates distributed
 ** exp(1/beta).  The algorithm for this generator is the inverse
 ** transform of the standard exponential distribution and may be
 ** found in Law and Kelton, Simulation Modeling & Analysis, Second
 ** Edition, pg 486.  The algorithm is as follows:
 **     Generate U ~ U(0,1).
 **     Return X = -Beta * ln(U).
 \**********************************************************************/

int exponential(double *, int, double *, float);

int geometric(double *, int, double *, float);

float scalit(float lower, float upper, sd_long *RUN_SEED);

float randUniform(sd_long *SEED);

#endif /* RVGEN_H_ */
