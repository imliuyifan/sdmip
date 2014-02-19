/**********************************************************************\
**	The random variate generators in this file all require similar
 **	function call with the first three parameters being:

 !!!!!!!!!!!!!! Note... I have changed the prob_array from float to double
 !!!!!!!!!!!!!! throughout this file, so it matches other modules...  I have
 !!!!!!!!!!!!!! tried to make sure this doesn't violate anything...
 !!!!!!!!!!!!!!    Jason. 25 March 92.

 **	1)	double *prob_array	array of percentiles for random variates
 **	2)	int num_probs		number of random variates to be ganerated
 **	3)	double *rv_array	array to hold random variates corresponding
 **							to the percentiles in prob_array
 **	Each function will take one or two additional parameters depending
 **	on the nature of the random variates to be generated.  
 **	example: normal  two additional parameters; mean and standard dev.
 **	int normal(float *prob_array, int num_probs, double *rv_array, 
 **			float mean, float sigma)
 \**********************************************************************/

#include "prob.h"
#include "rvgen.h"

/* What are these global variables for? */
/* They are replaced by local variables in every function. */
/* I'm commenting them out... (Jason)  22 Mar 92 */

/***********************************************************************\
**	The following inverse normal variate generator was published by 
 **	Micheal J. Wichura, University of Chicago in Applied Statistics, as 
 **	Algorithm AS 241.  The C function normal() was converted from the 
 **	Fortran function PPND7 and produces normal random variates for the 
 **	lower tail of a normal distribution accurate to approx. 7 significant 
 **	figures. 
 \***********************************************************************/

int normal(double *prob_array, int num_probs, double *rv_array, float mu,
		float sigma)
{
	int i;
	float zero, one, half, split1, split2, const1, const2, a0, a1, a2, a3, b1;
	float b2, b3, c0, c1, c2, c3, d1, d2, e0, e1, e2, e3, f1, f2, p, q, r;
	float endval;

	for (i = 0; i < num_probs; i++)
	{
		p = (float) prob_array[i];

		zero = 0.0;
		one = 1.0;
		half = one / 2.0;
		split1 = 0.425;
		split2 = 5.0;
		const1 = 0.180625;
		const2 = 1.6;

		/* coefficients for p close to 1/2 */
		a0 = 3.3871327179;
		a1 = 50.434271938;
		a2 = 159.29113202;
		a3 = 59.109374720;
		b1 = 17.895169469;
		b2 = 78.775757664;
		b3 = 67.18756360;

		/* coefficients for p neither close to 1/2 nor 0 or 1 */
		c0 = 1.4234372777;
		c1 = 2.7568153900;
		c2 = 1.3067284816;
		c3 = .17023821103;
		d1 = .73700164250;
		d2 = .12021132975;

		/* coefficients for p near 0 or 1 */
		e0 = 6.6579051150;
		e1 = 3.0812263860;
		e2 = .42868294337;
		e3 = .017337203997;
		f1 = .24197894225;
		f2 = .012258202635;

		q = p - half;

		if (fabs(q) <= split1)
		{
			r = const1 - q * q;
			endval = q * (((a3 * r + a2) * r + a1) * r + a0)
					/ (((b3 * r + b2) * r + b1) * r + one);
			rv_array[i] = mu + sigma * endval;
			continue;
		}

		if (q < 0.0)
		{
			r = p;
		}
		else
		{
			r = one - p;
		}

		if (r <= zero)
		{
			return 0;
		}

		r = sqrt(-log(r));

		if (r <= split2)
		{
			r = r - const2;
			endval = (((c3 * r + c2) * r + c1) * r + c0)
					/ ((d2 * r + d1) * r + one);
			rv_array[i] = endval;
		}
		else
		{
			r = r - split2;
			endval = (((e3 * r + e2) * r + e1) * r + e0)
					/ ((f2 * r + f1) * r + one);
			rv_array[i] = endval;
		}
		if (q < 0)
		{
			rv_array[i] = -1 * rv_array[i];
		}
		rv_array[i] = mu + sigma * rv_array[i];
	}

	return (1);
}

/**********************************************************************\
**	the function exponential() generates random variates distributed
 **	exp(1/beta).  The algorithm for this generator is the inverse 
 **	transform of the standard exponential distribution and may be 
 **	found in Law and Kelton, Simulation Modeling & Analysis, Second
 **	Edition, pg 486.  The algorithm is as follows:
 **		Generate U ~ U(0,1).
 **		Return X = -Beta * ln(U).
 \**********************************************************************/

int exponential(double *prob_array, int num_probs, double *rv_array, float beta)
{
	int i;

	for (i = 0; i < num_probs; i++)
	{
		rv_array[i] = (-1 * beta) * log(prob_array[i]);
	}

	return (1);
}

/* ------------------------------------------------------------------*/
/* ------------------------------------------------------------------*/

int geometric(double *prob_array, int num_probs, double *rv_array, float lambda)
{
	int i;

	for (i = 0; i < num_probs; i++)
	{
		rv_array[i] = (-1 / lambda) * log(prob_array[i]);
	}

	return (1);
}

float scalit(float lower, float upper, sd_long *RUN_SEED)
{
	float val, wide;

	wide = upper - lower;

	val = randUniform(RUN_SEED);

	/* Yifan 06/25/2012 batch mean */
	//val = randUniform(&config.RUN_SEED1);
	return ((wide * val) + lower);
}

/* ------------------------------------------------------------------*/
/* ------------------------------------------------------------------*/

float randUniform(sd_long *SEED)
{
	/* static int to static long int: modified by Yifan 2013.02.18 */
	static int lo_bits, hi_bits;

	lo_bits = ((*SEED) & 0xFFFFL) * 16807;
	hi_bits = (int) (((*SEED) >> 16) * 16807) + (lo_bits >> 16);
	*SEED = ((lo_bits & 0xFFFFL) - 0x7FFFFFFFL) + ((hi_bits & 0x7FFFL) << 16)
			+ (hi_bits >> 15);
	return ((*SEED) < 0 ? ((*SEED) += 0x7FFFFFFFL) : (*SEED)) * 4.656612875E-10;
}

