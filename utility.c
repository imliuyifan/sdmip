/***********************************************************************\
**
 ** utility.c
 **
 ** This file contains miscellaneous functions needed for coding the 
 ** SD algorithm.
 **
 ** one_norm()
 ** equal_arr()
 ** expand_vect()
 ** reduce_vect()
 ** print_vect()
 ** print_sparse_vect()
 ** print_sparse_matrix()
 ** print_num()
 ** PIxR()
 ** PIxT()
 ** TxX()
 ** compute_Mu()
 ** History:
 **   ?? Oct 1991 - <Jason Mai> - created.
 **   14 Dec 1991 - <Jason Mai> - prepared for testing / compilation.
 **   25 Feb 1992 - <Jason Mai> - changed structures.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include "prob.h"
#include "utility.h"
#include "solver.h"
#include "sdconstants.h"
#include "log.h"
#include "sdglobal.h"

/*
 ** This function calculates the 1-norm (here, the sum of the absolute
 ** value of all elements in the array) of any array of doubles passed to
 ** it.  It returns the 1-norm.
 */
double one_norm(double *a, int len)
{
	int cnt;
	double sum;

	sum = 0.0;
	for (cnt = 0; cnt < len; cnt++)
	{
		sum += DBL_ABS(a[cnt]);
	}

	return sum;
}

/*
 ** This function compares two arrays of doubles (assumed to be the same
 ** length) element by element to determine whether or not the are equal
 ** up to a certain tolerance.  It returns TRUE if the arrays are
 ** equal; FALSE otherwise.  Note that if the zeroth element is
 ** the 1-norm of the vector, it will be checked first.  (It assumes
 ** that the vector is one element larger than _len_, to accomodate
 ** the 1-norm).
 ** 
 ** The tolerance is relative to the size of the numbers... so it 
 ** represents the percentage of the 1-norm up to which j
 */
BOOL equal_arr(double *a, double *b, int len, double tolerance)
{
	int cnt;

	for (cnt = 0; cnt <= len; cnt++)
		if (DBL_ABS(a[cnt] - b[cnt]) > DBL_ABS(tolerance * a[cnt]))
			return FALSE;

	return TRUE;
}

/*
 ** Are functions passing in len or len+1 ???
 */

/*
 ** This function duplicates the array passed to it by copying it
 ** into a second array (assumed of to be the right length) passed to it.
 ** The 1-norms are assumed to take up the zeroth position in the
 ** first array, and are not counted in the _len_ parameter.
 void copy_arr(double *a, double *b, int len)
 {
 int 		i;

 for (i = 0; i <= len; i++)
 b[i] = a[i];
 }
 */

/*
 ** This function, like copy_arr, duplicates an array passed to it.
 ** However, it allocates space for the new copy, fills it, and returns
 ** a pointer to the copy, rather than just filling it.  As usual,
 ** the _len_ parameter does not include the 1-norm, assumed in position 0.
 */
double *duplic_arr(double *a, int len)
{
	int i;
	double *b;

	if ((b = arr_alloc(len+1, double)))
	{
		for (i = 1; i <= len; i++)
			b[i] = a[i];
		b[0] = one_norm(b + 1, len);
	}
	else
		err_msg("Allocation", "duplic_arr", "b");

	return b;
}

/*
 ** This function expands a sparse vector into a full vector.  Each
 ** (non-zero) element in the sparse vector is placed in its appropriate
 ** place in the full vector, and the rest of the positions are filled
 ** with zeroes.  The function returns a pointer to a dynamically 
 ** allocated array containing the full vector.  
 */
vector expand_vect(double *s_vect, int *row, int num_elem, int length)
{
	int cnt;
	vector new_vect;

	if (!(new_vect = arr_alloc(length+1, double)))
		err_msg("Allocation", "expand_vect", "new_vect");

	for (cnt = 0; cnt <= length; cnt++)
		new_vect[cnt] = 0;

	for (cnt = 1; cnt <= num_elem; cnt++)
		new_vect[row[cnt]] = s_vect[cnt];
	new_vect[0] = one_norm(new_vect + 1, length); /* = s_vect[0]; */

	return new_vect;
}

/*
 ** This function compresses a full vector into a sparse vector, 
 ** according to the array of row coordinates passed to it.  It pulls
 ** out only the values stored at the locations referenced by row. 
 ** It assumes that the final length of the vector is known ahead 
 ** of time! (this is always the case for SD).  It returns a pointer 
 ** to a dynamically allocated array containing the sparse vector.
 ** Note: the zeroth element still contains the 1-norm!
 */
double *reduce_vect(double *f_vect, int *row, int num_elem)
{
	int cnt;
	double *s_vect;

	if (!(s_vect = arr_alloc(num_elem+1, double)))
		err_msg("Allocation", "reduce_vect", "s_vect");

	for (cnt = 1; cnt <= num_elem; cnt++)
		s_vect[cnt] = f_vect[row[cnt]];
	s_vect[0] = one_norm(s_vect + 1, num_elem);

	return s_vect;
}

/* CALLER: DONT FORGET TO FREE THE ORIGINAL FULL VECTOR !! !! !! !! */

void print_num(sdglobal_type* sd_global, num_type *num)
{
	printf("mast_rows:%d\t\tmast_cols:%d\n", num->mast_rows, num->mast_cols);
	printf("sub_rows:%d\t\tsub_cols:%d\t\trv:%d\n", num->sub_rows,
			num->sub_cols, num->rv);
	printf("rv_rows:%d\t\trv_cols:%d\t\trv_R:%d\t\trv_T:%d\n", num->rv_rows,
			num->rv_cols, num->rv_R, num->rv_T);
	printf("nz_cols:%d\t\tmax_cuts:%d\t\tmin_iter:%d\n", num->nz_cols,
			num->max_cuts, sd_global->config.MIN_ITER);
	printf("MIN_ITER: %d\n\n", sd_global->config.MIN_ITER);
	printf(
			"Notation:\n'+':new incumbent solution found.\n'<':in-sample test fails.\n'>':in-sample test succeeds.\n");
}

/*
 ** Assumes 1-norm is stored at position 0
 */
void print_vect(vector X, int size, char *string)
{
	fprint_vect(stdout, X, size, string);
}

/* This one prints to a file / stream, not just to the screen */
void fprint_vect(FILE *fptr, vector X, int size, char *string)
{
	int cnt;
	fprintf(fptr, "Vector %s :: ", string);
	for (cnt = 0; cnt <= size; cnt++)
		fprintf(fptr, "%.20f ", X[cnt]);
	fprintf(fptr, "\n");
}

/*
 ** Assumes 1-norm is stored at position 0
 */
void print_sparse_vect(sparse_vect *V, char *string)
{
	int cnt;

	printf("\nSparse Vect %s (%d)::", string, V->cnt);
	for (cnt = 0; cnt <= V->cnt; cnt++)
		printf("%f ", V->val[cnt]);
	printf("\nSparse Vect Rows :: ");
	for (cnt = 0; cnt <= V->cnt; cnt++)
		printf("%d ", V->row[cnt]);
	printf("\n");
}

/*
 ** Assumes 1-norm is stored at position 0
 */
void print_sparse_matrix(sparse_matrix *V, char *string)
{
	int cnt;

	printf("\nSparse matrix %s (%d)::", string, V->cnt);
	for (cnt = 0; cnt <= V->cnt; cnt++)
		printf("%f ", V->val[cnt]);
	printf("\nSparse matrix Rows :: ");
	for (cnt = 0; cnt <= V->cnt; cnt++)
		printf("%d ", V->row[cnt]);
	printf("\nSparse matrix cols :: ");
	for (cnt = 0; cnt <= V->cnt; cnt++)
		printf("%d ", V->col[cnt]);
	printf("\n");
}

/* 
 ** This function multiplies a full vector (PI) times a sparse vector
 ** (R) to produce a scalar.  Thus, the first vector is assumed to be
 ** a row, and the second a column.  _num_ specifies the number of 
 ** elements in the second, sparse vector, while _row_ specifies the
 ** row number of the non-zero elements of the sparse vector.
 ** It is assumed that pi_k and R have 1-norms occupying position 0.
 */
double PIxR(vector pi_k, sparse_vect *R)
{
	int cnt;
	double pi_R;

	pi_R = 0.0;
	for (cnt = 1; cnt <= R->cnt; cnt++)
		pi_R += R->val[cnt] * pi_k[R->row[cnt]];

	return pi_R;
}

/*
 ** This function multiplies a full column vector (PI) times a sparse
 ** matrix (T) to produce a full vector.  _row_ specifies the row number
 ** of each (non-zero) element in the sparse matrix, while _num_ specifies
 ** the number of these elements.  It returns a pointer to a dynamically
 ** allocated array containing the answer (its size is specified with the 
 ** _length_ parameter, equal to the number of columns in PI).  Note that 
 ** the zeroth element of the array is set to the 1-norm of the vector. 
 ** It is also assumed that both pi_k and T have 1-norms at location 0.
 */
vector PIxT(vector pi_k, sparse_matrix *T, int length)
{
	int cnt;
	vector pi_T;

	if (!(pi_T = arr_alloc(length+1, double)))
		err_msg("Allocation", "PIxT", "pi_T");
    
    
//#pragma omp parallel for private(cnt, a) num_threads(2)
//    for (cnt = 0; cnt <= length; cnt++){
//        a[cnt] = 0;
//    }
    
//#pragma omp parallel for private(cnt, temp)
//#pragma omp parallel for private(cnt,temp) shared(pi_T) num_threads(1)
    for (cnt = 0; cnt <= length; cnt++){
        pi_T[cnt] = 0;}

//#pragma omp parallel for private(cnt)
	for (cnt = 1; cnt <= T->cnt; cnt++){
		pi_T[T->col[cnt]] += pi_k[T->row[cnt]] * T->val[cnt];
    }

	pi_T[0] = one_norm(pi_T + 1, length);

	return pi_T;
}

/*
 ** This function mulitplies a sparse matrix times a full column
 ** vector, to produce a full column vector.  The result of each
 ** multiplication is subtracted from the _ans_ vector in the 
 ** appropriate location (it is assumed to have already been allocated 
 ** and initialized with desired values).  It assumes the 0th element
 ** is reserved for the 1-norm, in the answer, the T sparse vector,
 ** and the X vector.
 */
vector TxX(sparse_matrix *T, vector X, vector ans)
{
	int cnt;

	for (cnt = 1; cnt <= T->cnt; cnt++)
		ans[T->row[cnt]] -= T->val[cnt] * X[T->col[cnt]];

	return ans;
}

/* This function will return the positive value of A * X */
vector TxX_plus(sparse_matrix *T, vector X, vector ans)
{
	int cnt;

	for (cnt = 1; cnt <= T->cnt; cnt++)
		ans[T->row[cnt]] += T->val[cnt] * X[T->col[cnt]];

	return ans;
}

/*
 ** This function multiplies two vectors to produce a scalar.
 ** (thus, it assumes a row then a column vector are passed)
 ** The 0th element, as usual, is reserved for the 1-norm; however,
 ** the routine does NOT fill this location.
 */
double CxX(vector c, vector x, int len)
{
	int cnt;
	double sum;

	sum = 0.0;
	for (cnt = 1; cnt <= len; cnt++)
		sum += c[cnt] * x[cnt];

	return sum;
}

BOOL encode(int *plain, one_key *key, int *cipher, int len)
{
	int i;

	for (i = 0; i < len; i++)
		cipher[key[i].element] |= plain[i] << key[i].shift;

	return TRUE;
}

BOOL decode(int *cipher, one_key *key, int *plain, int len)
{
	int i;

	for (i = 0; i < len; i++)
		plain[i] = (cipher[key[i].element] >> key[i].shift) & key[i].mask;

	return TRUE;
}

/*
 ** This function initializes a key to be used when encoding/decoding
 ** arrays whose (integer) elements have a discrete & finite range,
 ** specified within the array parameter _ranges_.  The created key 
 ** (which may be used by encode() and decode()) is assumed to be of
 ** length _num_ranges_.  The return value specifies the number of 
 ** integers which will be required for the cipher array.
 */
int form_key(one_key *key, int *ranges, int num_ranges)
{
	int cnt, elem; /* counters for the loops */
	int num_bits; /* # of bits needed to store a given range */
	int num_cipher; /* # of ints being used in cipher array */
	int *used_bits; /* # bits left in each int of cipher array */
	BOOL not_done;

	/* Note: There can't be more cipher ints than ranges to be encoded. */
	used_bits = arr_alloc(num_ranges, int); /* initialized to zero */

	/* 
	 ** For each range (corresponding to a position in the plain array)
	 ** initialize a key structure describing how to encode/decode it.
	 */

#ifdef CAL_CHECK
	printf("**********Start forming keys**********\n");
#endif

	for (cnt = 0; cnt < num_ranges; cnt++)
	{
		/* Find out how many bits this range is going to need */
		num_bits = get_num_bits(ranges[cnt]);

#ifdef CAL_CHECK
		printf("cnt = %d\n",cnt);
		printf("num_bits = %d\n",num_bits);
#endif

		/*
		 ** Find an int in the cipher array which has room for
		 ** the number of bits required for this range.
		 */
		not_done = TRUE;
		for (elem = 0; not_done && elem < num_ranges; elem++)
		{
#ifdef CAL_CHECK
			printf("elem = %d;\n",elem);
			printf("used_bits[%d] = %d\n",elem,used_bits[elem]);
#endif
			if (used_bits[elem] + num_bits <= MAX_BITS)
			{
				key[cnt].element = elem; /* range is in this int */
				key[cnt].shift = used_bits[elem]; /* range is displaced in int */
				key[cnt].mask = (1 << num_bits) - 1; /* sub-pattern within int */
				used_bits[elem] += num_bits; /* do accounting for int */
				not_done = FALSE;
#ifdef CAL_CHECK
				printf("key[%d].element = %d;\tkey[%d].shift = %d;\tkey[%d].mask = %d;\tused_bits[%d]=%d;\n",
						cnt,key[cnt].element,cnt,key[cnt].shift,cnt,key[cnt].mask,elem,used_bits[elem]);
#endif
			}
#ifdef CAL_CHECK
			printf("----------------\n");
#endif
		}
	}

	for (num_cipher = 0; used_bits[num_cipher] > 0; num_cipher++)
		; /* count the number of ints used to encode all the ranges */

#ifdef CAL_CHECK
	for (cnt=0; cnt<num_ranges; cnt++)
	{
		printf("used_bits[%d] = %d\n",cnt,used_bits[cnt]);
	}
	printf("num_cipher = %d;\n",num_cipher);
	printf("**********Stop forming keys**********\n");
#endif

	mem_free(used_bits);
	return num_cipher;
}

/*
 ** This function returns the minimum number of bits needed to represent
 ** a given number, passed as the parameter _num_.
 */
int get_num_bits(int num)
{
	int hi_bit = 1;
	int num_bits;

	for (num_bits = 0; hi_bit <= num; num_bits++)
		hi_bit = hi_bit << 1;
#ifdef CAL_CHECK
	printf("num_bits: %d\n",num_bits);
#endif
	return num_bits;
}

/*
 ** This function converts an integer into a character string and copies it
 ** into _string_ starting at location _beg_.   _max_ specifies the
 ** greatest integer expected for the _num_ parameter, and thus determines
 ** the number of characters needed to store the integer.
 */
void filename_number(char *string, int beg, int max, int num)
{
	int cnt;
	int place;

	cnt = 0;
	for (place = max; place >= 1; place /= 10)
	{
		string[beg + cnt] = '0' + (num) / place % 10;
		cnt++;
	}

	printf("\n\nfilename_number::%s.\n\n", string);
}

/***********************************************************************\
 ** This function compute the reduced cost of every second stage 
 ** variables. They will be used to calculate the Mu x R and then added  
 ** to the Pi x R.
 \***********************************************************************/
/*added by Yifan to update _PixR_*/
double compute_Mu(one_problem *p, int sub_cols)
{
	int i;
	double Mu_R;

	vector dj = NULL;
	dj = (double *) malloc((sub_cols + 1) * sizeof(double));
	if (dj == NULL)
	{
		fprintf(stderr, "No memory for solution dj.\n");
		return 1;
	}

	get_dual_slacks(dj, p, NULL, sub_cols);

	Mu_R = MuxR(p, sub_cols, dj);

	if (0)
	{
		for (i = 0; i <= sub_cols; i++)
		{
			printf("***This is reduced cost for y[%d]: %+f\n", i, dj[i]);
		}
	}mem_free(dj);
	return Mu_R;

}

/***********************************************************************\
 ** This function obtains the basis infomation of variables for example basis 
 ** or nonbasis at lowerbound or nonbasis at higher bound. Then return   
 ** the value of Mu x R.
 \***********************************************************************/
/*added by Yifan to update _PixR_*/
double MuxR(one_problem *p, int sub_cols, vector dj)
{
	int j;
	int *cstat = NULL;
	double *y = NULL;
	double MuR = 0;
	char *basismsg;
	cstat = (int *) malloc((sub_cols + 1) * sizeof(int));
	y = (double *) malloc((sub_cols + 1) * sizeof(double));

	get_basis(p, cstat + 1, NULL); /* 2011.10.30 */
	get_x(p, y + 1, 0, sub_cols - 1); /* 2011.10.30 */

	if (0)
	{
		/* Write out the solution */

		for (j = 1; j <= sub_cols; j++)
		{
			printf("y[%d] = %17.10g", j, y[j]);
			if (cstat != NULL)
			{
				switch (cstat[j])
				{
				case AT_LOWER:
					basismsg = "Nonbasic at lower bound";
					break;
				case BASIC:
					basismsg = "Basic";
					break;
				case AT_UPPER:
					basismsg = "Nonbasic at upper bound";
					break;
				case FREE_SUPER:
					basismsg = "Superbasic, or free variable at zero";
					break;
				default:
					basismsg = "Bad basis status";
					break;
				}
				printf("  %s   %d", basismsg, cstat[j]);
			}
			printf("\n");
		}

	}
    
	/*added by Yifan to enable parallel computation*/
//#pragma omp parallel for private(j) shared(MuR)
	for (j = 1; j <= sub_cols; j++)
	{
		if (cstat != NULL)
		{
			switch (cstat[j])
			{
			case AT_LOWER:
                //#pragma omp atomic
				MuR += dj[j] * y[j];
				break;
			case AT_UPPER:
                //#pragma omp atomic
				MuR += dj[j] * y[j];
				break;
			default:
				break;
			}
		}
	}

	mem_free(y);
	mem_free(cstat);

	return MuR;
}

/*
 ** This function calculate the variance of the 
 ** vector x.
 */
double calc_var(sdglobal_type* sd_global, double *x, double *mean_value,
		double *stdev_value, int batch_size)
{
	double mean, vari, temp;
	int count, length;
	double stdev;
	stdev = 10000000.0;
	temp = 0.0;
	mean = x[0];
	vari = 0.0;

	if (mean_value != NULL)
	{
		length = batch_size;
	}
	else
	{
		length = sd_global->config.SCAN_LEN;
	}

	for (count = 1; count < length; count++)
	{
		temp = mean;
		mean = mean + (x[count] - mean) / (double) (count + 1);
		vari = (1 - 1 / (double) count) * vari
				+ (count + 1) * (mean - temp) * (mean - temp);
	}

	if (mean_value != NULL)
	{
		*mean_value = mean;
	}
	if (stdev_value != NULL)
	{
		stdev = sqrt(vari / (double) count);
		*stdev_value = stdev;
	}

	return vari;

}

void calc_mean_stdev(vector *x, vector mean_value, vector stdev_value,
		int num_element, int batch_size)
{
	double ans = 0.0, temp = 0.0;
	double mean, vari = 0.0, stdev = 10000000;
	int count, i;

	for (i = 0; i <= num_element; i++)
	{
		mean = x[0][i];
		vari = 0.0;
		for (count = 1; count < batch_size; count++)
		{
			ans = x[count][i];
			temp = mean;
			mean = mean + (ans - mean) / (double) (count + 1);
			vari = (1 - 1 / (double) count) * vari
					+ (count + 1) * (mean - temp) * (mean - temp);
		}
		stdev = sqrt(vari);

		mean_value[i] = mean;
		stdev_value[i] = stdev;
	}

}

/* function calculate the supnorm of the input array */
double sup_norm(vector a, int size)
{
  int i;
  double max, temp;
  
  /* Here we assume that the zero-th element is saved for one-norm of the vector  */
  max = DBL_ABS(a[1]);
  
  for (i = 2; i <= size; i++) {
    temp = DBL_ABS(a[i]);
    if (temp > max)
      max = temp;
  }
  
  return max;
}

/* This function calculate the relative difference of two array, and the resulting array
 is stored in the last input array*/
void rdiff(sdglobal_type *sd_global, double *a, double *b, int size, vector rd)
{
  int i;
  for (i = 0; i <= size ; i++) {
    if (b[i]<sd_global->config.TOLERANCE && a[i]<sd_global->config.TOLERANCE) {
      rd[i] = .1 * sd_global->config.TOLERANCE;
    }
    else{rd[i] = (a[i]-b[i])/b[i];}
  }
}

