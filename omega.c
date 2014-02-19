/***********************************************************************\
**
 ** omega.c
 **
 **
 **
 ** generate_observ()
 ** equal_obs()
 ** valid_omega_idx()
 ** next_omega_idx()
 ** init_R_T_omega()
 ** get_R_T_omega()
 ** print_omega()
 ** new_omega()
 ** free_omega()
 **
 ** History:
 **   ?? Oct 1991 - <Jason Mai> - created.
 **   16 Dec 1991 - <Jason Mai> - fixed compilation errors.
 **   18 Feb 1992 - <Jason Mai> - updated with new structures.
 **   10 Oct 1992 - <Jason Mai> - added filter, most, last.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 **
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "supomega.h"
#include "omega.h"
#include "log.h"
#include "sdglobal.h"

/***********************************************************************\
** This function obtains a new vector of realizations of the random variables.
 ** It compares the new vector with all previous vectors, looking for
 ** a duplication.  If it finds a duplicate, it returns the index of that
 ** duplicate; otherwise, it adds the vector to the list of distinct
 ** realizations and returns the index of that realization.  It assumes
 ** that the get_omega_idx function returns the 1-norm!
 \***********************************************************************/
int generate_observ(sdglobal_type* sd_global, omega_type *omega, num_type *num,
		BOOL *new_omeg, sd_long *RUN_SEED)
{
	int *observ;
	int cnt, next;

#ifdef TRACE
	printf("Inside generate_observ\n");
#endif

	/*
	 ** Make room for the observations and their 1-norm, then request a 
	 ** new vector of realizations (indices) along with the 1-norm.
	 */
	if (!(observ = arr_alloc(num->cipher+1, int)))
		err_msg("Allocation", "generate_observ", "observ");
	observ[0] = get_omega_idx(sd_global, observ + 1, NULL, 1, RUN_SEED);
	/* Compare vector with all the previous observations */
	for (cnt = 0; cnt < omega->most; cnt++)
		if (valid_omega_idx(omega, cnt))
			if (equal_obs(observ, omega->idx[cnt], num->cipher))
			{
				(*new_omeg) = FALSE;
				omega->weight[cnt]++;
				mem_free(observ);
				//omega->batch_idx[omega->k] = omega->idx[cnt];
				return cnt;
			}

	/* Add the realization vector to the list */
	next = next_omega_idx(omega);
	omega->idx[next] = observ;
	omega->weight[next] = 1;
	omega->filter[next] = USED;
	(*new_omeg) = TRUE;

	/* Update counters */
	if (next >= omega->most)
		omega->most = next + 1;
	omega->cnt++;

#ifdef ENCODE
	if (omega->cnt == 0)
	print_omega(omega, num, 0);
#endif

#ifdef OMEG
	printf("cnt=%d, next=%d, most=%d.\n", omega->cnt, omega->next, omega->most);
#endif

#ifdef TRACE
	printf("Exiting generate_observ\n");
#endif
	//omega->batch_idx[omega->k] = omega->idx[next];
	return next;
}

/***********************************************************************\
 ** This function obtains a new vector of realizations of the random variables.
 ** It adds the vector to the list of realizations and returns the index of
 ** that realization.  It assumes that the get_omega_idx function returns the 1-norm!
 \***********************************************************************/
int get_observ(sdglobal_type* sd_global, omega_type *omega, num_type *num, BOOL *new_omeg)
{
  int next;
  
  next = next_omega_idx(omega);
  omega->weight[next] = 1;
  omega->filter[next] = USED;
  (*new_omeg) = TRUE;
  
  if (next >= omega->most)
    omega->most = next+1;
  
  omega->cnt++;
  
  return next;
}

/***********************************************************************\
** This function determines whether or not two vectors of realizations
 ** are identical (in which case, the contain exactly the same indices into
 ** the list of discrete observations).  It returns TRUE of they are equal;
 ** FALSE otherwise.  Note that if the 1-norms have been stored in the zeroth
 ** location, then they will be checked first.  (It assumes that the
 ** observation vector is one element longer than _len_, for the 1-norm).
 \***********************************************************************/
BOOL equal_obs(int *a, int *b, int len)
{
	int cnt;

#ifdef LOOP
	printf("Inside equal_obs\n");
#endif

	for (cnt = 0; cnt <= len; cnt++)
		if (a[cnt] != b[cnt])
			return FALSE;

	return TRUE;
}

/***********************************************************************\
** This function returns TRUE if the omega->idx array is storing an
 ** realization vector at the index passed to it.  It return FALSE
 ** if the index is out-of-bounds or does not contain a realization.
 \***********************************************************************/
BOOL valid_omega_idx(omega_type *omega, int idx)
{

#ifdef LOOP
	printf("Inside valid_omega_idx()");
#endif
	if (idx < 0 || idx >= omega->most)
		return FALSE;

	/* Check the range of valid indices */
	if (idx < 0 || idx >= omega->most)
		return FALSE;

	/* Check the contents of the idx array */
	if (omega->filter[idx] == UNUSED)
		return FALSE;

	return TRUE;
}

/***********************************************************************\
** This function returns the index of the next available location in
 ** the omega->idx array.  At first, these will just be consecutive
 ** positions in the array.  But, once some omegas get dropped, the
 ** next available position may be in the middle of the array.
 ** It assumes that an available location exists!
 \***********************************************************************/
int next_omega_idx(omega_type *omega)
{
	int cnt;

#ifdef LOOP
	printf("Inside next_omega_idx()");
#endif

	for (cnt = omega->next; omega->filter[cnt] == USED; cnt++)
		; /* Loop until an UNUSED position is found */

	omega->next = cnt;
	return cnt;
}

/***********************************************************************\
** This function sets the references in Romega and Tomega to point
 ** to their corresponding locations in the structure omega.
 ** When the values stored at omega->RT are changed, then the values
 ** of Romega and Tomega will also be changed.  Recall that the
 ** 1-norms are stored in the zeroth location of each array in omega.
 ** Romega and Tomega do NOT actually have 1-norms stored in position 0;
 ** however, their values start at position 1 for consistency, and the
 ** 0th location contains some value which shouldn't be changed.
 **
 ** Romega and Tomega should be moved into the prob data structure.
 ** They will be initialized only once....  no, they can't, because
 ** omega will alloc and free, so pointers will become invalid.
 ** Just keep it this way.
 \***********************************************************************/
void init_R_T_omega(sparse_vect *Romega, sparse_matrix *Tomega,
		omega_type *omega, num_type *num)
{

#ifdef TRACE
	printf("Inside init_R_T_omega\n");
#endif

	/* Assume the R(omega) values come first in all arrays of omega */
	Romega->cnt = num->rv_R;
	Romega->row = omega->row;
	Romega->val = omega->RT;

	/* Start T(omega) arrays where the R(omega) arrays left off */
	Tomega->cnt = num->rv_T;
	Tomega->row = omega->row + num->rv_R;
	Tomega->col = omega->col + num->rv_R;
	Tomega->val = omega->RT + num->rv_R;

#ifdef TRACE
	printf("Exiting init_R_T_omega\n");
#endif
}

/***********************************************************************\
** This function updates the values referenced by Romega and Tomega
 ** by changing the vector of realizations stored at omega->val.
 ** Romega and Tomega are presumed to have pointers into this array,
 ** so by changing it, the function is also changing Romega and Tomega.
 ** It's odd that Romega and Tomega are not required as parameters.
 ** Out of habit, the 1-norm is stored in 0th location.
 ** 
 ** Eventually, it must also decode the observations, which have been
 ** stored in a bitwise manner.
 \***********************************************************************/
void get_R_T_omega(sdglobal_type* sd_global, omega_type *omega, int obs_idx)
{

#ifdef LOOP
	printf("Inside get_R_T_omega");
#endif
#ifdef OMEGA_FILE
  omega->RT[0] = get_omega_vals_from_file(sd_global, obs_idx, omega->RT+1, omega->fidx);
#else
  omega->RT[0] = get_omega_vals(sd_global, omega->idx[obs_idx] + 1, omega->RT + 1);
#endif
}

/***********************************************************************\
**
 ** Must set omega->weight[drop] = 0 ! ! ! ! ! ! ! ! !
 **
 \***********************************************************************/

/***********************************************************************\
** This function displays the indices of a particular observation
 ** of omega.  It also displays the row and column coordinates of
 ** each variable.  It is intended for debugging purposes.
 \***********************************************************************/
void print_omega(omega_type *omega, num_type *num, int idx)
{
	int cnt;

	printf("\nOmega %d:: %d : ", idx, omega->weight[idx]);
	for (cnt = 0; cnt <= num->cipher; cnt++)
		printf("%d ", omega->idx[idx][cnt]);
	printf("\nOmega rows::");
	for (cnt = 0; cnt <= num->rv; cnt++)
		printf("%d ", omega->row[cnt]);
	printf("\nOmega cols::");
	for (cnt = 0; cnt <= num->rv; cnt++)
		printf("%d ", omega->col[cnt]);
	printf("\n");
}

/***********************************************************************\
** This function allocates memory for an omega structure.  It
 ** allocates the structure itself, the row, col, and RT arrays,
 ** and the array of pointers to observation vectors (idx).  It 
 ** allocates the weight and filter arrays, and intializes cnt and next.
 ** However, the actual arrays of indices for each observation are
 ** NOT allocated, since this is done as each realization is observed.
 \***********************************************************************/
omega_type *new_omega(int num_iter, int num_rv, coord_type *coord)
{
	omega_type *omega;

#ifdef TRACE
	printf("Inside new_omega\n");
#endif

	if (!(omega = (omega_type *) mem_malloc (sizeof(omega_type))))
		err_msg("Allocation", "new_omega", "omega");

	if (!(omega->RT = (double *) mem_calloc (num_rv+1, sizeof(double))))
		err_msg("Allocation", "new_omega", "omega->RT");

	if (!(omega->weight = (int *) mem_calloc (num_iter, sizeof(int))))
		err_msg("Allocation", "new_omega", "omega->weight");

	/* Calloc automatically initializes the array to zero */
	if (!(omega->filter = (int *) mem_calloc (num_iter, sizeof(int))))
		err_msg("Allocation", "new_omega", "omega->filter");

	/* Yifan 03/20/2012 Test for omega issues*/
	if (!(omega->used_opt_cut = (BOOL *) mem_calloc (num_iter, sizeof(BOOL))))
		err_msg("Allocation", "new_omega", "omega->filter");

	if (!(omega->idx = (int **) mem_calloc (num_iter, sizeof(int *))))
		err_msg("Allocation", "new_omega", "omega->idx");

	if (!(omega->batch_idx = (int **) mem_calloc (num_iter, sizeof(int *))))
		err_msg("Allocation", "new_omega", "omega->idx");
  
    if(!(omega->fidx = (sd_long *) mem_calloc (num_iter, sizeof(sd_long))))
      err_msg("Allocation", "new_omega", "omega->idx");

	omega->cnt = 0;
	omega->next = 0;
	omega->most = 0;
	omega->col = coord->omega_col;
	omega->row = coord->omega_row;
    omega->fidx[0] = 0;

	return omega;
}

/***********************************************************************\
** This function frees ALL the memory associated with
 ** the omega structure, including the structure itself.  
 ** Note that omega->row and omega->col are not freed, 
 ** since the arrays actually belong to the coord structure.
 \***********************************************************************/
void free_omega(omega_type *omega)
{
	int cnt;

#ifdef TRACE
	printf("Inside free_omega\n");
#endif

	for (cnt = 0; cnt < omega->most; cnt++)
		if (valid_omega_idx(omega, cnt))
			mem_free(omega->idx[cnt]);

	mem_free(omega->idx);
	mem_free(omega->weight);
	mem_free(omega->filter);
	mem_free(omega->RT);
	mem_free(omega->used_opt_cut);
	mem_free(omega->batch_idx);
    mem_free(omega->fidx);
	mem_free(omega);
}
/*  */
