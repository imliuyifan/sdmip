/*************************************************************************\
**
 ** log.c
 **
 **  This file contains functions to replace the normal calls to 
 **  alloc, malloc, calloc, realloc, and free with another set of functions,
 **  my_alloc, my_malloc, my_calloc, my_realloc, and my_free which
 **  store information about all dynamically allocated memory.  
 **
 **  The new functions maintain a list of all addresses which currently
 **  have memory allocated to them.  If any memory has not been freed by
 **  the end of the program, its address will still be in the list. 
 **
 **  To use the functions best, define macros like mem_alloc, mem_free, etc
 **  which immediately call my_alloc, my_free, etc when you want to debug,
 **  but call the regular alloc, free, etc. when you want to ship.
 **
 ** History: 
 **   20 Aug 92 - Jason Mai - created.
 **   Should NOT be used without the consent of either Suvrajeet Sen or
 **   Jason Mai
 ** 
 \*************************************************************************/

#include "prob.h"
#include "log.h"
#include "sdglobal.h"
#include <time.h>

sd_long MEM_USED = 0;

/*************************************************************************\
** This function records the information present in a call to one of the
 ** dynamic allocation routines.  It stores a string, intended to be the
 ** name of the pointer or the size of the allocation, and a pointer,
 ** intended to be the return value from the memory allocation.
 \*************************************************************************/
void *log_alloc(char *string, void *return_ptr, int size)
{
	/*
	 fprintf(sd_global->fptrALLOC, "%lld : %s : size %d\n", return_ptr, string, size);
	 */

	MEM_USED += size;

	return return_ptr;
}

/*************************************************************************\
** This function records information present in a call to free some
 ** dynamically allocated memory.  It records a string (presumed to be
 ** the variable name of the pointer passed to free()) and the value of 
 ** the pointer.  It returns its second pointer parameter, assumed to be
 ** the result of a free operation
 \*************************************************************************/

void *log_free(char *string, void *ptr)
{
	//fprintf(fptrFREE, "%d : %s\n", ptr, string); 
	//modified by Yifan to avoid wrong type warning July 26 2011

	free(ptr);
	return NULL;
}

/*************************************************************************\
** A call to realloc involves both a free (of the original pointer) and
 ** an alloc (for the new pointer).  This function accepts a string,
 ** an original pointer for the free log file, and the same string with 
 ** the realloc'ed pointer for the alloc log file.
 \*************************************************************************/
void *log_realloc(char *string, void *free_ptr, void *alloc_ptr, int size)
{
	/*
	 fprintf(fptrFREE, "%lld : %s\n", free_ptr, string);
	 fprintf(sd_global->fptrALLOC, "%lld : %s : size %d\n", alloc_ptr, string, size);
	 */

	MEM_USED += size;

	return alloc_ptr;
}

/*************************************************************************\
** This function opens the dynamic allocation log files, in preparation
 ** for the the logging of allocation events.  This includes both the
 ** alloc and the free log files.  It assumes the global variables
 ** fptrALLOC and fptrFREE have been defined somewhere.
 \*************************************************************************/
void log_start(sdglobal_type* sd_global)
{
#ifdef WRITE
	FILE *fptr;
#endif
	 /*sd_global->fptrALLOC = fopen("alloc.log", "w");
	 sd_global->fptrFREE  = fopen("free.log", "w");*/

	MEM_USED = 0;
	sd_global->LAST_CLOCK = clock();

#ifdef WRITE
	/* Erase any previous data files */
	fptr = fopen(ITER_DAT, "w");
	fclose(fptr);

	fptr = fopen(CNT_DAT, "w");
	fclose(fptr);
#endif

}

/*************************************************************************\
** This function closes the dynamic allocation log files, including both
 ** the alloc and the free log files.  It assumes the global variables
 ** fptrALLOC and fptrFREE have been defined and initialized somewhere.
 \*************************************************************************/
void log_stop(sdglobal_type* sd_global)
{
	/* FILE		*fptr;*/
#ifdef WRITE
	if (sd_global->fptrALLOC != NULL)
	{
		fclose(sd_global->fptrALLOC);
	}
	if (sd_global->fptrFREE != NULL)
	{
		fclose(sd_global->fptrFREE);
	}
#endif
}

/* modified by zl 
 void malloc_verify()
 {  }
 */

BOOL malloc_verify(void)
{
	return TRUE;
}
