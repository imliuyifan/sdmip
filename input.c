/**********************************************************************\
**	This file contains the functions for interpreting the data contained
 **	int the core, stoch and time files for stochastic decomposition
 **	linear programming problems.
 \**********************************************************************/

#include <string.h>
#include <ctype.h>
#include "prob.h"
#include "solver.h"
#include "supomega.h"
#include "parser.h"
#include "rvgen.h"
#include "input.h"
#include "log.h"
#include "sdconstants.h"
#include "sdglobal.h"

/**********************************************************************\
**	the function load_core calls the function get_line() which returns
 **	the data obtained from the input file, separated into fields, data-
 **	type, and number of fields.  The function then determines the nature
 **	of the data and stores the data in the struct original of type one_
 **	problem.
 \**********************************************************************/

int load_core(one_problem **original, identity **ident, char *fname, int objsen)
{
	FILE *core;
	int i, j; /* default increment counter */
	int result; /* return value from function get_line */
	int wordSize = WORDSIZE; /* default name length for row and col. arrays */
	int blockSize = BLOCKSIZE; /* base allocation for storage arrays */
//	int strSize = STRSIZE; /* base allocation length for name arrays */
	int category; /* category # for ROW, COLUMNS, BOUNDS and RHS */
	int num_blks;
	int rowblock1, colblock1, colblock3;
	int colblock4;
	int objcoef = 0;
	int field_idx;
	char field1[15], field2[60], field3[15];
	char field4[15], field5[15], field6[15], field7[15];
	/* char	fileType[6]; */
	char fieldType, dataType, dummy[1];
	char *nextval[4];
	char *next_fld[4];
	char *fld_ptr[8];
	char name[15];

	/* initialize arrays of strings */
	nextval[0] = field3;
	nextval[1] = field5;
	nextval[2] = field7;
	next_fld[0] = field2;
	next_fld[1] = field4;
	next_fld[2] = field6;
	fld_ptr[0] = field1;
	fld_ptr[1] = field2;
	fld_ptr[2] = field3;
	fld_ptr[3] = field4;
	fld_ptr[4] = field5;
	fld_ptr[5] = field6;
	fld_ptr[6] = field7;

	nextval[3] = dummy;
	next_fld[3] = dummy;
	fld_ptr[7] = dummy;

	dummy[0] = '\0';

	/* allocate structure one_problem original */
	if (!(*original = (one_problem *) mem_malloc(sizeof(one_problem))))
	{
		printf("memory allocation failure\n");
		return 0;
	}

	/* allocate structure identity ident */
	if (!(*ident = (identity *) mem_malloc(sizeof(identity))))
	{
		printf("memory allocation failure\n");
		return 0;
	}

	dummy[0] = '\0'; /* string of length zero */

	result = get_file(&core, fname, ".cor", name);
	if (result == 0)
		return 0;

	/* allocate and load problem name obtained from get_file */
	if (!((*original)->name = (char *) mem_calloc(15, sizeof(char))))
	{
		printf("memory allocation error\n");
		return 0;
	}

	strcpy((*original)->name, name);
	(*original)->objsen = objsen;

	/* Initialize all parameters that won't be used */
	(*original)->mae = 0;
	(*original)->maesz = 0;
	(*original)->enzsz = 0;
	(*original)->estorsz = 0;
	(*original)->rngcol = NULL;
	(*original)->nrowind = NULL;
	(*original)->etype = NULL;
	(*original)->ename = NULL;
	(*original)->estore = NULL;
	(*original)->enzbeg = NULL;
	(*original)->enzcnt = NULL;
	(*original)->enzind = NULL;
	(*original)->enzval = NULL;
	(*original)->dataname = NULL;
	(*original)->rhsname = NULL;
	(*original)->rngname = NULL;
	(*original)->bndname = NULL;

	while (get_line(&core, field1, field2, field3, field4, field5, field6,
			field7, &fieldType))
	{
		if (fieldType == 't')
		{
			num_blks = 1;
			dataType = field1[2]; /* get third char. in category name */

#ifdef LOAD_CORE
			printf("after get_line dataType = %c\n", dataType);
#endif

			/* determine current category and allocate appropriate arrays */
			switch (dataType)
			{
			case 'W': /* category name ROW */

				/* allocate and initialize relevent row variables for cplex */
				(*original)->mar = 0;
				rowblock1 = blockSize;
//				rowblock2 = strSize;

				/* allocate row name array for index reference */
				if (!((*ident)->row_name =
						(char**) mem_calloc(blockSize, sizeof (char *))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				for (i = 0; i < blockSize; i++)
				{
					if (!((*ident)->row_name[i] =
							(char *) mem_calloc(wordSize, sizeof (char))))
					{
						printf("memory allocation error\n");
						return 0;
					}
				}

				/* allocate sense array */
				if (!((*original)->senx =
						(char *) mem_calloc(blockSize, sizeof (char))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				(*original)->marsz = 0;

				/* allocate array for objective name */
				if (!((*original)->objname =
						(char *) mem_calloc(wordSize, sizeof (char))))
				{
					printf("memory allocation error\n");
					return 0;
				}

#ifdef LOAD_CORE
				printf("case_W ori->mar=%d, ori->marsz=%d, rowblock1=blockSize=%d, rowblock2=strSize=%d\n", (*original)->mar, (*original)->marsz, rowblock1, rowblock2);
#endif
				category = 1;
				break;

			case 'L': /* category name COLUMNS */

				copy_names((*ident)->row_name, &(*original)->rname,
						&(*original)->rstore, &(*original)->rstorsz,
						(*original)->mar);

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mar; i++)
				printf("case_L id-rname=%s, ori->rname=%s\n", (*ident)->row_name[i],
						(*original)->rname[i]);
#endif

				/* initialize relevent column arrays and variables for cplex */
				(*original)->mac = 0;
				(*original)->matsz = 0;
				(*original)->macsz = 0;

				colblock1 = blockSize;
//				colblock2 = strSize;
				colblock3 = blockSize;
				colblock4 = 3 * blockSize;

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mar; i++)
				printf("case_L2 id-rname=%s, ori->rname=%s, ori->senx=%c\n",
						(*ident)->row_name[i], (*original)->rname[i], (*original)->senx[i]);

				printf("name = %s, mac = %d, mar = %d, objsen = %d\n", (*original)->name,
						(*original)->mac, (*original)->mar, (*original)->objsen);
				printf("macsz = %d, marsz = %d, matsz = %d, maesz = %d, enzsz = %d\n",
						(*original)->macsz, (*original)->marsz, (*original)->matsz,
						(*original)->maesz, (*original)->enzsz);
#endif

				/* initialize objective coefficient array for cplex */
				if (!((*original)->objx =
						(double *) mem_calloc(colblock4, sizeof (double))))
				{
					printf("memory allocation error\n");
					return 0;
				}

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mar; i++)
				printf("case_L3 id-rname=%s, ori->rname=%s, ori->senx=%c\n",
						(*ident)->row_name[i], (*original)->rname[i], (*original)->senx[i]);

				printf("name = %s, mac = %d, mar = %d, objsen = %d\n", (*original)->name,
						(*original)->mac, (*original)->mar, (*original)->objsen);
				printf("macsz = %d, marsz = %d, matsz = %d, maesz = %d, enzsz = %d\n",
						(*original)->macsz, (*original)->marsz, (*original)->matsz,
						(*original)->maesz, (*original)->enzsz);
#endif

				if (!((*original)->bdl =
						(double *) mem_calloc(colblock4, sizeof (double))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				if (!((*original)->bdu =
						(double *) mem_calloc(colblock4, sizeof (double))))
				{
					printf("memory allocation error\n");
					return 0;
				}

				/* Initialize objective coefficients and upper/lower bounds */
				for (i = 0; i < colblock4; i++)
				{
					(*original)->objx[i] = 0.0;
					(*original)->bdl[i] = 0.0;
					/* Change from cplex70 to cplex81. zl 04/13/05. */
					(*original)->bdu[i] = INFBOUND;
				}
				/* initialize constraint arrays for cplex */
				if (!((*original)->matbeg =
						(int *) mem_calloc(blockSize, sizeof (int))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				if (!((*original)->matcnt =
						(int *) mem_calloc(blockSize, sizeof (int))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				if (!((*original)->matind =
						(int *) mem_calloc(blockSize, sizeof (int))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				if (!((*original)->matval =
						(double *) mem_calloc(blockSize, sizeof (double))))
				{
					printf("memory allocation error\n");
					return 0;
				}

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mar; i++)
				printf("case_L4 id-rname=%s, ori->rname=%s\n", (*ident)->row_name[i],
						(*original)->rname[i]);
#endif

				/* allocate column name array for index reference */
				if (!((*ident)->col_name =
						(char **) mem_calloc(blockSize, sizeof (char *))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				for (i = 0; i < blockSize; i++)
				{
					if (!((*ident)->col_name[i] =
							(char *) mem_calloc(wordSize, sizeof (char))))
					{
						printf("memory allocation error\n");
						return 0;
					}
				}

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mar; i++)
				printf("case_L5 id-rname=%s, ori->rname=%s\n", (*ident)->row_name[i],
						(*original)->rname[i]);
#endif

				category = 2;
				break;

			case 'U': /* category name BOUNDS */

				/* allocate appropriate arrays for constraint bounds */
				/* Don't allocate here, since it was done back with objx arrays                    */
				/*
				 if(!((*original)->bdl =	
				 (double *)mem_calloc((*original)->macsz, sizeof (double))))
				 {
				 printf("memory allocation error\n");
				 return 0;
				 }
				 if(!((*original)->bdu =
				 (double *)mem_calloc((*original)->macsz, sizeof (double))))
				 {
				 printf("memory allocation error\n");
				 return 0;
				 }
				 */
				/* Just set the category to bounds, and exit */
				category = 3;
				break;

			case 'S': /* category name RHS */
				copy_names((*ident)->col_name, &(*original)->cname,
						&(*original)->cstore, &(*original)->cstorsz,
						(*original)->mac);

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mac; i++)
				printf("case_S ori->cname[%d]=%s\n", i, (*original)->cname[i]);
#endif

				/* allocate appropriate arrays for rhs arrays */
				if (!((*original)->rhsx =
						(double *) mem_calloc((*original)->marsz, sizeof (double))))
				{
					printf("memory allocation error\n");
					return 0;
				}

				category = 4;
				break;

			case 'D': /* category name ENDATA */
				/* end of file obtained, return 1 to calling function */
				/* YOU CAN'T JUST RETURN HERE!!! YOU'VE GOT STUFF TO DO */
				/* Wait for get_line to fail ???  Do it here ??? */
				/* Changed from row_name and col_name to num_rows and num_cols..                   Think this is right */
				(*ident)->num_rows = (*original)->mar;
				(*ident)->num_cols = (*original)->mac;

#ifdef LOAD_CORE
				printf("zl_case_D id->num_rows = %d, id->num_cols = %d\n", (*ident)->num_rows,
						(*ident)->num_cols);
#endif

				fclose(core);
				return 1;

			default: /* problem title encountered */
				break;
			} /* End of switch. */
		} /* End of fieldType = 't'. */
		else if (fieldType == 'f')
		{
			switch (category)
			{
			case 1:

				/* load sense of row 'N', 'E', 'L', 'G' */
				if (field1[0] == 'N')
				{
					strcpy((*original)->objname, field2);

#ifdef LOAD_CORE
					printf("case1_N ori->objname = %s\n", (*original)->objname);
#endif
				}
				else
				{
					/* Verify array lengths not exceeded. */
#ifdef LOAD_CORE
					printf("case1 ori->mar = %d\n", (*original)->mar);
#endif

					if ((*original)->mar >= rowblock1) /* senx array exceeded */
					{
						num_blks = (*original)->mar / blockSize;
						++num_blks;
						rowblock1 = num_blks * blockSize;
						if (!((*original)->senx =
								(char *) mem_realloc((*original)->senx,
										rowblock1 * sizeof(char))))
						{
							printf("memory allocation error\n");
							return 0;
						}
						if (!((*ident)->row_name =
								(char **) mem_realloc((*ident)->row_name, rowblock1 * sizeof(char *))))
						{
							printf("memory allocation error\n");
							return 0;
						}
						for (i = (*original)->mar; i < rowblock1; i++)
						{
							if (!((*ident)->row_name[i] =
									(char *) mem_calloc(wordSize, sizeof(char))))
							{
								printf("memory allocation error\n");
								return 0;
							}
						}
					}

					/* Load sense of row */
					(*original)->senx[(*original)->marsz++] = field1[0];

#ifdef LOAD_CORE
					printf("case1_2 ori->senx[%d] = %c\n", (*original)->marsz-1,
							(*original)->senx[(*original)->marsz-1]);
#endif

					/* Load row name into row_name and rstore. */
					strcpy((*ident)->row_name[(*original)->mar++], field2);

#ifdef LOAD_CORE
					printf("case1_3 id->row_name[%d] = %s\n", (*original)->mar-1,
							(*ident)->row_name[(*original)->mar-1]);
#endif
				}

				break;

			case 2:

#ifdef LOAD_CORE
				for (i = 0; i < (*original)->mar; i++)
				printf("case2 id-rname=%s, ori->rname=%s\n", (*ident)->row_name[i],
						(*original)->rname[i]);
#endif

				/* Compare column name for new column */
				if ((*original)->mac == 0
						|| strcmp(field1,
								(*ident)->col_name[(*original)->mac - 1]) != 0) /* new column */
				{
					/* Verify array lengths not exceeded. */
					if ((*original)->mac >= colblock1) /* cname array exceeded */
					{
						num_blks = (*original)->mac / blockSize;
						++num_blks;
						colblock1 = num_blks * blockSize;

						if (!((*original)->matbeg =
								(int *) mem_realloc((*original)->matbeg, colblock1 * sizeof(int))))
						{
							printf("memory allocation error\n");
							return 0;
						}

						if (!((*original)->matcnt =
								(int *) mem_realloc((*original)->matcnt, colblock1 * sizeof(int))))
						{
							printf("memory allocation error\n");
							return 0;
						}

						if (!((*ident)->col_name =
								(char **) mem_realloc((*ident)->col_name, colblock1 * sizeof(char *))))
						{
							printf("memory allocation error\n");
							return 0;
						}

						for (i = (*original)->mac; i < colblock1; i++)
						{
							if (!((*ident)->col_name[i] =
									(char *) mem_calloc(wordSize, sizeof (char))))
							{
								printf("memory allocation error\n");
								return 0;
							}
						}
					}

					if ((*original)->mac >= colblock4) /* objx array exceeded */
					{
						num_blks = (*original)->mac / blockSize;
						++num_blks;
						colblock4 = num_blks * blockSize;

						if (!((*original)->objx =
								(double *) mem_realloc((*original)->objx, colblock4 * sizeof (double))))
						{
							printf("memory allocation error\n");
							return 0;
						}
						if (!((*original)->bdl =
								(double *) mem_realloc((*original)->bdl, colblock4 * sizeof (double))))
						{
							printf("memory allocation error\n");
							return 0;
						}
						if (!((*original)->bdu =
								(double *) mem_realloc((*original)->bdu, colblock4 * sizeof (double))))
						{
							printf("memory allocation error\n");
							return 0;
						}

						for (i = (*original)->mac; i < colblock4; i++)
						{
							(*original)->objx[i] = 0.0;
							(*original)->bdl[i] = 0.0;
							/* Change from cplex70 to cplex81. zl 04/13/05. */
							(*original)->bdu[i] = INFBOUND;
						}
					}

					/* update cname, matbeg */
					(*original)->matbeg[(*original)->mac] = (*original)->matsz;
					/*   ??????  */
					(*original)->matcnt[(*original)->mac] = 0;
					(*original)->macsz++;

#ifdef LOAD_CORE
					for (i = 0; i < (*original)->mar; i++)
					printf("case2_2 ori->rname=%s\n", (*original)->rname[i]);
#endif

					/* Load new column name into col_name and cstore. */
					strcpy((*ident)->col_name[(*original)->mac], field1);
					(*original)->mac++;

#ifdef LOAD_CORE
					printf("case2_3 ori->cname[%d] = %s\n", (*original)->mac-1,
							(*ident)->col_name[(*original)->mac-1]);
#endif
				}
				field_idx = 0;

				/* Identify and load row indices and coefficient values. */

				while (strlen(next_fld[field_idx]) > 0)
				/* repeat while additional coeffients remain */
				{
					/* Verify array lengths not exceeded. */
					if ((*original)->matsz >= (colblock3 - 3))
					{
						/* matind and matval array lengths exceeded */
						num_blks = (*original)->matsz / blockSize;
						num_blks += 2;
						colblock3 = num_blks * blockSize;
						if (!((*original)->matind =
								(int *) mem_realloc((*original)->matind, colblock3 * sizeof(int))))
						{
							printf("memory allocation error\n");
							return 0;
						}
						if (!((*original)->matval =
								(double *) mem_realloc((*original)->matval, colblock3 * sizeof(double))))
						{
							printf("memory allocation error\n");
							return 0;
						}
					}

#ifdef LOAD_CORE
					for (i=0; i <(*original)->mar; i++)
					printf("case2_4 ori-rname[%d] = %s\n", i, (*original)->rname[i]);
#endif

					/* Identify row index number. */
					i = 0;
					while (i < (*original)->mar
							&& strcmp((*original)->rname[i],
									next_fld[field_idx]) != 0)
					{

#ifdef LOAD_CORE
						printf("case2_5 ori->rname[%d]=%s, next_fld[%d]=%s\n", i,
								(*original)->rname[i], field_idx, next_fld[field_idx]);
#endif

						i++;
					}
					objcoef = 0;
					if (i >= (*original)->mar)
					{

#ifdef LOAD_CORE
						printf("case2_6 next_fld[%d]=%s, ori->objname=%s\n", field_idx,
								next_fld[field_idx], (*original)->objname);
#endif

						if (!(strcmp(next_fld[field_idx], (*original)->objname)))
							objcoef = 1;
						else
						{
							printf(
									"undefined row name in constraint matrix: %s\n",
									next_fld[field_idx]);
							printf(
									"Randomness in Cost Coefficients and W metrix will be supported in later release\n");
							/*added by Yifan Oct 12 2011*/
							return (0);
						}
					}
					/* Load row indices and constraint coefficients. */
					if (objcoef == 0)
					{
						(*original)->matind[(*original)->matsz] = i;
						(*original)->matval[(*original)->matsz] = str_to_float(
								nextval[field_idx]);
						++((*original)->matcnt[(*original)->mac - 1]);
						(*original)->matsz++;
					}
					else
					{
						(*original)->objx[(*original)->mac - 1] = str_to_float(
								nextval[field_idx]);
					}
					++field_idx;
				}
				break;

			case 3:

				/* identify bound coordinates and load bounds */
				/* Temporarily disconnected, since no input files have bounds... */
				/* Use default bounds of Xi >= 0 */
				break;

			case 4:

				/* get row number by name */
				i = 0;
				while (i < (*original)->mar
						&& strcmp((*ident)->row_name[i], field2) != 0)
				{
					i++;
				}

				if (i >= (*original)->mar)
				{
					printf("undefined row name in constraint matrix\n");
					printf(
							"Randomness in Cost Coefficients and W metrix will be supported in later release\n");
					/*added by Yifan Oct 12 2011*/
					return (0);
				}
				/* convert rhs value from string to float */
				(*original)->rhsx[i] = str_to_float(field3);
				break;

			default:
				break;
			} /* End of switch. */

			/* reset all fields to length zero */
			for (i = 0; i < 6; i++)
				for (j = 0; j < 15; j++)
					fld_ptr[i][j] = '\0';
		} /* End of fieldType = 'f'. */
	} /* End of while (get_line()). */

	/* This won't ever happen !!!! unless get_line fails !!!! */

	/* Changed from row_name and col_name to num_rows and num_cols... 
	 Think this is right */
	(*ident)->num_rows = (*original)->mar;
	(*ident)->num_cols = (*original)->mac;

	fclose(core);
	printf("finished loading core file\n");
	return (1);
}

/**********************************************************************\
**	the function load_stoch calls the function get_line() which returns
 **	the data obtained from the input file, separated into fields, data-
 **	type, and number of fields.  The function then determines the nature
 **	of the data and stores the data in the struct omega of type omega.
 **
 **	for DISCRETE values, the fields are as follows:
 **		field1 = ' '
 **		field2 = column name
 **		field3 = row name
 **		field4 = value
 **		field6 = prob of value
 **	Options other than discrete are: normal, uniform, exponential,
 **	gamma and geometric.  All distributions share the same field values
 **	for fields 1-3 but use field 4 and 6 differently.  The uses are as
 **	follows:
 **		NORMAL		field4 = mean			field6 = var
 **		UNIFORM		field4 = lowerbound		field6 = upperbound
 **		EXPONENTIAL	field4 = mean
 **		GEOMETRIC	field4 = mean
 \**********************************************************************/

int load_stoch(sdglobal_type* sd_global, one_problem *original, identity *ident,
		char *fname)
{
	FILE *stoch;
	int idx, group_idx, *coeff_idx, i, cnt;
	int num_rvs = 101; /* default # of vals for each continuous rv */
	int blockSize = BLOCKSIZE; /* default blk size for initial omegas arrays */
	int category;
	double percentile = 0.0, base;
	double val1, val2, **mean_dest;
	char field1[15] =
	{ 0 }, field2[60] =
	{ 0 }, field3[15] =
	{ 0 };
	char field4[15] =
	{ 0 }, field5[15] =
	{ 0 }, field6[15] =
	{ 0 }, field7[15] =
	{ 0 };
	char fieldType;
	/* char	fileType[6];*/
	char last_col[16] =
	{ 0 }, last_row[16] =
	{ 0 };
	char name[NAME_SIZE] =
	{ 0 };

	/* modified by Yifan 2013.05.20 */

	/********************\
    1.open stoch file
	 \********************/
	if (!(get_file(&stoch, fname, ".sto", name)))
	{
		printf("failure opening stoch file\n");
		return 0;
	}

	last_col[0] = '\0';
	last_row[0] = '\0';

	idx = 0;
	group_idx = 0;
	coeff_idx = 0;

	/******************************************************************\
		2.initialize data structures for storage of stochastic elements.
	 \******************************************************************/
	idx = 0;
	sd_global->omegas.mapping = NULL;
	if (!(sd_global->omegas.omega_probs =
			(double **) mem_calloc(blockSize, sizeof(double *))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	if (!(sd_global->omegas.omega_vals =
			(double **) mem_calloc(blockSize, sizeof(double *))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	if (!(sd_global->omegas.row = (int *) mem_calloc(blockSize, sizeof(int))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	if (!(sd_global->omegas.col = (int *) mem_calloc(blockSize, sizeof(int))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	if (!(sd_global->omegas.num_vals =
			(int *) mem_calloc(blockSize, sizeof(int))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	if (!(coeff_idx = (int *) mem_calloc(blockSize, sizeof(int))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	if (!(mean_dest = (double **) mem_calloc(blockSize, sizeof(double *))))
	{
		printf("memory allocation error\n");
		return 0;
	}
	sd_global->omegas.num_omega = 0;

	/******************************************************************\
		While data remains in stoch file continue reading and storing
	 of values in struct omegas
	 \******************************************************************/
	while (get_line(&stoch, field1, field2, field3, field4, field5, field6,
			field7, &fieldType))
	{
		if (fieldType == 't')
		{
			if (strcmp(field2, "DISCRETE") == 0)
			{
				category = 1;
			}
			else if (strcmp(field2, "NORMAL") == 0)
			{
				category = 2;
			}
			else if (strcmp(field2, "EXPONENTIAL") == 0)
			{
				category = 3;
			}
			else if (strcmp(field2, "UNIFORM") == 0)
			{
				category = 4;
			}
			else if (strcmp(field2, "GAMMA") == 0)
			{
				category = 5;
			}
			else if (strcmp(field2, "GEOMETRIC") == 0)
			{
				category = 6;
			}
			else if (strcmp(field1, "ENDATA") == 0)
			{
				break;
			}
			else
			{
				continue;
			}
		}

		else if (fieldType == 'f')
		{

			/****************************************************************\
			**	3. New stochastic element encountered, update row and col names,
			 **	verify arrays are sufficiently allocated and check for
			 **	continuous distribution.  If continuous, load next array
			 **	of percentiles.
			 \****************************************************************/
			if (strcmp(field1, last_col) != 0 || strcmp(field2, last_row) != 0)
			{
				idx = sd_global->omegas.num_omega++;
				/* allocate array in struct omegas for next set of vals */
				if (idx == blockSize) /* allocation exceeded, reallocate */
				{
					blockSize = blockSize + BLOCKSIZE;
					if (!(sd_global->omegas.omega_probs =
							(double **) mem_realloc(sd_global->omegas.omega_probs,
									blockSize * sizeof(double *))))
					{
						printf("memory allocation error\n");
						return 0;
					}
					if (!(sd_global->omegas.omega_vals =
							(double **) mem_realloc(sd_global->omegas.omega_vals,
									blockSize * sizeof(double *))))
					{
						printf("memory allocation error\n");
						return 0;
					}
					if (!(sd_global->omegas.row =
							(int *) mem_realloc(sd_global->omegas.row, blockSize * sizeof(int))))
					{
						printf("memory allocation error\n");
						return 0;
					}
					if (!(sd_global->omegas.col =
							(int *) mem_realloc(sd_global->omegas.col, blockSize * sizeof(int))))
					{
						printf("memory allocation error\n");
						return 0;
					}
					if (!(sd_global->omegas.num_vals =
							(int *) mem_realloc(sd_global->omegas.num_vals, blockSize * sizeof(int))))
					/* BETTER INIT TO ZERO!!! THIS ISN'T CALLOC */
					{
						printf("memory allocation error\n");
						return 0;
					}
					if (!(coeff_idx =
							(int *) mem_realloc(coeff_idx, blockSize * sizeof(int))))
					{
						printf("memory allocation error\n");
						return 0;
					}
					if (!(mean_dest =
							(double **) mem_realloc(mean_dest, blockSize * sizeof(double *))))
					{
						printf("memory allocation error\n");
						return 0;
					}

					/* Initialize new counts to zero, after realloc */
					for (i = idx; i < blockSize; i++)
						sd_global->omegas.num_vals[i] = 0;

				}

				/*  4. allocate space for idividual stochastic element */
				if (!(sd_global->omegas.omega_probs[idx] =
						(double *) mem_calloc(num_rvs, sizeof(double))))
				{
					printf("memory allocation error\n");
					return 0;
				}
				if (!(sd_global->omegas.omega_vals[idx] =
						(double *) mem_calloc(num_rvs, sizeof(double))))
				{
					printf("memory allocation error\n");
					return 0;
				}

				/* initialize new arrays to zero */
				for (i = 0; i < num_rvs; i++)
				{
					sd_global->omegas.omega_probs[idx][i] = 0.0;
					sd_global->omegas.omega_vals[idx][i] = 0.0;
				}

				/* update row and column names */
				strcpy(last_col, field1);
				strcpy(last_row, field2);

				/***********************************************\
					5. identify and record row and column numbers
				 \***********************************************/
				/* get row number by name */
				i = 0;
				while (i < original->mar
						&& strcmp(ident->row_name[i], field2) != 0)
				{
					i++;
				}
				if (i >= original->mar)
				{
					if (strcmp(original->objname, field2) == 0)
					{
						i = -1;
					}
					else
					{
						printf("undefined row name in constraint matrix: %s\n",
								field2);
						printf(
								"Randomness in Cost Coefficients and W metrix will be supported in later release\n");
						/*added by Yifan Oct 12 2011*/
						return (0);
					}
				}
				sd_global->omegas.row[idx] = i;

				/* get column number by name */
				i = 0;
				while (i < original->mac
						&& strcmp(ident->col_name[i], field1) != 0)
				{
					i++;
				}

				if (i >= original->mac)
				{
					if (strcmp(field1, "rhs") == 0
							|| strcmp(field1, "RHS") == 0)
						i = -1;
					else
					{
						printf("undefined column name in constraint matrix\n");
						printf(
								"Randomness in Cost Coefficients and W metrix will be supported in later release\n");
						/*added by Yifan Oct 12 2011*/
						return (0);
					}
				}
				sd_global->omegas.col[idx] = i;

				/*********************************************************\
					6. identify and store location of stochastic element
				 in original matval according to col and row numbers
				 obtained above.  Locations are stored in coeff_idx[].
				 \*********************************************************/
				if (i != -1)
				{
					if (sd_global->omegas.row[idx] == -1)
					{
						mean_dest[idx] =
								&(original->objx[sd_global->omegas.col[idx]]);
#ifdef CAL_CHECK
						printf("mean_dest[%d] is %f\n",idx,*mean_dest[idx]);
#endif
					}
					else
					{
						/*original->mac is changed to original->matsz  modified by Yifan 02/07/12*/
						group_idx =
								original->matbeg[sd_global->omegas.col[idx]];

#ifdef CAL_CHECK
						printf("%d\n",idx);
						printf("gropu_idx:%d\n",group_idx);
						printf("sd_global->omegas.col[%d]:%d\n",idx,sd_global->omegas.col[idx]);
#endif

						while (group_idx
								< (original->matbeg[sd_global->omegas.col[idx]]
										+ original->matcnt[sd_global->omegas.col[idx]])
								&& original->matind[group_idx]
										!= sd_global->omegas.row[idx])
						{
							++group_idx;
						}
						coeff_idx[idx] = group_idx;
						mean_dest[idx] = &(original->matval[coeff_idx[idx]]);
					}
				}
				else
					mean_dest[idx] =
							&(original->rhsx[sd_global->omegas.row[idx]]);
#ifdef CAL_CHECK
				printf("coeeff_idx[%d] is %d\n",idx,coeff_idx[idx]);
#endif
				/*****************************************************\
					check for continuous distribution. If continuous,
				 set new array of percentiles.
				 \*****************************************************/
				if (category > 1) /* continuous distribution encountered */
				{
					/* set percentiles in sd_global->omegas.omega_probs */
					base = 1.0 / (float) num_rvs;
					for (i = 0; i < num_rvs - 1; i++)
					{
						percentile += base;
						sd_global->omegas.omega_probs[idx][i] = percentile;
					}
				}
			}

			/****************************************************************\
				7. load values for stochastic elements into array sd_global->omegas.omega_vals.
			 If element is discrete, increment sd_global->omegas.num_vals[idx], else
			 set sd_global->omegas.num_vals[idx] to num_rvs - 1.

			 NOTE: This does not consider continuous r.v.'s for the rhs.
			 \****************************************************************/
			switch (category)
			{
			case 1: /* discrete distribution */
				sd_global->omegas.omega_vals[idx][sd_global->omegas.num_vals[idx]] =
						str_to_float(field3);
				sd_global->omegas.omega_probs[idx][sd_global->omegas.num_vals[idx]] =
						str_to_float(field4);
				++sd_global->omegas.num_vals[idx];
				break;

			case 2: /* Normal distribution */
				val1 = str_to_float(field3);
				val2 = str_to_float(field4);
				normal(sd_global->omegas.omega_probs[idx], num_rvs - 1,
						sd_global->omegas.omega_vals[idx], val1, val2);
				sd_global->omegas.num_vals[idx] = num_rvs - 1;

				/* set value to mean in lp prob. struct */
				*mean_dest[idx] = val1;
				break;

			case 3: /* Exponential distribution */
				val1 = str_to_float(field3);
				exponential(sd_global->omegas.omega_probs[idx], num_rvs - 1,
						sd_global->omegas.omega_vals[idx], val1);
				sd_global->omegas.num_vals[idx] = num_rvs - 1;

				/* set value to mean in lp prob. struct */
				*mean_dest[idx] = val1;
				break;

			case 4: /* Uniform distribution */
				val1 = str_to_float(field3);
				val2 = str_to_float(field4);
				uniform(num_rvs - 1, sd_global->omegas.omega_vals[idx], val1,
						val2);
				sd_global->omegas.num_vals[idx] = num_rvs - 1;

				/* set value to mean in lp prob. struct */
				*mean_dest[idx] = val1;
				break;
				/*
				 case 5:		 Gamma distribution
				 val1 = str_to_float(field3);
				 val2 = str_to_float(field4);
				 gamma(sd_global->omegas.omega_probs[idx], num_rvs - 1,
				 sd_global->omegas.omega_vals[idx], val1, val2);
				 sd_global->omegas.num_vals[idx] = num_rvs - 1;
				 set value to mean in lp prob. struct
				 *mean_dest[idx] = val1;	
				 break; 
				 */

			case 6: /* Geometric distribution */
				val1 = str_to_float(field3);
				geometric(sd_global->omegas.omega_probs[idx], num_rvs - 1,
						sd_global->omegas.omega_vals[idx], val1);
				/* call chg_to_deviation */
				sd_global->omegas.num_vals[idx] = num_rvs - 1;

				/* set value to mean in lp prob. struct */
				*mean_dest[idx] = val1;
				break;
			default:
				break;
			}
		}
	}

	/**********************************************************************************\
   TESTING ... STARTS
	 \**********************************************************************************/
#ifdef CAL_CHECK
	printf("ident->col_name[-1] is %s\n",ident->col_name[0]);
	printf("*************************************\n");
	printf("num_omega = %d;\n",sd_global->omegas.num_omega);
	for (idx=0; idx<sd_global->omegas.num_omega; idx++)
	{
		printf("num_vals[%d] = %d;\n",idx,sd_global->omegas.num_vals[idx]);
		if (sd_global->omegas.col[idx]==-1)
		{
			printf("omega_cname[%d] = \"%s\";\n", idx, "RHS");
		}
		else
		{
			printf("omega_cname[%d] = \"%s\";\n", idx, ident->col_name[sd_global->omegas.col[idx]]);}

		if (sd_global->omegas.row[idx]==-1)
		{
			printf("omega_rname[%d] = \"%s\";\n", idx, "FOBJ");
		}
		else
		{
			printf("omega_rname[%d] = \"%s\";\n", idx, ident->row_name[sd_global->omegas.row[idx]]);}

		for (i=0; i<sd_global->omegas.num_vals[idx]; i++)
		{
			printf("omega_vals[%d][%d] = %f; omega_probs[%d][%d] = %f;\n",idx,i,sd_global->omegas.omega_vals[idx][i],idx,i,sd_global->omegas.omega_probs[idx][i]);
		}
	}
	printf("*************************************\n\n\n");
#endif

	/**********************************************************************************\
   TESTING ... ENDS
	 \**********************************************************************************/

	/**************************************************************\
		8. determine and load means of discrete distribution into
	 lp prob. struct array original.matval[].
	 \**************************************************************/
	for (idx = 0; idx < sd_global->omegas.num_omega; idx++)
	{
		if (sd_global->omegas.num_vals[idx] < num_rvs - 1) /* discrete distribution */
		{
			*mean_dest[idx] = get_mean(sd_global->omegas.omega_vals[idx],
					sd_global->omegas.omega_probs[idx],
					sd_global->omegas.num_vals[idx]);
		}
		/*** Took out a brace   }   right there  -- (Jason) 25 Apr 92 ***/
		/*** Changed the sd_global->omegas.omega_vals to contain deviations from the mean, rather than actual values -- (Jason) 07 Apr 92 ***/
		for (cnt = 0; cnt < sd_global->omegas.num_vals[idx]; cnt++)
			sd_global->omegas.omega_vals[idx][cnt] -= *mean_dest[idx];
	}

	/**********************************************************************************\
   TESTING ... STARTS
	 \**********************************************************************************/
#ifdef CAL_CHECK
	printf("*************************************\n");
	printf("num_omega = %d;\n",sd_global->omegas.num_omega);
	for (idx=0; idx<sd_global->omegas.num_omega; idx++)
	{
		printf("num_vals[%d] = %d;\n",idx,sd_global->omegas.num_vals[idx]);
		for (i=0; i<sd_global->omegas.num_vals[idx]; i++)
		{
			printf("omega_vals[%d][%d] = %f; omega_probs[%d][%d] = %f;\n",idx,i,sd_global->omegas.omega_vals[idx][i],idx,i,sd_global->omegas.omega_probs[idx][i]);
		}
	}
	printf("*************************************\n\n\n");
#endif

	/**********************************************************************************\
   TESTING ... ENDS
	 \**********************************************************************************/

	/**********************************************************************************\
		convert probabilities for discrete distributions into cumulative probs.
	 \**********************************************************************************/
#ifdef CAL_CHECK
	for (idx=0; idx<sd_global->omegas.num_omega; idx++)
	{
		printf("sd_global->omegas.num_vals[%d] is: %d\n",idx,sd_global->omegas.num_vals[idx]); /* modified by Yifan 11/04/2011*/
	}
#endif

	for (idx = 0; idx < sd_global->omegas.num_omega; idx++)
	{
		if (sd_global->omegas.num_vals[idx] < num_rvs) /* discrete distribution */
		{
#ifdef CAL_CHECK
			printf("sd_global->omegas.omega_probs[%d][%d]: %f\n",idx,0,sd_global->omegas.omega_probs[idx][0]);
#endif
			for (i = 1; i < sd_global->omegas.num_vals[idx]; i++)
			{
				sd_global->omegas.omega_probs[idx][i] +=
						sd_global->omegas.omega_probs[idx][i - 1];
#ifdef CAL_CHECK
				printf("sd_global->omegas.omega_probs[%d][%d]: %f\n",idx,i,sd_global->omegas.omega_probs[idx][i]);
#endif
			}
		}
	}

	sort_omegas(sd_global);
	cipher_omegas(sd_global);
  
  /* Yifan 2012.05.21 */
  
  /*  Recording mean value and external omega file's name */
  if(!(sd_global->omegas.mean = (double *)
       mem_calloc(sd_global->omegas.num_omega, sizeof(double))))
  {
    printf("memory allocation error\n");
    return 0;
  }
  if(!(sd_global->omegas.file_name = (char *)
       mem_calloc(NAME_SIZE, sizeof(char))))
  {
    printf("memory allocation error\n");
    return 0;
  }
  for (i = 0; i < sd_global->omegas.num_omega; i++) {
    sd_global->omegas.mean[i] = *mean_dest[i];
  }
  
  
  /* Yifan 2012.05.21 */

	/*added by Yifan to fix the memory leaks 09/30/2011*/
	mem_free(coeff_idx);
	mem_free(mean_dest);

	return (1);
}

int copy_names(char **x_name, char ***xname, char **xstore, unsigned *xstorsz,
		int max)
{
	int i;

	*xstorsz = 0;

	if (!(*xname = arr_alloc(max, string)))
	{
		printf("Allocation error in function copy_names, allocating xname\n");
		return 0;
	}

	if (!(*xstore = arr_alloc(max * WORDSIZE , char)))
	{
		printf("Allocation error in function copy_names, allocating xstore\n");
		return 0;
	}

	for (i = 0; i < max; i++)
	{
		*(*xname + i) = *xstore + *xstorsz;
		strcpy((*xstore + *xstorsz), x_name[i]);
		*xstorsz += strlen(x_name[i]) + 1;
	}

	*xstore = (char *) mem_realloc (*xstore, *xstorsz * sizeof(char));

	return 1;
}

/**********************************************************************\
**	struct identity is used to maintain a record of the order of
 **	columns and rows within the constraint matrix.  the names of the
 **	columns and rows are stored in the arrays col_name and row_name
 **	by index number.
 **
 **	This structure will be loaded by load_core_cpx during the reading of
 **	the core file.  This structure will be used by load omega for
 **	addressing the stochastic elements of the constraint matrix.  The
 **	structure will also be used for the decomposition of the constraint
 **	matrix following the running of the main problem in cplex.
 \**********************************************************************/

int load_time(sd_small *row_num, sd_small *col_num, identity *ident, char *fname)
{
	FILE *time;
	int i;
	char name[15], fieldType, *p, *q;
	char field1[15], field2[60], field3[12];
	char field4[12], field5[12], field6[12];
	char field7[15];

	if (!(get_file(&time, fname, ".tim", name)))
		return 0;
	while (get_line(&time, field1, field2, field3, field4, field5, field6,
			field7, &fieldType))
	{
		if (fieldType == 'f')
		{

			p = field3;
			q = p;
			while (*p != '\0' && *p != '\n')
				*q++ = toupper(*p++);

			if (!(strcmp(field3, "TIME2")))
			{
				/* find row idx */
				i = 0;
				while (i < ident->num_rows && strcmp(field2, ident->row_name[i]))
					i++;
				if (i < ident->num_rows)
					*row_num = i;
				else
				{
					printf("unidentified row name in time file\n");
					printf(
							"Randomness in Cost Coefficients and W metrix will be supported in later release\n");
					/*added by Yifan Oct 12 2011*/
					fclose(time);
					return 0;
				}

				/* find col idx */
				i = 0;
				while (i < ident->num_cols && strcmp(field1, ident->col_name[i]))
					i++;
				if (i < ident->num_cols)
					*col_num = i;
				else
				{
					printf("unidentified col name in time file\n");
					printf(
							"Randomness in Cost Coefficients and W metrix will be supported in later release\n");
					/*added by Yifan Oct 12 2011*/
					fclose(time);
					return 0;
				}

				fclose(time);
				return 1;
			}
		}
	}
	fclose(time);
	return 0;
}

/**********************************************************************\
** This function reads configuration parameters from the configuration
 ** file "config.sd".  These include all the constants required for 
 ** solution by SD.  These are stored in a global structure, config.
 ** It does only very simple error-checking by counting the number of
 ** parameters read.
 \**********************************************************************/
int load_config(sdglobal_type* sd_global, BOOL read_seeds, BOOL read_iters)
{
	FILE *f_in;
	char param[NAME_SIZE * 2];
	char comment[80];
	sd_long long_dummy;
	int int_dummy; /* added by zl. 06/18/02 */
	int sum;
	int x,status=0;

	sum = 0;

    /* Establish default settings for all parameters */
    sd_global->config.EPSILON      = 0.001;
    sd_global->config.SCAN_LEN     = 256;
    
    sd_global->config.MULTIPLE_REP = 0;
    
    sd_global->config.PRINT_CYCLE  = 100;
    sd_global->config.EVAL_RUN_FLAG = 0;
    sd_global->config.EVAL_FLAG    = 0;         	/* JH 3/20/98 */
    sd_global->config.EVAL_ERROR   = 0.01;      	/* JH 3/20/98 */
    sd_global->config.MEAN_DEV     = 0.05;
    
    sd_global->config.OVERRIDE     = 0;
    
    sd_global->config.MASTER_TYPE  = 1;        /* 1 in config for QP master */
    sd_global->config.LB_TYPE      = 1;
    
    sd_global->config.R            = 0.2;
    sd_global->config.R2           = 0.95;       	/* zl 06/04/02 */
    sd_global->config.R3           = 2.0;       	/* zl 06/04/02 */
    
    sd_global->config.MIN_QUAD_SCALAR = 0.001; /*Yifan 03/27*/
    sd_global->config.MAX_QUAD_SCALAR = 10000.0; 	/* zl 06/20/02 */
    
    sd_global->config.TAU          = 2;
    
    sd_global->config.PRE_EPSILON  = 0.01;     	/* for pre_tests. zl */
    sd_global->config.ITER_FACT    = 0;       	/* JH 3/13/98 */
    sd_global->config.PERCENT_PASS = 0.95;
    sd_global->config.M            = 50;         	/* zl trial. */
    sd_global->config.CUT_MULT     = 5;
    
    
    sd_global->config.CONFID_HI    = 1.0;
    sd_global->config.CONFID_LO    = 1.45;
    sd_global->config.TOLERANCE    = 0.001;
    sd_global->config.FEA_TOLER    = 0.05;
    sd_global->config.THIN_TOLER   = 0.001;
    sd_global->config.START_THIN   = 9001;
    sd_global->config.THIN_CYCLE   = 200;
    sd_global->config.DROP_TIME    = 16300;
    sd_global->config.TEST_TYPE    = 1;        /* 1 in config for full test */
    /* 0 in config for LP master */
    
    sd_global->config.PI_EVAL_START = 1;
    sd_global->config.DETAILED_SOLN = 1;
    sd_global->config.PI_CYCLE      = 1;
    sd_global->config.SUB_LB_CHECK  = 0;			/* Subprob LB check. 0 for no check.
                                       zl 09/20/05 */
    sd_global->config.AUTO_SEED     = 0;
    sd_global->config.BOX_SIZE      = 1;    /* modified by Yifan 2013.9.13 Default box size */

    
    sd_global->config.SMOOTH_I     = 50;
    sd_global->config.SMOOTH_PARM  = 0.25;
    sd_global->config.SUBPROB_LB   = 0.0;		/* LB on subprob obj values,
                                       zl 07/01/04. */
    if (read_iters)                  	/* added by zl. 06/18/02. */
    {
      sd_global->config.MIN_ITER     = 0;
      sd_global->config.MAX_ITER     = 5000;
    }

	if (read_seeds)
	{
		sd_global->config.RUN_SEED1 = -1;
		sd_global->config.RUN_SEED2 = -1;
		sd_global->config.EVAL_SEED1 = -1;
	}

	f_in = fopen("config.sd", "r");

	if (f_in == NULL)
	{
		printf("Unable to open configuration file config.sd.\n");
		printf("Using default parameters.\n");
	}
	else
	{
		/* Read in any settings from the configuration file */
		while ((x = (fscanf(f_in, "%s", param) != EOF)))
		{
			if (!strcmp(param, "MASTER_TYPE"))
				status = fscanf(f_in, "%d", &(sd_global->config.MASTER_TYPE));
			else if (!strcmp(param, "TEST_TYPE"))
				status = fscanf(f_in, "%d", &(sd_global->config.TEST_TYPE));
			else if (!strcmp(param, "CONFID_HI"))
				status = fscanf(f_in, "%lf", &(sd_global->config.CONFID_HI));
			else if (!strcmp(param, "CONFID_LO"))
				status = fscanf(f_in, "%lf", &(sd_global->config.CONFID_LO));
			else if (!strcmp(param, "MASTER_TYPE"))
				status = fscanf(f_in, "%d", &(sd_global->config.MASTER_TYPE));
			else if (!strcmp(param, "LB_TYPE"))
				status = fscanf(f_in, "%d", &(sd_global->config.LB_TYPE));
			else if (!strcmp(param, "PERCENT_PASS"))
				status = fscanf(f_in, "%lf", &(sd_global->config.PERCENT_PASS));
			else if (!strcmp(param, "EPSILON"))
				status = fscanf(f_in, "%lf", &(sd_global->config.EPSILON));
			else if (!strcmp(param, "PRE_EPSILON"))
				status = fscanf(f_in, "%lf", &(sd_global->config.PRE_EPSILON));
			else if (!strcmp(param, "SMOOTH_PARM"))
				status = fscanf(f_in, "%lf", &(sd_global->config.SMOOTH_PARM));
			else if (!strcmp(param, "SMOOTH_I"))
				status = fscanf(f_in, "%d", &(sd_global->config.SMOOTH_I));
			else if (!strcmp(param, "TOLERANCE"))
				status = fscanf(f_in, "%lf", &(sd_global->config.TOLERANCE));
			else if (!strcmp(param, "FEA_TOLER"))
				status = fscanf(f_in, "%lf", &(sd_global->config.FEA_TOLER));
			else if (!strcmp(param, "THIN_TOLER"))
				status = fscanf(f_in, "%lf", &(sd_global->config.THIN_TOLER));
			else if (!strcmp(param, "R"))
				status = fscanf(f_in, "%lf", &(sd_global->config.R));
			else if (!strcmp(param, "R2"))
				status = fscanf(f_in, "%lf", &(sd_global->config.R2));
			else if (!strcmp(param, "R3"))
				status = fscanf(f_in, "%lf", &(sd_global->config.R3));
			else if (!strcmp(param, "MIN_QUAD_SCALAR"))
				status = fscanf(f_in, "%lf", &(sd_global->config.MIN_QUAD_SCALAR));
			else if (!strcmp(param, "MAX_QUAD_SCALAR"))
				status = fscanf(f_in, "%lf", &(sd_global->config.MAX_QUAD_SCALAR));
			else if (!strcmp(param, "CUT_MULT"))
				status = fscanf(f_in, "%d", &(sd_global->config.CUT_MULT));
			else if (!strcmp(param, "ITER_FACT"))
				status = fscanf(f_in, "%lf", &(sd_global->config.ITER_FACT)); /* JH 3/13/98 */
			else if (!strcmp(param, "EVAL_FLAG"))
				status = fscanf(f_in, "%d", &(sd_global->config.EVAL_FLAG)); /* JH 3/20/98 */
			else if (!strcmp(param, "EVAL_ERROR"))
				status = fscanf(f_in, "%lf", &(sd_global->config.EVAL_ERROR)); /* JH 3/20/98 */
			else if (!strcmp(param, "MEAN_DEV"))
				status = fscanf(f_in, "%lf", &(sd_global->config.MEAN_DEV)); /* JH 3/20/98 */
			else if (!strcmp(param, "M"))
				status = fscanf(f_in, "%d", &(sd_global->config.M));
			else if (!strcmp(param, "TAU"))
				status = fscanf(f_in, "%d", &(sd_global->config.TAU));
			else if (!strcmp(param, "MIN_ITER"))
				if (read_iters)
					status = fscanf(f_in, "%d", &(sd_global->config.MIN_ITER)); /* zl 06/18/02 */
				else
					status = fscanf(f_in, "%d", &(int_dummy));
			else if (!strcmp(param, "MAX_ITER"))
				if (read_iters)
				{
					status = fscanf(f_in, "%d", &(sd_global->config.MAX_ITER)); /* zl 06/18/02 */
					sd_global->config.START_THIN = sd_global->config.MAX_ITER;
				}
				else
					status = fscanf(f_in, "%d", &int_dummy);
            else if (!strcmp(param, "OVERRIDE"))
                status = fscanf(f_in, "%d", &(sd_global->config.OVERRIDE));
			else if (!strcmp(param, "START_THIN"))
				status = fscanf(f_in, "%d", &(sd_global->config.START_THIN));
			else if (!strcmp(param, "THIN_CYCLE"))
				status = fscanf(f_in, "%d", &(sd_global->config.THIN_CYCLE));
			else if (!strcmp(param, "PRINT_CYCLE"))
				status = fscanf(f_in, "%d", &(sd_global->config.PRINT_CYCLE));
			else if (!strcmp(param, "EVAL_RUN_FLAG"))
				status = fscanf(f_in, "%d", &(sd_global->config.EVAL_RUN_FLAG));
			else if (!strcmp(param, "DROP_TIME"))
				status = fscanf(f_in, "%d", &(sd_global->config.DROP_TIME));
			else if (!strcmp(param, "RUN_SEED1"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED1));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED2"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED2));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED3"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED3));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED4"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED4));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED5"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED5));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED6"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED6));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED7"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED7));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED8"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED8));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED9"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED9));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED10"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED10));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED11"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED11));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED12"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED12));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED13"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED13));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED14"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED14));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED15"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED15));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED16"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED16));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED17"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED17));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED18"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED18));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED19"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED19));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED20"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED20));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED21"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED21));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED22"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED22));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED23"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED23));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED24"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED24));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED25"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED25));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED26"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED26));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED27"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED27));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED28"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED28));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED29"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED29));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "RUN_SEED30"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.RUN_SEED30));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "EVAL_SEED1"))
				if (read_seeds)
					status = fscanf(f_in, "%lld", &(sd_global->config.EVAL_SEED1));
				else
					status = fscanf(f_in, "%lld", &long_dummy);
			else if (!strcmp(param, "SUB_LB_CHECK"))
				status = fscanf(f_in, "%d", &(sd_global->config.SUB_LB_CHECK));
			else if (!strcmp(param, "PI_EVAL_START"))
				status = fscanf(f_in, "%d", &(sd_global->config.PI_EVAL_START));
			else if (!strcmp(param, "PI_CYCLE"))
				status = fscanf(f_in, "%d", &(sd_global->config.PI_CYCLE));
			else if (!strcmp(param, "SCAN_LEN"))
				status = fscanf(f_in, "%d", &(sd_global->config.SCAN_LEN));
			else if (!strcmp(param, "DETAILED_SOLN"))
				status = fscanf(f_in, "%d", &(sd_global->config.DETAILED_SOLN));
			else if (!strcmp(param, "MULTIPLE_REP"))
				status = fscanf(f_in, "%d", &(sd_global->config.MULTIPLE_REP));
			else if (!strcmp(param, "AUTO_SEED"))
				status = fscanf(f_in, "%d", &(sd_global->config.AUTO_SEED));
            else if (!strcmp(param, "BOX_SIZE"))
				status = fscanf(f_in, "%lf", &(sd_global->config.BOX_SIZE));
			else if (!strcmp(param, "//"))
			{
              if (fgets(comment, 80, f_in) != NULL) {
                 sum--;
              }              
			}
			else
			{
				printf("Unrecognized configuration parameter: %s.\n", param);
				return 0;
			}
            if (status < 0) {
              printf("Error with config file, no parameter value");
            }
			sum++;
		}

		if (sum < 17)
		{
			/* modified by zl for online experimental runs, because we only allow changes 
			 on limited number of parameters for online users. 12/10/02. 

			 printf("Warning: Only %d configuration parameters read from config.sd.\n",
			 sum);
			 */
			printf("Using default parameters.\n");
		}

		fclose(f_in);
	}

	if (sd_global->config.RUN_SEED1 < 0 || sd_global->config.RUN_SEED2 < 0
			|| sd_global->config.EVAL_SEED1 < 0)
	{
		printf("Random number seeds not specified in config.sd.\n");
		printf("Please enter run seeds 1 and 2, evaluation seed 1: ");
		status = scanf("%lld %lld %lld", &(sd_global->config.RUN_SEED1), &(sd_global->config.RUN_SEED2), &(sd_global->config.EVAL_SEED1));
        if (status) {
          printf("Error reading seed from config file.");
          exit(1);
        }
	}
	return 1;
}

/**********************************************************************\
**	the function get_file receives a string representing the variety of
 **	file to be obtained and a FILE pointer.  The function then obtains
 **	the name of the file to be opened from the user and attempts to 
 **	open the file.  If the file opens successfully, the function returns
 **	a value of 1, 0 otherwise.
 **
 ** 16 May 92
 ** Changed this function so that it opens the specified filename and
 ** extension (for suite-testing purposes).  
 \**********************************************************************/

int get_file(FILE **fptr, char *fname, char *extension, char *probname)
{
	char file_name[STRSIZE], *p, *q;

	strcpy(file_name, "./sdinput/");
	strcat(file_name, fname);
	strcat(file_name, "/");
	strcat(file_name, fname);
	strcat(file_name, extension);

	*fptr = fopen(file_name, "r");
	if (*fptr == NULL)
	{
		printf("Unable to open problem file: %s.\n", file_name);
		printf("Please check for presence of file and try again.");
		return 0;
	}

	/* remove any .*** extension from file name and copy into probname */
	p = file_name;
	q = probname;
	while (*p != '.' && *p != '\0')
		*q++ = *p++;
	*q = '\0';

	return 1;
}

/**********************************************************************\
**	The function str_to_float converts a string passed to it as an 
 **	argument and returns a double.  The string may contain an exponent
 **	part as well as a real part.
 \**********************************************************************/

double str_to_float(char *string)
{
	double val;

	sscanf(string, "%lf", &val);

	return val;
}

/************************************************************************\
**	The function get_mean() receives an array of probabilities and an
 **	array of corresponding values for a discrete distribution along
 **	with an integer representing the number of values contained in each
 **	array.  The function then calculates and returns the mean of the 
 **	distribution.
 \************************************************************************/

double get_mean(double *vals, double *probs, int num)
{
	int i;
	double mean;

	mean = 0.0;
	for (i = 0; i < num; i++)
	{
		mean += vals[i] * probs[i];
	}

	return (mean);
}

void uniform(int num_rvs, double *val_array, double val1, double val2)
{
	int i;
	double increment, val;

	val = val1;
	increment = (val2 - val1) / (float) num_rvs;

	for (i = 0; i <= num_rvs; i++)
	{
		val += increment;
		val_array[i] = val;
	}
}

void free_ident(identity *ident)
{
	int r, c;
	for (r = 0; r < ident->num_rows; r++)
		mem_free(ident->row_name[r]);

	for (c = 0; c < ident->num_cols; c++)
		mem_free(ident->col_name[c]);

	mem_free(ident->row_name);
	mem_free(ident->col_name);

	mem_free(ident);
}

void free_omegas(sdglobal_type* sd_global)
{
	int r;

	mem_free(sd_global->omegas.row);
	mem_free(sd_global->omegas.col);
	mem_free(sd_global->omegas.num_vals);
	mem_free(sd_global->omegas.indices);
	mem_free(sd_global->omegas.key);

	for (r = 0; r < sd_global->omegas.num_omega; r++)
	{
		mem_free(sd_global->omegas.omega_vals[r]);
		mem_free(sd_global->omegas.omega_probs[r]);
	}

	mem_free(sd_global->omegas.omega_vals);
	mem_free(sd_global->omegas.omega_probs);
    mem_free(sd_global->omegas.mean);
    mem_free(sd_global->omegas.file_name);
}

/*************************************************************************\
** The function load_core_cpx reads the input file in mps format.  It then 
 ** stores the data in the struct 'original' of type one_problem.  zl
 \*************************************************************************/
int load_core_cpx(one_problem **original, identity **ident, char *fname,
		int objsen)
{
	/* int	i;*//* default increment counter */
	/* int	wordSize = WORDSIZE;*//* default name length for row and col. arrays */
	char name[STRSIZE];

	/* Allocate structure one_problem 'original'. */
	if (!(*original = (one_problem *) mem_malloc(sizeof(one_problem))))
	{
		printf("memory allocation failure, one_problem, in load_core_cpx\n");
		return 0;
	}

	/* Allocate structure identity 'ident'. */
	if (!(*ident = (identity *) mem_malloc(sizeof(identity))))
	{
		printf("Memory allocation failure, ident, in load_core_cpx\n");
		return 0;
	}

	/* allocate and load problem name. */
	if (!((*original)->name = (char *) mem_calloc(NAME_SIZE, sizeof(char))))
	{
		printf("Memory allocation failure, original, in load_core_cpx\n");
		return 0;
	}

	strcpy((*original)->name, fname);
	(*original)->objsen = objsen;

#ifdef LOAD_CORE_CPX
	printf("zl_ldcore ori->name = %s, ori->objsen = %d\n", (*original)->name,
			(*original)->objsen);
#endif

	/* read COR file from input folder modified by Yifan 2013.05.20 */
	/* the folder structure looks like this ./sdinput/prob_name/prob_name.cor */

	strcpy(name, "./sdinput/");
	strcat(name, fname);
	strcat(name, "/");
	strcat(name, fname);
#ifdef CPLEX
	strcat(name, ".cor");
#else
  strcat(name, ".mps");
#endif
	printf("Reading problems from %s \n", name);

	/* Read the problem in external Solver. */
	(*original)->lp = read_problem((*original), name, "MPS"); /* 2011.10.30 */
	if ((*original)->lp == NULL)
	{
		printf(" load_core_cpx: rats, original wasn't created! \n");
		err_msg("Loading COR file", "load_core_cpx", "read_problem");
	}

	/* Initialize all parameters that won't be used */
	(*original)->mae = 0;
	(*original)->maesz = 0;
	(*original)->enzsz = 0;
	(*original)->estorsz = 0;
	(*original)->rngcol = NULL;
	(*original)->nrowind = NULL;
	(*original)->etype = NULL;
	(*original)->ename = NULL;
	(*original)->estore = NULL;
	(*original)->enzbeg = NULL;
	(*original)->enzcnt = NULL;
	(*original)->enzind = NULL;
	(*original)->enzval = NULL;
	(*original)->dataname = NULL;
	(*original)->rhsname = NULL;
	(*original)->rngname = NULL;
	(*original)->bndname = NULL;

	allocate_arrays(original, ident);

#if 0
	(*ident)->num_rows = (*original)->mar;
	(*ident)->num_cols = (*original)->mac;
#endif

	printf("finished loading core file via external Solver.\n");
	return (1);
}

/******************************************************************************\
** This function allocates spaces for the arrays in structs 'orig' (one_problem
 ** type) and 'ident' (identity type), then obtain values for them from external Solver. 
 \******************************************************************************/
void allocate_arrays(one_problem **orig, identity **ident)
{
	int cols, rows, nzcnt, rstorsz, cstorsz;
	int i, status, surplus;

#ifdef LOAD_CORE_CPX
	printf("zl_allocate_arr ori->name = %s, ori->objsen = %d\n", (*orig)->name,
			(*orig)->objsen);
#endif

	/* 2011.10.30 */
	cols = get_numcols((*orig));
	rows = get_numrows((*orig));
	nzcnt = get_numnz((*orig));

	(*orig)->mac = cols;
	(*orig)->mar = rows;
	(*orig)->macsz = nzcnt;
	(*orig)->marsz = nzcnt;
	(*orig)->matsz = nzcnt;

#ifdef LOAD_CORE_CPX
	printf("zl_alloc_arr cols = %d, rows = %d, nzcnt = %d\n", (*orig)->mac,
			(*orig)->mar, (*orig)->matsz);
#endif

	/* Allocate spaces for objx, rhsx, senx, bdl, and bdu of struct 'orig'. */
	(*orig)->objx = (double *) calloc(cols, sizeof(double));
	(*orig)->rhsx = (double *) calloc(rows, sizeof(double));
	(*orig)->senx = (char *) calloc(rows, sizeof(char));
    (*orig)->typex = (char *) calloc(cols, sizeof(char));
	(*orig)->bdl = (double *) calloc(cols, sizeof(double));
	(*orig)->bdu = (double *) calloc(cols, sizeof(double));

	/* Obtain objective function coefficiences, objx, from external Solver. */
	status = get_obj((*orig), (*orig)->objx, 0, cols - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_obj failed in load_core_cpx\n");

#ifdef LOAD_CORE_CPX
	for (i=0; i<cols; i++)
	printf("ori->objx[%d] = %f\n", i, (*orig)->objx[i]);
#endif

	/* Obtain rhs values, rhsx, from external Solver. */
	status = get_rhs((*orig), (*orig)->rhsx, 0, rows - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_rhs failed in load_core_cpx\n");
    status = get_type(*orig, (*orig)->typex, 0, cols-1); /* Yifan 2013.09.12 */
	if (status != 0)
		printf("get_type failed in load_core_cpx\n");

#ifdef LOAD_CORE_CPX
	for (i=0; i<rows; i++)
	printf("ori->rhsx[%d] = %f\n", i, (*orig)->rhsx[i]);
#endif

	/* Obtain the senses of constraints, senx, from external Solver. */
	status = get_sense((*orig), (*orig)->senx, 0, rows - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_sense failed in load_core_cpx\n");

#ifdef LOAD_CORE_CPX
	for (i=0; i<rows; i++)
	printf("ori->sense[%d] = %c\n", i, (*orig)->senx[i]);
#endif

	/* Obtain the bounds of variables, bdl and bdu, from external Solver. */
	status = get_lbound((*orig), (*orig)->bdl, 0, cols - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_lbound failed in load_core_cpx\n");

	status = get_ubound((*orig), (*orig)->bdu, 0, cols - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_ubound failed in load_core_cpx\n");

#ifdef LOAD_CORE_CPX
	for (i=0; i<rows; i++)
	printf("ori->bdl[%d] = %f, ori->bdu[%d] = %f\n", i, (*orig)->bdl[i], i,
			(*orig)->bdu[i]);
#endif

	/* Allocate spaces for columns of constraints. */
	(*orig)->matbeg = (int *) calloc(cols, sizeof(int));
	(*orig)->matcnt = (int *) calloc(cols, sizeof(int));
	(*orig)->matval = (double *) calloc(nzcnt, sizeof(double));
	(*orig)->matind = (int *) calloc(nzcnt, sizeof(int));

	/* Obtain matbeg, matval, and matind from external Solver. */
	status = get_cols((*orig), &nzcnt, (*orig)->matbeg, (*orig)->matind,
			(*orig)->matval, nzcnt, &surplus, 0, cols - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_cols failed in load_core_cpx\n");

	/* Calculate matcnt. */
	for (i = 0; i < cols - 1; i++)
		(*orig)->matcnt[i] = (*orig)->matbeg[i + 1] - (*orig)->matbeg[i];
	(*orig)->matcnt[cols - 1] = nzcnt - (*orig)->matbeg[cols - 1];

	/* Allocate space for objective name and then obtain it from external Solver. */
	(*orig)->objname = (char *) calloc(WORDSIZE, sizeof(char));
	status = get_objname((*orig), (*orig)->objname, WORDSIZE, &surplus); /* 2011.10.30 */
	if (status != 0)
		printf("get_objname failed in load_core_cpx\n");

#ifdef LOAD_CORE_CPX
	printf("objname = %s\n", (*orig)->objname);
#endif

	/* Obtain rstorsz and cstorsz from external Solver. ??what is the algorithm for namespace-I will use 18*numrow,18*nuncol-asked by Yifan 11/1/2011??*/
#ifdef CPLEX
	getnamespace((*orig), &rstorsz, &cstorsz, rows, cols);
#else
    rstorsz = 18*rows;
    cstorsz = 18*cols;  
#endif
	(*orig)->rstorsz = rstorsz;
	(*orig)->cstorsz = cstorsz;

#ifdef LOAD_CORE_CPX
	printf("rstorsz = %d, cstorsz = %d\n", (*orig)->rstorsz, (*orig)->cstorsz);
#endif

	/* Allocate spaces for rname, rstore, cname, cstore of struct 'orig'. */
	if ((cstorsz > 0) && (rstorsz > 0))
	{
		(*orig)->rname = (char **) calloc(rows, sizeof(char*));
		(*orig)->rstore = (char *) calloc(rstorsz, sizeof(char));
		(*orig)->cname = (char **) calloc(cols, sizeof(char*));
		(*orig)->cstore = (char *) calloc(cstorsz, sizeof(char));
		if ((*orig)->rname == NULL || (*orig)->rstore == NULL)
			printf("Row Name memory allocation failed in load_core_cpx.\n");
		if ((*orig)->cname == NULL || (*orig)->cstore == NULL)
			printf("Column Name memory allocation failed in load_core_cpx\n");
	}

	/* Obtain row names (rname, rstore) from external Solver. */
	status = get_rowname((*orig), (*orig)->rname, (*orig)->rstore, rstorsz,
			&surplus, 0, rows - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_rowname failed in load_core_cpx\n");

	/* Obtain column names (cname, cstore) from external Solver. */
	status = get_colname((*orig), (*orig)->cname, (*orig)->cstore, cstorsz,
			&surplus, 0, cols - 1); /* 2011.10.30 */
	if (status != 0)
		printf("get_colname failed in load_core_cpx\n");

#ifdef LOAD_CORE_CPX
	for (i=0; i<rows; i++)
	printf("ori-rname[%d] = %s\n", i, (*orig)->rname[i]);
	for (i=0; i<cols; i++)
	printf("ori-cname[%d] = %s\n", i, (*orig)->cname[i]);
#endif

	/* Copy num_rows and num_cols of ident from orig. ??why need ident--asked by Yifan 11/1/2011??*/
	(*ident)->num_rows = (*orig)->mar;
	(*ident)->num_cols = (*orig)->mac;

#ifdef LOAD_CORE_CPX
	printf("id->num_rows = %d, id->num_cols = %d\n", (*ident)->num_rows,
			(*ident)->num_cols);
#endif

	/* Allocate spaces for row_name in struct 'ident', and then copy names 
	 from struct 'orig'. */
	(*ident)->row_name = (char **) calloc(rows, sizeof(char *));
	if ((*ident)->row_name == NULL)
		printf("ident->row_name memory allocation failed in load_core_cpx.\n");
	for (i = 0; i < rows; i++)
	{
		(*ident)->row_name[i] = (char *) calloc(WORDSIZE, sizeof(char));
		if ((*ident)->row_name[i] == NULL)
			printf("memory allocation error, ident->row_name[%d].\n", i);
		strcpy((*ident)->row_name[i], (*orig)->rname[i]);

#ifdef LOAD_CORE_CPX
		printf("id->rname[%d] = %s, ori->rname[%d] = %s\n", i, (*ident)->row_name[i], i, (*orig)->rname[i]);
#endif
	}

	/* Allocate spaces for col_name in struct 'ident', and then copy names 
	 from struct 'original'. */
	(*ident)->col_name = (char **) calloc(cols, sizeof(char *));
	if ((*ident)->col_name == NULL)
		printf("ident->col_name memory allocation failed in load_core_cpx.\n");
	for (i = 0; i < cols; i++)
	{
		(*ident)->col_name[i] = (char *) calloc(WORDSIZE, sizeof(char));
		if ((*ident)->col_name[i] == NULL)
			printf("memory allocation error, ident->col_name[%d].\n", i);
		strcpy((*ident)->col_name[i], (*orig)->cname[i]);

#ifdef LOAD_CORE_CPX
		printf("id->cname[%d] = %s, ori->cname[%d] = %s\n", i, (*ident)->col_name[i],
				i, (*orig)->cname[i]);
#endif
	}
} /*** END allocate_arrays ***/

void getnamespace(one_problem *p, int *rs, int *cs, int numrows, int numcols) /* 2011.10.30 */
{
	int status;
	/* int namespace; */
	int surplus;
	int errors = 0;

	status = get_rowname(p, NULL, NULL, 0, &surplus, 0, numrows - 1); /* 2011.10.30 */
	if ((status != NEGATIVE_SURPLUS) && (status != 0))
	{
		printf(" Could not determine amount of space necessary \n");
		++errors;
	}
	*rs = -surplus;

	/* added by Yifan to get the rstorsz 11/1/2011*/
	/*
	 printf(" In getname space for row\n");
	 printf(" This is  numrows:%d\n",numrows);
	 printf(" This is  *rs:%d\n",*rs);
	 */

	status = get_colname(p, NULL, NULL, 0, &surplus, 0, numcols - 1); /* 2011.10.30 */
	if ((status != NEGATIVE_SURPLUS) && (status != 0))
	{
		printf(" Could not determine amount of space necessary \n");
		++errors;
	}
	*cs = -surplus;

	/* added by Yifan to get the rstorsz 11/1/2011*/
	/*
	 printf(" In getname space for column\n");
	 printf(" This is  numcols:%d\n",numcols);
	 printf(" This is  *cs:%d\n",*cs);
	 */

} /** End getnamespace **/
