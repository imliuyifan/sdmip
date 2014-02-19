/*************************************************************************\
**	This file contains the function parser() which separates a 
 **	string from an input file into it's appropriate fields according
 **	to column values:  
 **
 **		Title fields:	columns 1 to 14 and 15 to 72
 **		Info fields:	field1  columns	2 and 3     code field
 **                              field2  columns	5-12        first name field
 **                              field3	columns	15-22       second name field
 **                              field4	columns 25-36       first numeric field
 **                              field5	columns 40-47       third name field
 **                              field6	columns 50-61       second numeric field
 **
 **	Functions in this file are called by all input functions for reading
 **	in core, stoch and cell files and do not contain any information 
 **	particular to any function.  All special characteristics are handled
 **	by the calling functions.
 \*************************************************************************/

#include <string.h>
#include "prob.h"
#include "parser.h"
#include "sdglobal.h"

/*************************************************************************\
**	The function get_line() reads in the next input line, checks its length,
 **	and determines the appropriate parsing function to call.  In the case of
 **	a line consisting of only a single carriage return, the function 
 **	continues to read in lines until a string of non_zero length or an EOF
 **	is encountered. 
 **
 **	Upon return of fields from the parsing functions, each field will be 
 **	sent to remove_spaces for the removal of all blank spaces.  The get_line
 **	function then returns control to the calling function.
 \*************************************************************************/

int get_line(FILE **input, char *field1, char *field2, char *field3,
		char *field4, char *field5, char *field6, char *field7, char *type)
{
	char input_str[100], *strptr;
	int num_fields, i;
    sd_long len;

	for (i = 0; i < 100; i++)
		input_str[i] = '\0';

	len = 0;
	input_str[0] = '*';

	/*  test for and skip over empty lines and comments  */
	while (len <= 1 || input_str[0] == '*')
	{
		strptr = fgets(input_str, 100, *input);
		if (strptr == NULL)
		{
			return 0;
		}

		len = strlen(input_str);
	}

	/*  identify type of input string (title or field)  */
	if (input_str[0] >= '0' && input_str[0] <= 'Z')
	{
		num_fields = title(input_str, field1, field2);
		type[0] = 't'; /* input string is title string */
	}
	else
	{
		num_fields = get_fields(input_str, field1, field2, field3, field4,
				field5, field6, field7);
		type[0] = 'f'; /* input string is field string */
	}

	return (num_fields);
}

/*************************************************************************\
**	The function get_fields receives a pointer to the input string and
 **	pointers to the input fields field1 to field6.  The function returns
 **	an integer representing tyhe last field loaded.
 **
 **	Internally, the function parses the input string into six strings 
 **	representing the input fields as defined by the International Institute
 **	for Applied Systems Analysis WP-87-118, December 1987.
 \*************************************************************************/

int get_fields(char *input_str, char *field1, char *field2, char *field3,
		char *field4, char *field5, char *field6, char *field7)
{
	int end, loc, loc2, last_field;
	int field_num;
    sd_long len;

	char *fields[7];
	char *current, lastCh;

	/* initialize field number to be loaded */
	field_num = 0;

	/* initialize last field loaded to zero */
	last_field = 0;

	/* initialize array of string pointers */
	fields[0] = field1;
	fields[1] = field2;
	fields[2] = field3;
	fields[3] = field4;
	fields[4] = field5;
	fields[5] = field6;
	fields[6] = field7;

	/* initialize location of first character in input string to be read */
	end = 0;
	loc = 0;
	loc2 = 0;
	lastCh = ' ';

	/* obtain pointer to first field to be loaded */
	current = fields[field_num];

	/* get length of string to be parsed */
	len = strlen(input_str);

	while (loc < len)
	{
		/* check for end of string */
		if (input_str[loc] == '\0' || input_str[loc] == '\n')
		{
			end = 1;
		}

		/*  Load characters into field */

		/*  if tab encountered, mark end of string */
		else if (input_str[loc] == '\t' && lastCh != ' ')
		{
			end = 1;
			lastCh = ' ';
		}
		/* if blank space follows valid character, mark end of string*/
		else if (lastCh >= 33 && input_str[loc] == ' ')
		{
			end = 1;
			lastCh = ' ';
		}
		/* else load next character */
		else if (input_str[loc] >= 33)
		{
			current[loc2] = input_str[loc];
			lastCh = current[loc2++];
		}

		/* check for end of current field and begin next */
		if (end == 1)
		{
			end = 0;
			last_field++;
			current[loc2] = '\0';
			++field_num;
			if (field_num < 7)
			{
				loc2 = 0;
				current = fields[field_num];
			}
		}

		/* increment location index to input string */
		++loc;
	}

	/* eliminate additional spaces from fields */
	for (field_num = 0; field_num < last_field; field_num++)
	{
		remove_spaces(fields[field_num]);
	}

	return (last_field);
}

/*************************************************************************\
**	The function title() receives a pointer to the input string and
 **	pointers to the input fields field1 and field2.  The function returns
 **	an integer representing the last field loaded.
 **
 **	Internally, the function parses the input string into two strings 
 **	representing the input fields as defined by the International Institute
 **	for Applied Systems Analysis WP-87-118, December 1987.
 \*************************************************************************/

int title(char *input_str, char *field1, char *field2)
{
	/* Yifan 03/24/2012 endpt is not correctlly defined!*/
	int i, last_field=0, end, end1 = 0, end2 = 0, endpt = 0;
	end = 0;
	for (i = 0; i < 72; i++)
	{
		if (input_str[i] == '\0' || input_str[i] == '\n')
		{
			end = 1;
			end1 = 1;
			end2 = 1;
		}

		/*  Load first string into field 1  */

		if (!end1)
		{
			/*  if tab encountered, set end1 to 1 */
			if (input_str[i] == '\t' || input_str[i] == ' ')
			{
				field1[i] = '\0';
				end1 = 1;
				endpt = i + 1;
			}
			else if (i == 13)
			{
				/* field1[i] = '\n'; modified by zl */
				field1[i] = '\0';
				end1 = 1;
				endpt = i + 1;
			}
			else
			{
				field1[i] = input_str[i];
			}
			last_field = 1;
		}

		/*  Load second string in title into field 2  */
		if (end1)
			if (!end2)
			{
				/*  if tab encountered, replace with space */
				if (input_str[i] == '\t')
				{
					if (i >= endpt)
					{
						field2[i - endpt] = ' ';
					}
				}
				else
				{
					if (i >= endpt)
					{
						field2[i - endpt] = input_str[i];
					}
				}
				last_field = 2;
			}

		/*  end of string obtained, set conditions to exit loop.  */
		if (end == 1)
		{
			endpt = i - endpt;
			i = 72;
		}
	}

	/* Terminate field2 */
	if (last_field == 2)
	{
		field2[endpt] = '\0';
	}

	/* remove additional spaces from fields */
	if (last_field > 0)
	{
		remove_spaces(field1);
	}
	if (last_field > 1)
	{
		remove_spaces(field2);
	}

	return (last_field);
}

/******************************************************************\
**	The function remove_space() removes additional spaces from a
 **	string.  This function will be used after an input string has
 **	been broken into its appropriate fields according to column
 **	values.
 \******************************************************************/

void remove_spaces(char *field)
{
	char *p, *q;

	p = field;
	q = p;
	while (*p != '\0' && *p != '\n')
	{
		if (*p >= 33)
		{
			*q++ = *p++;
		}
		else
		{
			p++;
		}
	}
	*q = '\0';
	return;
}
