/*
 * parser.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef PARSER_H_
#define PARSER_H_

/*************************************************************************\
**  The function get_line() reads in the next input line, checks it's length,
 ** and determines the appropriate parsing function to call.  In the case of
 ** a line consisting of only a single carriage return, the function
 ** continues to read in lines until a string of non_zero length or an EOF
 ** is encountered.
 **
 ** Upon return of fields from the parsing functions, each field will be
 ** sent to remove_spaces for the removal of all blank spaces.  The get_line
 ** function then returns control to the calling function.
 \*************************************************************************/

int get_line(FILE **, char *, char *, char *, char *, char *, char *, char *,
		char *);

/*************************************************************************\
**  The function get_fields receives a pointer to the input string and
 ** pointers to the input fields field1 to field6.  The function returns
 ** an integer representing tyhe last field loaded.
 **
 **  Internally, the function parses the input string into six strings
 ** representing the input fields as defined by the International Institute
 ** for Applied Systems Analysis WP-87-118, December 1987.
 \*************************************************************************/

int get_fields(char *, char *, char *, char *, char *, char *, char *, char *);

/*************************************************************************\
**  The function title() receives a pointer to the input string and
 ** pointers to the input fields field1 and field2.  The function returns
 ** an integer representing the last field loaded.
 **
 ** Internally, the function parses the input string into two strings
 ** representing the input fields as defined by the International Institute
 ** for Applied Systems Analysis WP-87-118, December 1987.
 \*************************************************************************/

int title(char *, char *, char *);

/**********************************************************************\
**  The function remove_space() removes additional spaces from a
 ** string.  This function will be used after an input string has
 ** been broken into its appropriate
 \**********************************************************************/

void remove_spaces(char *);

#endif /* PARSER_H_ */
