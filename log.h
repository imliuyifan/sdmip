/*
 * log.h
 *
 *  Created on: Jun 10, 2013
 *      Author: lian
 */

#ifndef LOG_H_
#define LOG_H_
#include "sdglobal.h"

BOOL malloc_verify(void);
void *log_alloc(char *string, void *return_ptr, int size);
void *log_free(char *string, void *ptr);
void *log_realloc(char *string, void *free_ptr, void *alloc_ptr, int size);
void log_start(sdglobal_type* sd_global);
void log_stop(sdglobal_type* sd_global);

#endif /* LOG_H_ */
