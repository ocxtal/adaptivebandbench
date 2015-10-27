
/**
 * @file util.c
 */

#include <stdlib.h>
#include "log.h"

void *
aligned_malloc(
	size_t size,
	size_t align)
{
	void *ptr = NULL;
	posix_memalign(&ptr, align, size);
	debug("posix_memalign(%p)", ptr);
	return(ptr);
}

/**
 * end of util.c
 */
