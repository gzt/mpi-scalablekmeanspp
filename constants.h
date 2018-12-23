#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

/* suss out C standard enforced */
#if defined(__STDC__)
#	define C89
#	if defined(__STDC_VERSION__)
#		define C90
#		if (__STDC_VERSION__ >= 199409L)
#			define C94
#		endif
#		if (__STDC_VERSION__ >= 199901L)
#			define C99
#		endif
#		if (__STDC_VERSION__ >= 201112L)
#			define C11
#		endif
#	endif
#endif

/*
#define SIZE_T int
#define SIZE_T_FMT "%d"
*/
#include <stddef.h>
#define SIZE_T size_t
#define SIZE_T_FMT "%zu"

#endif
