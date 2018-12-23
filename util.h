#ifndef _UTIL_H
#define _UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include "constants.h"

/**
 * Simple mathematical functions.
 */
#define SQ(x) ((x) * (x)) 
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQ(x) ((x) * (x))
#define LTINDEX(k, l) ((MAX(k, l) * (MAX(k, l) + 1))/2 + MIN(k, l))
#define LLTINDEX(k, l) (LTINDEX(k - 1, l))
enum {NO_ERROR,
	MEMORY_ALLOCATION_ERROR,	/*!< failure to allocate memory */
	INVALID_COMMAND_LINE,		/*!< user misuse of command line */
	FILE_OPEN_ERROR,		/*!< failure to open file */
	FILE_FORMAT_ERROR,		/*!< unexpected file format */
	FILE_WRITE_ERROR,		/*!< failure to write to file */
	INTERNAL_ERROR,			/*!< any other error */
	NUMBER_OF_ERRORS
};

#define NO_EXIT 0	/* do not let called function exit */
#define DO_EXIT 1	/* induce a called function to exit */

extern const char *error_string[];


/**
 * Convert upper triangle symmetric matrix indices to vector index.
 *
 * @param n size of matrix
 * @param i row index
 * @param j index, must have j>i
 * @return index of same position in row-order vector
 */
#define VINDEX(n, i, j) (                                     \
	(i) < (j)                                             \
		? (j) + (i) * (n) - ((i) + 1) * ((i) + 2) / 2 \
		: (i) + (j) * (n) - ((j) + 1) * ((j) + 2) / 2 \
)


/**
 * Convert vector index to upper triangle indices.
 *
 * @param n size of matrix
 * @param i row index
 * @param j column index
 * @param ij vector index
 */
#define RINDEX(n, i, j, ij) do {                \
	unsigned int NADD = n - 1, NINDEX = 0;           \
	for (i = 0;                             \
		ij >= NINDEX + NADD;             \
		i++, NINDEX += NADD, NADD--);   \
	j = i + ij - NINDEX + 1;                       \
} while (0)

/**
 * Generic matrix object.  When you need to operate on a matrix using generic
 * functions, like random samplers, that use call-backs to check user-required
 * conditions, can pass this object around as a void pointer.  It currently
 * just provides the data and the number of columns.
 */
struct data_obj {
	double **x;	/*!< data matrix */
	unsigned int p;	/*!< number of columns in matrix */
};



#endif
