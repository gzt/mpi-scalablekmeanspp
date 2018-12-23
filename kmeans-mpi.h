/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         kmeans-mpi.h   (an OpenMP version)                            */
/*   Description:  header file for a simple k-means clustering program       */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Copyright (c) 2005 Wei-keng Liao
// Copyright (c) 2017 Geoffrey Z Thompson
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// -----------------------------------------------------------------------------

#ifndef _H_KMEANSMPI
#define _H_KMEANSMPI

#include "array.h"
#include "util.h"
/** Remember to define this variable before including Rmath.h.
 * See the <a href="http://cran.r-project.org/doc/manuals/R-exts.html#Standalone-Mathlib">r-project documentation</a> for more information.
 */
#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include <assert.h>

#define msg(format, ...) do { fprintf(stderr, format, ##__VA_ARGS__); } while (0)
#define err(format, ...) do { fprintf(stderr, format, ##__VA_ARGS__); exit(1); } while (0)

#define malloc2D(name, xDim, yDim, type) do {               \
    name = (type **)malloc(xDim * sizeof(type *));          \
    assert(name != NULL);                                   \
    name[0] = (type *)malloc(xDim * yDim * sizeof(type));   \
    assert(name[0] != NULL);                                \
    for (size_t i = 1; i < xDim; i++)                       \
        name[i] = name[i-1] + yDim;                         \
} while (0)

#ifdef __CUDACC__
inline void checkCuda(cudaError_t e) {
    if (e != cudaSuccess) {
        // cudaGetErrorString() isn't always very helpful. Look up the error
        // number in the cudaError enum in driver_types.h in the CUDA includes
        // directory for a better explanation.
        err("CUDA Error %d: %s\n", e, cudaGetErrorString(e));
    }
}

inline void checkLastCudaError() {
    checkCuda(cudaGetLastError());
}
#endif

double** omp_kmeans(int, double**, int, int, int, double, int*);
double** seq_kmeans(double**, int, int, int, double, int*, int*);
double** cuda_kmeans(double**, int, int, int, double, int*, int*);

double** file_read(int, char*, int*, int*);
int     file_write(char*, int, int, int, double**, int*);


double  wtime(void);

extern int _debug;


/** Default maximum number of K-means iterations.
 * The user can specify any number of maximum iterations when calling 
 * #kmeans(), but for consistency across the code, this define can be 
 * used.
 */
#define MAX_KMEANS_ITERATIONS	1000

/** Errors producible by #kmeans().  */
enum {
	KMEANS_NO_ERROR,		/*!< error-free condition */
	KMEANS_NULL_CLUSTER_ERROR,	/*!< K-means produced null cluster */
	KMEANS_EXCEED_ITERATIONS_ERROR,	/*!< K-means exceeded max iterations */
	KMEANS_CALLER_INPUT_ERROR,	/*!< invalid caller input */
	KMEANS_UNUSED_ERROR,		/*!< not used */
	KMEANS_NO_SEEDS,		/*!< initialization gives 0 seeds */
	KMEANS_NUMBER_ERRORS		/*!< number of K-means errors */
};

/** Human-friendly character strings describing each K-means error.
 */
extern const char *KMEANS_ERROR_STRING[];


#endif
