#ifndef KMEANS_H
#define KMEANS_H

#include "array.h"

/** Remember to define this variable before including Rmath.h.
 * See the <a href="http://cran.r-project.org/doc/manuals/R-exts.html#Standalone-Mathlib">r-project documentation</a> for more information.
 */
#define MATHLIB_STANDALONE 1
#include <Rmath.h>

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

void kmeans(double **a, unsigned int m, unsigned int n, double **c, unsigned int k, unsigned int *ic1,
	unsigned int *nc, unsigned int iter, double *wss, int *ifault); 

/* 
   a = input dataset of m rows (observations) and n columns
   m = number of observations
   n = number of variables (coordinates)
   c = input initializing k-means and output terminating centers
   k = number of classes
   ic1 = the class id for each observation (upon convergence)
   nc = number of observations in each of k classes
   iter = total number of maximum iterations
   wss = within sums of squares for each group
   ifault = 0 (successful termination)
            1 (empty cluster at some stage)
	    2 (maximum number of iterations reached)
	    3 (number of groups is either less than 2 or greater 
	    than or equal tp the total number of observations)
*/

int srswor(unsigned int n, unsigned int k, unsigned int *y);
int srswor_subset(unsigned int n, unsigned int k, unsigned int *y, unsigned int *z);
/* above are in srswor.c  */
void supkmeans(unsigned int *l, double **a, unsigned int m, unsigned int n, double **c, unsigned int k, unsigned int *ic1,
	unsigned int *nc, unsigned int iter, double *wss, int *ifault); 
void supoptra(unsigned int *l, double **a, unsigned int m, unsigned int n, double **c, unsigned int k, unsigned int *ic1, unsigned int *ic2,
	   unsigned int *nc, double *an1, double *an2, unsigned int *ncp, double *d, unsigned int *itran,
	   unsigned int *live, unsigned int *indx);

void supqtran(unsigned int *l, double **a,unsigned int m,unsigned int n,double **c,unsigned int k,unsigned int *ic1,unsigned int *ic2,unsigned int *nc,
	   double *an1,double *an2,unsigned int *ncp,double *d,unsigned int *itran,unsigned int *indx);


void sampleoptra(unsigned int *l, unsigned int *s,  double **a, unsigned int m, unsigned int n, double **c, unsigned int k, unsigned int *ic1, unsigned int *ic2,
	   unsigned int *nc, double *an1, double *an2, unsigned int *ncp, double *d, unsigned int *itran,
	   unsigned int *live, unsigned int *indx);

void sampleqtran(unsigned int *l, unsigned int *s, double **a,unsigned int m,unsigned int n,double **c,unsigned int k, unsigned int *ic1,unsigned int *ic2,unsigned int *nc,
	   double *an1,double *an2,unsigned int *ncp,double *d,unsigned int *itran,unsigned int *indx);

void samplekmeans(unsigned int *l, unsigned int *s, double **a, unsigned int m, unsigned int n, double **c, unsigned int k, unsigned int *ic1,unsigned int *nc,
	       unsigned int iter,double *wss,int *ifault);




int compare(const void *f1, const void *f2);
double vectormedian(double *vec,int n);
double vectormean(double *vec,int n);
double vecpercentile(double *vec, int n, double p);
void triminput(double **o, double **i, unsigned int h, unsigned int width, int blocksize);
void prepareblockmatrix(double **image, unsigned int h, unsigned int width, int w, double **b, double **m);
void restoremeans(double **b, double **m, double **mean, unsigned int *class, unsigned int num, int w, unsigned int k);
void unrollvectors(unsigned char **x, double **b, unsigned int h, unsigned int num, int w);
void bridgeoutput(unsigned char **x, double **b, double **m, double **mean, unsigned int *class, unsigned int num, unsigned int h, int w, unsigned int k);
unsigned int krza(double *vec, int n, int p);
int update_subset( unsigned int *labeled, unsigned int *sample, unsigned int numblocks, unsigned int k);
unsigned int optimk(unsigned int *label, unsigned int *subset, double **x, unsigned int n, unsigned int p, unsigned int k, double **means, unsigned int *iclass, unsigned int max);
/*above are in blockmeansutils.c  */




#endif /* KMEANS_H */
