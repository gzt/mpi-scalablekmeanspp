/**
   @file: kmeans++.c
   
   Perform kmeans++ algorithm.

   Reference: Arthur, David, and Sergi Vassilvitskii. "K-means++: The Advantages of Careful Seeding." SODA ‘07: Proceedings of the Eighteenth Annual ACM-SIAM Symposium on Discrete Algorithms. 2007, pp. 1027–1035.

   Author:

   Israel Almodovar-Rivera and Ranjan Maitra
   Iowa State University
   Department of Statistics
   israel.almodovar@upr.edu and maitra@iastate.edu

   Copyright 2015-2016

   edited by Geoffrey Thompson to reduce memory and number of loops 12/30/16
   Iowa State University
   Department of Statistics
   gzt@iastate.edu

**/

#include "array.h"
#include "util.h"

unsigned int PPSsample(unsigned int m, double *prob)
{
  double uni = runif(0., 1.);
  // double *cumprob;
  double cumprob;
	unsigned int i;
	/* this could probably be rewritten to not require a cumprob vector */
	/* MAKE_VECTOR(cumprob, m); */
	/* cumprob[0] = prob[0]; */

	/* for (i = 1; i < m; i++)  */
	/* 	cumprob[i] = cumprob[i - 1] + prob[i]; */

	
	cumprob = 0.;
	for ( i = 0; ( (i < m) && ((cumprob += prob[i]) <= uni) ) ; i++);
	/* for (i = 0; ( (i < m) && (cumprob[i] <= uni)); i++); */
	/* FREE_VECTOR(cumprob); */
	return i;
}

void kmeanspp(double **x, unsigned int n, unsigned int p, unsigned int k, double **centers){
	unsigned int i, j, l, y, *z,  ans , chosen = 0;
	double distMin, sum = 0., *prob, *D ;
  
	if(k == 1){
		/* If k = 1, just choose any center with Probability =1/n */
		for(i = 0; i < k; i++)
			for(j = 0; j < p; j++)
				centers[i][j] = 0.;
		y = n * runif(0., 1.);
		for (i = 0;i < p; i++)
			centers[0][i] = x[y][i];
	} else {
		for(i = 0; i < k; i++)
			for(j = 0; j < p; j++)
				centers[i][j] = 0.;
		MAKE_VECTOR(prob, n);
		/* rewriting this to have distMin as scalar - don't need to store a matrix! */
		/* keeping the name, though */
		MAKE_VECTOR(z, k);
		MAKE_VECTOR(D, n);
    
		/* choose first center at random with Prob = 1/n */
		y = n * runif(0., 1.);
		for(i = 0; i < p; i++)
			centers[0][i] = x[y][i];
    
		z[0] = y; /* store index of observation chosen */
		/* Set distance to INFINITY to find minimum */
		for(i = 0; i < n; i++)
			D[i] = INFINITY;  
    
		for(l = 0; l < (k-1); l++){
			chosen = l + 1;
			sum = 0.;
			for(i = 0; i < n; i++){
				prob[i] = 0.;

			}
			for(i = 0; i < n; i++){
			   /* collapsing a few loops into one loop! */
			  distMin = 0.;
			  for(j = 0; j < p; j++){
					distMin += SQ(x[i][j]-centers[l][j]);
					if (distMin > D[i]) j = p; /*go to the end of the loop*/
			  }
					if(distMin < D[i])
						D[i] = distMin;
			
			/* Compute denominator */
				sum += D[i]; 
			/* Compute probability note is define as, Prob = ||x_i - bmu||^2/\sum ||x_i -\bmu ||^2 */
			}
			/* Set those already chosen Distance = 0. since we don't want to repeat centers */
			for(i = 0; i < l; i++)
			  D[z[i]] = 0.0;
			/* This should guarantee Prob = 0 for those points */
			for(i = 0; i < n; i++)
			  prob[i] = D[i]/sum;

      

			/* Choose one center with probability defined from above */
      
			ans = PPSsample(n, prob);
			z[chosen] = ans;
			D[ans] = 0.;
			for(j = 0; j < p; j++)
				centers[chosen][j] = x[ans][j];
      
		}
		FREE_VECTOR(z);
		FREE_VECTOR(prob);
		FREE_VECTOR(D);
	}
}


