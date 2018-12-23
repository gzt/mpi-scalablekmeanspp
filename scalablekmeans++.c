/**
   @file: scalablekmeans++.c

   Scalable kmeans++

   Reference: Bahman Bahmani, Benjamin Moseley, Andrea Vattani, Ravi Kumar, and Sergei Vassilvitskii. 2012. Scalable k-means++. Proc. VLDB Endow. 5, 7 (March 2012), 622-633. DOI=http://dx.doi.org/10.14778/2180912.2180915

   A non-parallel C implementation of k-means|| (or scalable k-means++). Requires wkmeans.c.

   Author:

   Geoffrey Thompson
   Iowa State University
   Department of Statistics
   gzt@iastate.edu

   1/2017

**/

#include "array.h"
#include "util.h"

unsigned int PPSsample(unsigned int m, double *prob);

void wkmeans(double **a, unsigned int m, unsigned int n, double *weights, double **c, unsigned int k, unsigned int *ic1,double *nc,
	     unsigned int iter,double *wss,int *ifault);

/* selected is an indicator for whether the point is chosen this round or not, should be all zeros, length m */
/* returns the count of those selected */
int Multsample(unsigned int m, double *prob, unsigned int *selected){
  int count = 0;
  for(int i = 0; i<m; ++i){
    selected[i] = 0;
    double uni = runif(0., 1.);
    if(uni <= prob[i]){
      count++;
      selected[i] = 1;
    }
  }
  return count;
}


void findneighbors(double **x, unsigned int n, unsigned int p, unsigned int k, double *neighbors, double **centers){
  /* initialize number of neighbors to zero */
  /* chose neighbors to be double rather than int or anything because input weights for wkmeans can be doubles */
  for(int l = 0; l < k; ++l) neighbors[l] = 0.;

  for(int i = 0; i < n; ++i){
    /* for each point, find closest center */
    unsigned int marker = 0;
    double mindist = INFINITY;
    for(int l = 0; l < k; ++l){
      double distance = 0.;
      for(int j = 0; j < p; j++){
	distance += SQ(x[i][j]-centers[l][j]);
	if (distance > mindist) j = p; /*go to the end of the loop*/
      }
      if(distance < mindist){
	mindist = distance;
	marker = l;
      }
    }
    neighbors[marker]++;
  }
}



void wkmeanspp(double **x, unsigned int n, unsigned int p, unsigned int k, double *neighbors, double **centers){
	unsigned int i, j, l, y, *z,  ans , chosen = 0;
	double distMin, sum = 0., *prob, *D ;

	/*  could use weighting as some kind of error code or NA or whatever so not going to error out if find one <0 */

	/* weighted kmeans++ for use in step 8 of scalable kmeans++*/
	if(k == 1){
		/* If k = 1, just choose any center with Probability prop to neighbors[i] */
		for(i = 0; i < k; i++)
			for(j = 0; j < p; j++)
				centers[i][j] = 0.;
		MAKE_VECTOR(prob,n);
		for(i = 0; i < n; i++){
		  if(neighbors[i] <= 0.){
		    prob[i] = 0.;
		  } else {
		    sum += neighbors[i];
		  }
		}
		for(i = 0; i < n; i++){
		  if(neighbors[i] > 0) prob[i] = neighbors[i]/sum;
		}
		/* y = n * runif(0., 1.); */
		y = PPSsample(n, prob);
		for (i = 0;i < p; i++)
			centers[0][i] = x[y][i];
		FREE_VECTOR(prob);
	} else {
		for(i = 0; i < k; i++)
			for(j = 0; j < p; j++)
				centers[i][j] = 0.;
		MAKE_VECTOR(prob, n);
		MAKE_VECTOR(z, k);
		MAKE_VECTOR(D, n);
		sum = 0.;
		/* choose first center at random with prob proportional to neighbors[i] */
		for(i = 0; i < n; i++){
		  if(neighbors[i] <= 0.){
		    prob[i] = 0.;
		  } else {
		    sum += neighbors[i];
		  }
		}
		for(i = 0; i < n; i++){
		  prob[i] = 0.;
		  if(neighbors[i] > 0) prob[i] = neighbors[i]/sum;
		}
		/* y = n * runif(0., 1.); */
		y = PPSsample(n, prob);
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

			  distMin = 0.;
			  for(j = 0; j < p; j++){
					distMin += SQ(x[i][j]-centers[l][j]);
					if (distMin > D[i]) j = p; /*go to the end of the loop*/
			  }
					if(distMin < D[i])
						D[i] = distMin;
			
			/* Compute denominator */
					if(neighbors[i] > 0.)
					  sum += D[i]*neighbors[i]; 
			/* Compute probability note is define as, Prob = n_i ||x_i - bmu||^2/\sum n_i ||x_i -\bmu ||^2 */
			}
			/* Set those already chosen Distance = 0. since we don't want to repeat centers */
			for(i = 0; i < l; i++)
			  D[z[i]] = 0.0;
			/* This should guarantee Prob = 0 for those points */
			for(i = 0; i < n; i++)
			  prob[i] = (neighbors[i] > 0.) ? D[i]*neighbors[i]/sum : 0.;

      

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


/* not done yet  */

void scalablekmeanspp(double **x, unsigned int n, unsigned int p, unsigned int k, unsigned int rounds, double L, double **centers){
  /*
    x : data matrix, n x p
    n : number of rows of x
    p : number of columns of x
    k : number of clusters
    rounds : number of rounds of kmeans|| - paper suggests 5 is sufficient
    L : oversampling factor, should be a constant times k, paper chose between .5k and 2k
    centers: the cluster centers that get returned, KxP matrix
   */
	unsigned int i, j, y, *z;
	unsigned int *selected;


	/*
	  selected: length N vector of 0 and 1 indicating whether the point has been selected this round
	 */
	
	double **chosenpoints;
	double sum = 0., *prob, *D ;
	
	/* Strategy for storing point selections and the points themselves:
	   Make a vector and matrix large enough to store all the points (with high probability)
	   But hopefully not so large that it causes problems.
	   Insert error checking and cap if it ends up too small
	   This is easier/safer than realloc()-ing each round
	   Just have to be sure not to read out into uninitialized memory!
	   For large N, could use poisson approx to poisson-binomial distribution for tighter bound
	   see http://projecteuclid.org/download/pdf_1/euclid.aoms/1177705799
	   though that is overkill
	*/
	
	int pointsbound =  ((rounds * L * 5)/4 + 50);

	// debugging output: 
	// printf("pointsbound: %d \n", pointsbound);
	
	if(pointsbound > n) pointsbound = n;
	/* the multiplicative factor is overkill for decent N, rounds, and L 
	   but I want to make sure works even for small N, rounds, or L (hence adding the constant) */


		for(i = 0; i < k; i++)
			for(j = 0; j < p; j++)
				centers[i][j] = 0.;
		MAKE_VECTOR(prob, n);
		MAKE_VECTOR(D, n);
		/* D is vector of smallest distance */
		MAKE_VECTOR(z, pointsbound);
		/* z is the vector of point selections */
		MAKE_MATRIX(chosenpoints,pointsbound,p);
		/* need to check that don't go past the bounds */

		
		/* choose first center at random with Prob 1/n */
		y = n * runif(0., 1.);
		for(i = 0; i < p; i++)
			chosenpoints[0][i] = x[y][i];

		z[0] = y; /* store index of observation chosen */


		// assign initial probabilities
		// calc initial distances
		sum = 0.;
   
		for(i = 0; i < n; i++){
		  prob[i] = 0.;
		  D[i] = INFINITY;
		  }

		  /* Compute probability note is defined as, Prob = L ||x_i - bmu||^2/\sum  ||x_i -\bmu ||^2 */
		  

			/* Now do steps 3:6 repeated rounds of sampling */
			MAKE_VECTOR(selected,n);
			  int centercount = 0;
			  int newcenters = 1;
			  double currentDist = 0.;
			  int currentpoint = 0;
			   /* doing the rounds */

			  for(int step = 0; step < rounds; ++step){

			    // debugging:
			    //printf("centercount: %d new centers: %d round: %d \n",centercount, newcenters,step);

			    

			      sum = 0.;
			      /* do distance calculation */
			      /* use old centercount and newcenters to get the right part of the matrix */
			      if(newcenters+centercount > pointsbound){
				 /* if we pick too many points, just cap it off */
				 /* this is not good behavior */
				newcenters = pointsbound - centercount;
				printf("ERROR: too many points at beginning of round %d\n",step);
				rounds = step;
			      }
			    
			      for(i = 0; i < n; i++){
				for(currentpoint = centercount; currentpoint < (centercount + newcenters) && currentpoint < pointsbound; currentpoint++){
				   /* look only at rows centercount:(centercount + newcenters) in chosenpoints[][] */
				  currentDist = 0.;
				  for(j = 0;j<p;++j){
				    currentDist += SQ(x[i][j]-chosenpoints[currentpoint][j]);
				    if(currentDist > D[i]) j = p;
				  }
				  if(currentDist < D[i]) D[i] = currentDist;
				}
				sum += D[i];
			      }
			    
			    
			     /* set centercount to the current number of used rows in chosenpoints[][] */
			    centercount = centercount + newcenters;
			    newcenters = 0;
			    // debugging:
			    //printf("centercount: %d new centers: %d round: %d \n",centercount, newcenters,step);

			    
			    /* compute probabilities */
			    for(i = 0; i < n; i++)
			      prob[i] = (1.0*L) * D[i]/sum;

			     /* sample new centers */
			    newcenters = Multsample(n, prob, selected);

			    // debugging output
			    // printf("%d new centers in round %d \n", newcenters,step);
			    
			     /* error handling, should be unnecessary */


			    if(newcenters+centercount > pointsbound){
			       /* if we pick too many points, just cap it off */
			       /* this is not good behavior */
			      newcenters = pointsbound - centercount;
			      printf("ERROR: too many points chosen in round %d\n",step);
			      rounds = step;								
			    }

			    
			     
			      unsigned int tmpstorage=0;
			      for(i = 0; i<n && tmpstorage < newcenters; ++i){
				if(selected[i]==1){
				  D[i] = 0.;
				  z[centercount+tmpstorage] = i;
				  for(j = 0; j < p; ++j)
				    chosenpoints[centercount+tmpstorage][j] =  x[i][j];
				  tmpstorage++;
				}
			      }	      
			  

			    } /* end rounds loop */
			    centercount = centercount + newcenters;

			    // debugging output
			    printf("Final number of centers: %d\n", centercount);
			    
		FREE_VECTOR(z);
		FREE_VECTOR(prob);
		FREE_VECTOR(D);
		FREE_VECTOR(selected);

		/* move so we don't have problems */
		if(centercount > pointsbound){
		  printf("ERROR: too many points chosen for bounds \n");
		  FREE_MATRIX(chosenpoints);
		  return;
		} else {
		  double **tempcenters;
		  MAKE_MATRIX(tempcenters, centercount, p);
		  for(i = 0; i < centercount; ++i)
		    for(j = 0; j < p; ++j) tempcenters[i][j] = chosenpoints[i][j];
		  FREE_MATRIX(chosenpoints);

		
		
		  double  *neighbors;
		  /* neighbors: vector of number of neighbors of the chosen centers - double because goes into wkmeans */
		  MAKE_VECTOR(neighbors,centercount);	
		/* find the neighbors, then run weighted kmeans++ to initialize, then run weighted kmeans  */
		findneighbors( x, n, p, centercount, neighbors, tempcenters);
		
		
		wkmeanspp(tempcenters,centercount,p,k,neighbors, centers);
		int ifault = 0;
		unsigned int *ic1;
	        double *nc;
		double *wss;
		MAKE_VECTOR(ic1,centercount);
		MAKE_VECTOR(nc,k);
		MAKE_VECTOR(wss,k);

		wkmeans(tempcenters, centercount, p, neighbors, centers, k, ic1,nc,1000,wss,&ifault);
		if(ifault != 0) {
		  printf("ERROR: weighted k-means failed to converge, IFAULT = %d \n",ifault);
		}
		 /* FREE everything at this point. CHECK IF DONE */
		FREE_VECTOR(ic1);
		FREE_VECTOR(nc);
		FREE_VECTOR(wss);
		
		
		FREE_MATRIX(tempcenters);
		
		FREE_VECTOR(neighbors);
		return;
		}
}

