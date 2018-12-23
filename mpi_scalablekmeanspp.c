/**
   @file: mpi_scalablekmeanspp.c

   Scalable kmeans++

   Reference: Bahman Bahmani, Benjamin Moseley, Andrea Vattani, Ravi Kumar, and Sergei Vassilvitskii. 2012. Scalable k-means++. Proc. VLDB Endow. 5, 7 (March 2012), 622-633. DOI=http://dx.doi.org/10.14778/2180912.2180915
   Author:

   Geoffrey Thompson
   Iowa State University
   Department of Statistics
   gzt@iastate.edu or gzthompson@gmail.com

   1/2017

**/

#include <mpi.h>
#include "array.h"
#include "util.h"
#include <Rmath.h>
#include <assert.h>


/* NOTE: need to fix number generation strategy - parallelization is a problem! 
   Presumes caller has initiated correctly! */


unsigned int PPSsample(unsigned int m, double *prob);

void wkmeans(double **a, unsigned int m, unsigned int n, double *weights, double **c, unsigned int k, unsigned int *ic1,double *nc,
	     unsigned int iter,double *wss,int *ifault);

void wkmeanspp(double **x, unsigned int n, unsigned int p, unsigned int k, double *neighbors, double **centers);
int Multsample(unsigned int m, double *prob, unsigned int *selected);
void findneighbors(double **x, unsigned int n, unsigned int p, unsigned int k, double *neighbors, double **centers);

unsigned int PPSsample(unsigned int m, double *prob)
{
  double uni = runif(0., 1.);
  double cumprob;
	unsigned int i;
	cumprob = 0.;
	for ( i = 0; ( (i < m) && ((cumprob += prob[i]) <= uni) ) ; i++);
	return i;
}

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




void mpi_findneighbors(double **x, unsigned int n, unsigned int p, unsigned int k, double *neighbors, double **centers, int root, MPI_Comm comm){
  /* initialize number of neighbors to zero */
  /* chose neighbors to be double rather than int or anything because input weights for wkmeans can be doubles */
  /* entering the function, everybody is presumed to have the exact same copy of "centers" */
  double *localneighbors;
  extern int _debug;
  int rank;
  MPI_Comm_rank(comm,&rank);
  MAKE_VECTOR(localneighbors,k);
  for(int l = 0; l < k; ++l) localneighbors[l] = 0.;

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
      if(distance <= mindist){
	mindist = distance;
	marker = l;
      }
    }
    localneighbors[marker]++;
  }
  // if(_debug) printf("Center %d has %1.0f neighbors on process %d\n", 0,localneighbors[0],rank);
  MPI_Allreduce(localneighbors,neighbors,k,MPI_DOUBLE,MPI_SUM,comm);
   for(int i = 0; i < k; ++i)
    if(neighbors[i]<1) neighbors[i] = 1;
   // if(_debug && rank == root) printf("Center 0 has %1.0f neighbors\n",neighbors[0]);
  FREE_VECTOR(localneighbors);
  //  if(_debug && rank == root)    for(int i = 0; i < k; ++i) printf("%1.0f ",neighbors[i]);
  
  return;
}



/* stuff to do: really need to have some robust error handling in here 
   - there are a lot of things that can go wrong */




void mpi_scalablekmeanspp(double **x, unsigned int n, unsigned int p, unsigned int k, unsigned int rounds, double L, double **centers, int root, MPI_Comm comm){

  /*
    x: data matrix, n x p
    n: number of rows of x
    p: number of columns of x
    k: number of clusters
    rounds: number of rounds of kmeans|| - paper suggests 5 is sufficient
    L: oversampling factor, should be a constant times k, paper chose between .5k and 2k
    centers: the cluster centers that get returned, KxP matrix
   */
  int rank, nproc;
  unsigned int i, j, y, *localz; // don't actually need z
	unsigned int *selected;
	/*
	  selected: length N vector of 0 and 1 indicating whether the point has been selected this round
	 */
	extern int _debug;
	MPI_Comm_rank(comm,&rank);
	MPI_Comm_size(comm,&nproc);

	/* zero-ing out the output matrix */
	for(i = 0; i < k; i++)
	  for(j = 0; j < p; j++)
	    centers[i][j] = 0.;

	
	double **chosenpoints;
	double **localpoints;
	double sum = 0., localsum=0., *prob, *D ;
	
	/* Strategy for storing point selections and the points themselves:
	   Make a vector and matrix large enough to store all the points (with high probability)
	   But hopefully not so large that it causes problems.
	   Insert error checking and cap if it ends up too small
	   This is easier/safer than realloc()-ing each round
	   Just have to be careful not to read out into uninitialized memory!
	   For large N, could use poisson approx to poisson-binomial distribution for tighter bound
	   see http://projecteuclid.org/download/pdf_1/euclid.aoms/1177705799
	   though that is overkill
	*/
	
	int pointsbound =  ((rounds * L * 5)/4 + 50);
	// if(pointsbound > nproc*n) pointsbound = nproc*n;
	if(_debug && rank == root) printf("Points bound = %d \n", pointsbound);

	chosenpoints    = (double**) malloc(pointsbound *             sizeof(double*));
	assert(chosenpoints != NULL);
	chosenpoints[0] = (double*)  malloc(pointsbound * p * sizeof(double));
	assert(chosenpoints[0] != NULL);
	for (i=1; i<pointsbound; i++)
	  chosenpoints[i] = chosenpoints[i-1] + p;
	
	int localbound = (L * 2)/(nproc)+50;
	// if(_debug) printf("Rank = %d local bound = %d\n",rank,localbound);
	/* localbound is more generous because we can't assume uniform distribution among processors */
	
	// if(localbound > n) localbound = n;
	// not doing this limit because n varies by processor and may have to send something of size localbound!
	/* the multiplicative factor for pointsbound is overkill for decent N, rounds, and L 
	   but I want to make sure works even for small N, rounds, or L (hence adding the constant) */


	localpoints    = (double**) malloc(localbound *             sizeof(double*));
	assert(localpoints != NULL);
	localpoints[0] = (double*)  malloc(localbound * p * sizeof(double));
	assert(localpoints[0] != NULL);
	for (i=1; i<localbound; i++)
	  localpoints[i] = localpoints[i-1] + p;


	/*
	  Strategy:
	  0. each processor has n points (n varies by processor).  
	  1. in master node, select at random one point to initialize. *Should* do random from all.
	  2. communicate point to all 	  
	  3. put point in "chosenpoints" matrix - every processor has a copy 
	  4. compute distances 
	  5. summarize and broadcast total 
	  6. sample points 
	  7. put points in "localpoints" matrix - bigger than needed, everybody has one of roughly same size 
	  8. send number of points chosen, localpoints matrix back to master node
	  9. in master node, put localpoints into chosenpoints, then send chosenpoints, total # selected to all 
	  10. do 5-9 for R rounds
	  11. find neighbors of points, send back to master
	  12. on master node, finish it off as if it were all on one processor.

	 */
	  MAKE_VECTOR(localz,localbound);
	  	
	  /* vector of locally-selected means, written over each round */
	  /* Keeping it just in case it comes in handy */
	  
	/* Step 1: select initial point */
	if(rank == root){
	  //	  MAKE_VECTOR(z, pointsbound); // don't actually need z vector
	  /* choose first center at random with Prob 1/n */
	  y = n * runif(0., 1.);
	  for(i = 0; i < p; i++)
	    chosenpoints[0][i] = x[y][i];
	  localz[0] = y; /* store index of observation chosen */
	}

	
	  /* step 4 */
	
	  MAKE_VECTOR(D, n);
	  /* D is vector of smallest distance */
	  MAKE_VECTOR(prob, n);
	  localsum = 0.;
	  
	  for(i = 0; i < n; i++){
	    prob[i] = 0.;
	    D[i] = INFINITY;
	    }
	  
	  /* Step 5 */
	  MAKE_VECTOR(selected,n);
	  unsigned int centercount = 0;
	  unsigned int localcenters = 0;
	  unsigned int localnewcenters=0;
	  unsigned int newcenters = 1;
	  double currentDist = 0.;
	  unsigned int currentpoint = 0;
	  int* gathersample ;
	  MAKE_VECTOR(gathersample,nproc);
	  for(i = 0; i < nproc; ++i) gathersample[i]=0;
	  int* cumcount ;
	  MAKE_VECTOR(cumcount,nproc);
	  for(i = 0; i < nproc; ++i) cumcount[i]=0;

	  /* Enter the loop of rounds here: */
	  for(int step = 0; step < rounds; ++step){
	    	/* step 2/3 */
    			   
			      localsum = 0.;
			      /* do distance calculation */
			      /* use old centercount and newcenters to get the right part of the matrix */
			      if(newcenters+centercount > pointsbound){
				 /* if we pick too many points, just cap it off */
				 /* this is not good behavior */
				newcenters = pointsbound - centercount;
				if(rank == root)
				  printf("ERROR: too many points at beginning of round %d\n",step);
				
				rounds = step;
			      }
			      /* trying to get away with sending less info! 
			         Safer would be to send whole matrix or from [0][0] */
			      if(_debug && rank == root) printf("centercount: %d newcenters: %d step: %d",centercount,newcenters,step);
			      MPI_Bcast(chosenpoints[centercount], (newcenters)*p,MPI_DOUBLE,root,comm);
			      // if(_debug && rank == root) printf("Successful broadcast in step %d\n",step);

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
				localsum += D[i];
			      }
			   
			    MPI_Allreduce(&localsum,&sum,1,MPI_DOUBLE,MPI_SUM,comm);	  
			    	if(_debug && rank == root) printf("Sum = %f \n", sum);

			    centercount = centercount + newcenters;
			    newcenters = 0;
			    // debugging:
			    	if(_debug && rank == root) 
				  printf("centercount: %d new centers: %d round: %d \n",centercount, newcenters,step);

			    /* Step 6 */

			    /* compute probabilities */
			    for(i = 0; i < n; i++)
			      prob[i] = (1.0*L) * D[i]/sum;

			     /* sample new centers */
			    localnewcenters = Multsample(n, prob, selected);
			    // if(_debug && rank == root) printf("Successfully sampled %d points in %d \n", localnewcenters, rank);

			    if(localnewcenters > localbound) localnewcenters = localbound;
			    /* capping off number chosen with localbound - this is bad behavior */
			    /* should throw an error - will recode later */

			    /* Step 7 */
			      unsigned int tmpstorage=0;
			      for(i = 0; i<n && tmpstorage < localnewcenters; ++i){
				if(selected[i]==1){
				  D[i] = 0.;
				  localz[tmpstorage] = i;
				  /* why am I storing this if it gets overwritten each round? */
				  for(j = 0; j < p; ++j)
				    localpoints[tmpstorage][j] =  x[i][j];
				  tmpstorage++;
				}
			      }

			      /* Step 8 */
			      /* Step 9 */
			      /* putting tempcenters into chosenpoints */
	    if(rank == root){
	      cumcount[0] = 0;
	      /* this indicates the starting point of the entries for processor i */
	      gathersample[0] = localnewcenters;
	      /* this vector stores the number of new points per processor */

	      // if(_debug) printf("cumulative count: %d \n", cumcount[nproc-1]);
	      for(i = 0; i < localnewcenters; ++i){
		for(j = 0; j < p; ++j)
		  chosenpoints[centercount+i][j] = localpoints[i][j];
	      }
	      /* writing to assume master node is 0 - not sure how to avoid */
	      for(int node = 1; node < nproc; node++){
		//if(_debug) printf("gather sample = %d node %d \n", gathersample[node], node);
		MPI_Status Message;
		int number_amount;
		// if(_debug) printf("get info node %d ",node);
		MPI_Probe(node,step,comm,&Message);
		MPI_Get_count(&Message, MPI_DOUBLE, &number_amount);
		// if(_debug) printf("number: %d\n", number_amount);
		gathersample[node] = number_amount/p;
		cumcount[node] = cumcount[node-1] + gathersample[node-1];
		MPI_Recv(localpoints[0], number_amount, MPI_DOUBLE, node,step , comm,MPI_STATUS_IGNORE);
		for(i = 0; i < gathersample[node] && (centercount + cumcount[node]+i < pointsbound); ++i){
		  for(j = 0; j < p; ++j)
		    chosenpoints[centercount + cumcount[node]+i][j] = localpoints[i][j];
		}
	      }
	      newcenters = 0;
	      for(int node = 0; node < nproc; ++node) newcenters += gathersample[node];
	    } else {
	      MPI_Send(localpoints[0], localnewcenters*p, MPI_DOUBLE, root,step,comm);
	      // printf("node %d \n",rank);
	      // if(_debug) printf("Node: %d. %d points sent.\n", rank,localnewcenters);
	    }
	    MPI_Bcast(&newcenters,1,MPI_INT,root,comm); 
	  } /* Step 10: do for R rounds */

	  MPI_Bcast(chosenpoints[centercount], (newcenters)*p,MPI_DOUBLE,root,comm);
	  centercount = centercount + newcenters;
	  
	  // debugging output
	  if(_debug && rank == root)
	    printf("Final number of centers: %d\n", centercount);

	  if(centercount > pointsbound ){
	    
	    if(rank == root && _debug)
	      printf("ERROR: too many points chosen for bounds \n");
	    
	    free(chosenpoints[0]);
	    free(chosenpoints);
	    FREE_VECTOR(selected);
	    FREE_VECTOR(prob);
	    FREE_VECTOR(D);
	    FREE_VECTOR(localz);
	    FREE_VECTOR(gathersample);
	    FREE_VECTOR(cumcount);
	    free(localpoints[0]);
	    free(localpoints);
	    /* Should probably do more gracefully in MPI than in single processor... */
	    return;
	  } else {
	    double **tempcenters;
    
	    tempcenters    = (double**) malloc(centercount *             sizeof(double*));
	    assert(tempcenters != NULL);
	    tempcenters[0] = (double*)  malloc(centercount * p * sizeof(double));
	    assert(tempcenters[0] != NULL);
	    for (i=1; i<centercount; i++)
	      tempcenters[i] = tempcenters[i-1] + p;

	    
	    for(i = 0; i < centercount; ++i)
	      for(j = 0; j < p; ++j) tempcenters[i][j] = chosenpoints[i][j];
	    // if(_debug) printf("Node %d point %f\n",rank,chosenpoints[centercount-1][p-1]);
	    if(_debug && rank == root){
	      printf("Final centers, center count = %d rank %d ",centercount,rank);
	      //for(i = 0; i < centercount; ++i) printf(" %.2f",tempcenters[i][p-1]);
	      printf("\n");
	    }
	     

	    free(chosenpoints[0]);
	    free(chosenpoints);
	    
	  
	  /* Step 11: Find neighbors */
	  double  *neighbors;
	  /* neighbors: vector of number of neighbors of the chosen centers - double because goes into wkmeans */
	  MAKE_VECTOR(neighbors,centercount);	
	  /* find the neighbors, then run weighted kmeans++ to initialize, then run weighted kmeans  */
	

	  FREE_VECTOR(selected);
	  FREE_VECTOR(prob);
	  FREE_VECTOR(D);
	  FREE_VECTOR(localz);
	  free(localpoints[0]);
	  free(localpoints);
	  /* entering this function, everybody needs to agree on tempcenters contents */

	  if(_debug && rank == root) printf("Finding neighbors \n");
	  // MPI_Barrier(comm);
	  mpi_findneighbors(x, n, p, centercount, neighbors, tempcenters, root, comm);
	  if(_debug && rank == root) printf("Neighbors found \n");
	  if(rank == root){
	    if(_debug) printf("entering weighted kmeans++\n");
		wkmeanspp(tempcenters,centercount,p,k,neighbors, centers);
		int ifault = 0;
		unsigned int *ic1;
	        double *nc;
		double *wss;
		MAKE_VECTOR(ic1,centercount);
		MAKE_VECTOR(nc,k);
		MAKE_VECTOR(wss,k);
		if(_debug) printf("entering weighted kmeans \n");
		wkmeans(tempcenters, centercount, p, neighbors, centers, k, ic1,nc,1000,wss,&ifault);
		if(ifault != 0) {
		  printf("ERROR: weighted k-means failed to converge, IFAULT = %d \n",ifault);
		}
		 /* FREE everything at this point. CHECK IF DONE */
		FREE_VECTOR(ic1);
		FREE_VECTOR(nc);
		FREE_VECTOR(wss);
		if(_debug) printf("done with weighted kmeans! \n");
	  }
	  
	  free(tempcenters[0]);
	  free(tempcenters);
	  FREE_VECTOR(neighbors);
	  

	  FREE_VECTOR(gathersample);
	  FREE_VECTOR(cumcount);

	 
	  return;
	  }
}
