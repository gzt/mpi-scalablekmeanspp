#define BIG 1e+40
#include "array.h"



/* Code based on Applied Statistics algorithms (C) Royal Statistical Society
   1979. Adapted for C by Ranjan Maitra, Baltimore, 07/12/02 maitra@iastate.edu  */
/* Adapted into weighted k-means by Geoffrey Thompson, Ames, IA 12/31/16  gzthompson@gmail.com         */

/* This currently only works with strictly positive weights. */
/* Could do a lloyd's implementation */

void woptra(double **a, unsigned int m, unsigned int n, double *weights, double **c, unsigned int k, unsigned int *ic1, unsigned int *ic2,
	   double *nc, unsigned int *ncp, double *d, unsigned int *itran,
	   unsigned int *live, unsigned int *indx);

void wqtran(double **a,unsigned int m,unsigned int n,double *weights, double **c,unsigned int k,unsigned int *ic1,unsigned int *ic2,double *nc,
	   unsigned int *ncp,double *d,unsigned int *itran,unsigned int *indx);

void wkmeans(double **a, unsigned int m, unsigned int n, double *weights, double **c, unsigned int k, unsigned int *ic1,double *nc,
	   unsigned int iter,double *wss,int *ifault)
{
	unsigned int i, j, l, ij, il, indx, *ic2, *ncp, *itran, *live;
	double temp, da, db, dc, dt[2], *d;
	
	/* ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1 
	   Divide M prints in N-dimensional space unsigned into K clusters so that 
	   the within cluster sum of squares is minimized. */ 
	
	MAKE_VECTOR(ic2, m);
	
	*ifault = 3;
	if ((k <=1) || (k >= m)) {  /* check for number of clusters */
		FREE_VECTOR(ic2);
		return;
	}
	*ifault = 0;

	/* For each print i, find its two closest centres, ic1(i) and ic2(i).
	   Assign it to ic1(i). */
	
	double minweight = BIG;
	/* used for finding empty clusters */
	
	for (i = 0; i < m; i++) {
	  /* for this bit, the weights don't really matter */
	  if(weights[i] < minweight && weights[i]>0) minweight = weights[i];
		ic1[i] = 0;
		ic2[i] = 1;
		for (il = 0; il < 2; ++il) {
			dt[il] = 0.;
			for (j = 0; j < n; ++j) {
				da = a[i][j] - c[il][j];
				dt[il] += da * da ;
			}
		}
		if (dt[0] >= dt[1]) {
			ic1[i] = 1;
			ic2[i] = 0;
			temp = dt[0];
			dt[0] = dt[1];
			dt[1] = temp;
		}
		for (l = 2; l < k; ++l) {
			db = 0.;
			for (j = 0; j < n; ++j) {
				dc = a[i][j] - c[l][j];
				db += dc * dc ;
				if (db > dt[1]) j = n; /*go to the end of the loop -- end of story*/
			}
			if (db < dt[1]) {
				if (db >= dt[0]) {
					dt[1] = db;
					ic2[i] = l;
				}
				else {
					dt[1] = dt[0];
					ic2[i] = ic1[i];
					dt[0] = db;
					ic1[i] = l;
				}
			}
		}
	}



	
/* Update cluster centres to be the average of prints contained within them. */
	
	MAKE_VECTOR(d, m);
	MAKE_VECTOR(itran, k);
	MAKE_VECTOR(live, k);
	MAKE_VECTOR(ncp, k);
	
	for (l = 0; l < k; ++l) {
		nc[l] = 0.;
		for (j = 0; j < n; ++j) c[l][j] = 0.;
	}
	/* only use points with weight > 0 - can use negative numbers as error/NA/etc */
	for (i = 0; i < m; ++i) {
	  if(weights[i] >=minweight){
		nc[ic1[i]] += weights[i];
		for (j = 0; j < n; ++j) c[ic1[i]][j] += weights[i]*a[i][j];
	  }
	}
	
	/*  Check to see if there is any empty cluster at this stage */
	
	for (l = 0; l < k; ++l) {
		if (nc[l] < minweight) {
			*ifault = 1;
			FREE_VECTOR(ic2);
			//FREE_VECTOR(an1);
			//FREE_VECTOR(an2);
			/* weighted kmeans doesn't use an1 and an2 */
			FREE_VECTOR(d);
			FREE_VECTOR(itran);
			FREE_VECTOR(live);
			FREE_VECTOR(ncp);
			return;
		}
		for (j = 0; j < n; ++j) c[l][j]/=nc[l];
				
		/* Initialize itran & ncp: 
		   itran[l] = 1 if cluster l is updated in the quick-transfer 
		              stage, 0 ow.
	           In the optimal-transfer stage, ncp(l) stores the step at 
		   which cluster l is last updated. 
		   In the quick-transfer stage, ncp(l) stores the step 
		   at which cluster l is last updated plus m. */
    
		itran[l] = 1;
		ncp[l] = -1;
	}
	
	indx = 0;
	for (ij = 0;ij < iter; ++ij) {
		/* In this stage, there is only one pass through the data. 
		   Each print is re-allocated, if necessary, to the cluster 
		   that will induce the maximum reduction in within-cluster 
		   sum of squares. */
	  woptra(a, m, n, weights, c, k, ic1, ic2, nc, ncp, d, itran,
		      live, &indx);
		/* Stop if no transfer took place in the last m optimal 
		   transfer steps. */
		if (indx == m) 	ij=iter;
		else { /* Each print is tested in turn to see if it should be 
			  re-allocated to the cluster to which it is most 
			  likely to be transferred, ic2[i], from its present 
			  cluster, ic1[i]. Loop through the data until no 
			  further change is to take place. */

		  wqtran(a, m, n, weights, c, k, ic1, ic2, nc, ncp, d,
			      itran, &indx);
      			if (k==2) ij = iter;/* k=2: no need to re-enter the optimal transfer stage*/
			else {	/* ncp has to be set to 0 before entering optra. */
			  for (l=0;l<k; ++l) ncp[l] = 0;
			}
		}
	} 
	if ((indx != m) && (k != 2)) *ifault = 2; /* iterations exceeded: may indicate unforeseen looping */
		
	for (l = 0; l < k; ++l) { /* Compute within-cluster SS for each cluster. */
		wss[l] = 0.;
		for (j = 0; j < n; ++j) c[l][j] = 0.;
	}
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j)
		  if(weights[i]>=minweight) c[ic1[i]][j] += a[i][j] * weights[i];
	}
	/* has this always been this inefficient?! */
	for (j = 0; j < n; ++j) {
		for (l = 0; l < k; ++l) c[l][j] /= (double)nc[l];
		for (i = 0; i < m; ++i) {
			da = a[i][j] - c[ic1[i]][j];
			if(weights[i]>=minweight) wss[ic1[i]] += weights[i] * da * da;
		}
	}
	FREE_VECTOR(ic2);
	//FREE_VECTOR(an1);
	FREE_VECTOR(d);
	FREE_VECTOR(itran);
	FREE_VECTOR(live);
	FREE_VECTOR(ncp);
	//FREE_VECTOR(an2);
	return;
} 

void woptra(double **a, unsigned int m, unsigned int n, double *weights, double **c, unsigned int k, unsigned int *ic1, unsigned int *ic2,
	   double *nc, unsigned int *ncp, double *d, unsigned int *itran,
	   unsigned int *live, unsigned int *indx)
{
	unsigned int i, j, l, l1, l2, ll, flag;
	double r2, da, db, dc, dd, de, df, rr, al1, al2, alt, alw;
	double an2l2;
	double an2l1;
	double an1l2;
	double an1l1;
	
	/* ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1 
	   This is the optimal transfer stage. Each print is re-allocated, if 
	   necessary, to the cluster that will induce a maximum reduction in 
	   the within-cluster sum of squares. 
	   If cluster L is updated in the last quick-transfer stage, it belongs
	   to the live set throughout this stage. Otherwise, at each step, it 
	   is not in the live set if it has not been updated in the last M 
	   optimal transfer steps. */


	

	for (l = 0; l < k; ++l) {
		if (itran[l] == 1) {
			live[l] = m + 1;
		}
	}
	for (i = 0; i < m; ++i) {

		(*indx)++;
		l1 = ic1[i];
		l2 = ic2[i];
		ll = l2;
		 an2l2 = nc[l2]/(nc[l2] + weights[i]);
		 an2l1 = nc[l2]/(nc[l1] + weights[i]);
		 an1l2 = BIG;
		 an1l1 = BIG;
	  if(weights[i] > 0.){		
		 if (nc[l2] > weights[i]){
		  an1l2 = nc[l2] / (nc[l2]-weights[i]);
		}
		if (nc[l1] > weights[i]){
		  an1l1 = nc[l1] / (nc[l1]-weights[i]);
		}

		/* If print I is the only member of cluster L1, no transfer. */
		if (nc[l1] > weights[i]){
			/* If L1 has not yet been updated, no need to re-compute D(I).*/
			if (ncp[l1] != 0) {
				de = 0.;
				for (j = 0; j < n; ++j) {
					df = a[i][j] - c[l1][j];
					de += df*df;
				}
				d[i] = de * an1l1 * weights[i];
			}
			/* Find the cluster with minimum R2. */
			da=0.;
			for (j = 0; j < n; ++j) {
				db = a[i][j] - c[l2][j];
				da += db * db;
			}
			r2 = da * an2l2 * weights[i];
			for (l = 0; l < k; ++l) {
				/* If I >= LIVE(L1), then L1 is not in the live
				   set.   If this is true, we only need to 
				   consider clusters that are in the live set 
				   for possible transfer of print I. Otherwise,
				   we need to consider all possible clusters. */
				if (((i >= ( live[l1] - 1)) && (i >= (live[l] - 1))) || (l == l1) || (l == ll)) {
				}
				else {
				  // rr = r2 / an2[l];
					rr = r2 / (nc[l] / (nc[l] + weights[i]));
					dc = 0.;
					flag = 0;
					j = 0;
					while ((!flag) && (j < n)) {
						dd = a[i][j] - c[l][j];
						dc += dd * dd * weights[i];
						if (dc >= rr) {
							flag = 1;
						}
						j++;
					}
					if (!flag) {
						r2 = dc * (nc[l] / (nc[l] + weights[i]));
						l2 = l;
					}
				}
			}
			if (r2 >= d[i]) {
                       /* If no transfer is necessary, L2 is the new IC2(I). */
				ic2[i] = l2;
			}
			else {/* Update cluster centres, LIVE, NCP 
				 for clusters L1 and L2, and update IC1(I) & 
				 IC2(I). */
				(*indx) = 0;
				live[l1] = m + i + 1;
				live[l2] = m + i + 1;
				ncp[l1] = i + 1;
				ncp[l2] = i + 1;
				al1 = (double)nc[l1];
				alw = al1 - weights[i];
				al2 = (double)nc[l2];
				alt = al2 + weights[i];
				for (j = 0; j < n; ++j) {
					c[l1][j] = (c[l1][j] * al1 - weights[i]*a[i][j]) / alw;
					c[l2][j] = (c[l2][j] * al2 + weights[i]*a[i][j]) / alt;
				}
				nc[l1] = nc[l1] - weights[i];
				nc[l2] = nc[l2] + weights[i];
				ic1[i] = l2;
				ic2[i] = l1;
			}
		}
		if ( (*indx) == m) return;
	  }
	}
	for (l = 0; l < k; ++l) { /* ITRAN(L) = 0 before entering QTRAN. Also, 
				     LIVE(L) has to be decreased by M before 
				     re-entering OPTRA. */ 
		itran[l] = 0;
		live[l] -= m;
	}
	return;
}

void wqtran(double **a,unsigned int m,unsigned int n,double *weights, double **c,unsigned int k,unsigned int *ic1,unsigned int *ic2,double *nc,
	  unsigned int *ncp,double *d,unsigned int *itran,unsigned int *indx)
{
  unsigned int i,j,l1,l2,icoun,istep,flag = 0,iflag;
	double r2,da,db,dd,de,al1,al2,alt,alw;
  
	/* ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1 
	   This is the quick transfer stage. IC1(I) is the cluster which print 
	   I belongs to. IC2(I) is the cluster which print I is most likely to 
	   be transferred to. For each print I, IC1(I) & IC2(I) are switched, 
	   if necessary, to reduce within-cluster sum of squares.  The cluster 
	   centres are updated after each step. In the optimal transfer stage, 
	   NCP(L) indicates the step at which cluster L is last updated. In the
	   quick transfer stage, NCP(L) is equal to the step at which cluster 
	   L is last updated plus M. */
	icoun = 0;
	istep = 0;
	iflag=0;
	double an2l2 ;
	double an2l1 ;
	double an1l2 ;
	double an1l1 ;
	while (iflag==0) {
		/*Remove this as iflag creates an endless loop since it never gets updated. No!*/
		for (i = 0; i < m; i++) {

			++icoun;
			++istep;
			
		    if(weights[i] > 0.){
			l1 = ic1[i];
			l2 = ic2[i];
			an2l2 = nc[l2]/(nc[l2] + weights[i]);
			an2l1 = nc[l2]/(nc[l1] + weights[i]);
			an1l2 = BIG;
			an1l1 = BIG;
			if (nc[l2] > weights[i]){
			  an1l2 = nc[l2] / (nc[l2]-weights[i]);
			}
			if (nc[l1] > weights[i]){
			  an1l1 = nc[l1] / (nc[l1]-weights[i]);
			}
			/*  If print I is the only member of cluster L1, no 
			    transfer. */
			if (nc[l1] > weights[i]) {
				/* If ISTEP > NCP(L1), no need to re-compute distance 
				   from print I to cluster L1.   Note that if cluster
				   L1 is last updated exactly M steps ago, we still
				   need to compute the distance from print I to 
				   cluster L1. */
				if (istep <= ncp[l1]) {
					da=0.;
					for (j = 0; j < n; ++j) {
						db = a[i][j] - c[l1][j];
						da += db * db * weights[i];
					}
					d[i] = da * an1l1;
					
					/* If ISTEP >= both NCP(L1) & NCP(L2) 
					   there will be no transfer of print 
					   I at this step. */
				}
				/*if ((istep<ncp[l1]) || (istep < ncp[l2])) {*/
				if (((istep < ncp[l1]) || (istep < ncp[l2])) && (!flag)) {
					r2 = d[i] / an2l2;
					dd = 0.;
					flag = 0;
					j = 0;
					while ((j < n) && (!flag)) {
						de = a[i][j] - c[l2][j];
						dd += de*de*weights[i];
						if (dd >= r2) flag=1;
						j++;
					}
					if (!flag) {
						/* Update cluster centres, NCP,
						   NC, ITRAN, AN1 & AN2 for 
						   clusters L1 & L2. Also 
						   update IC1(I) & IC2(I). Note
						   that if any updating occurs 
						   in this stage, INDX is set 
						   back to 0. */ 
						icoun = 0;
						(*indx) = 0;
						itran[l1] = 1;
						itran[l2] = 1;
						ncp[l1] = istep+m;
						ncp[l2] = istep+m;
						al1 = (double)nc[l1];
						alw = al1-weights[i];
						al2 = (double)nc[l2];
						alt = al2+weights[i];
						for (j = 0;j < n; ++j) {
							c[l1][j] = (c[l1][j] * al1 - weights[i]*a[i][j]) / alw;
							c[l2][j] = (c[l2][j] * al2 + weights[i]*a[i][j]) / alt;
						}
						nc[l1] = nc[l1] - weights[i];
						nc[l2] = nc[l2] + weights[i];
						ic1[i] = l2;
						ic2[i] = l1;
					}
				}
			}
			if (icoun==m) return;
		  }
		}
	}
}
