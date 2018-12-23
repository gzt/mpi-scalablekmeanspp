/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         mpi_main.c   (an MPI version)                             */
/*   Description:  This program shows an example on how to call a subroutine */
/*                 that implements a simple k-means clustering algorithm     */
/*                 based on Euclid distance.                                 */
/*   Input file format:                                                      */
/*                 ascii  file: each line contains 1 data object             */
/*                 binary file: first 4-byte integer is the number of data   */
/*                 objects and 2nd integer is the no. of features (or        */
/*                 coordinates) of each object                               */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strtok() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     /* getopt() */

#include <mpi.h>
int      _debug;
#include "kmeans-mpi.h"
#include "util.h"

int     mpi_kmeans(double**, int, int, int, double, int*, double**, MPI_Comm);
double** mpi_read(int, char*, int*, int*, MPI_Comm);
int     mpi_write(int, char*, int, int, int, double**, int*, int, MPI_Comm);
void mpi_scalablekmeanspp(double **, unsigned int , unsigned int , unsigned int , unsigned int , double , double **, int , MPI_Comm );


/*---< usage() >------------------------------------------------------------*/
static void usage(char *argv0, double threshold) {
    char *help =
        "Usage: %s [switches] -i filename -n num_clusters\n"
        "       -i filename    : file containing data to be clustered\n"
        "       -b             : input file is in binary format (default no)\n"
        "       -r             : output file in binary format (default no)\n"
        "       -n num_clusters: number of clusters (K must > 1)\n"
        "       -t threshold   : threshold value (default %.4f)\n"
        "       -o             : output timing results (default no)\n"
        "       -d             : enable debug mode\n";
    fprintf(stderr, help, argv0, threshold);
}

/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
           int     opt;
    extern char   *optarg;
    extern int     optind;
           int     i, j;
           int     isInFileBinary, isOutFileBinary;
           int     is_output_timing, is_print_usage;

           int     numClusters, numCoords, numObjs, totalNumObjs;
           int    *membership;    /* [numObjs] */
           char   *filename;
           double **objects;       /* [numObjs][numCoords] data objects */
           double **clusters;      /* [numClusters][numCoords] cluster center */
           double   threshold;
	   double  lmult;
	   double  Lfactor;
           double  timing, io_timing, clustering_timing, initialization_timing;
	   unsigned int rounds;
           int        rank, nproc, mpi_namelen;
           char       mpi_name[MPI_MAX_PROCESSOR_NAME];
           MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Get_processor_name(mpi_name,&mpi_namelen);

    /* some default values */
    _debug           = 0;
    threshold        = 0.001;
    numClusters      = 0;
    isInFileBinary   = 0;
    isOutFileBinary  = 0;
    is_output_timing = 0;
    is_print_usage   = 0;
    filename         = NULL;
    rounds = 5;
    lmult = 0.;
    Lfactor = 0.;

    while ( (opt=getopt(argc,argv,"p:i:n:t:abdl:L:orR:h"))!= EOF) {
        switch (opt) {
            case 'i': filename=optarg;
                      break;
            case 'b': isInFileBinary = 1;
                      break;
            case 'r': isOutFileBinary = 1;
                      break;
            case 't': threshold=atof(optarg);
                      break;
            case 'n': numClusters = atoi(optarg);
                      break;
            case 'o': is_output_timing = 1;
                      break;
            case 'd': _debug = 1;
                      break;
            case 'h': is_print_usage = 1;
                      break;
	    case 'l': lmult=atof(optarg);
	              break;
	    case 'L': Lfactor=atof(optarg);
	              break;		      
	    case 'R': rounds = (unsigned int) atoi(optarg);
	              break;
            default: is_print_usage = 1;
                      break;
        }
    }

    if (filename == 0 || numClusters <= 1 || is_print_usage == 1) {
        if (rank == 0) usage(argv[0], threshold);
        MPI_Finalize();
        exit(1);
    }

    if (_debug) printf("Proc %d of %d running on %s\n", rank, nproc, mpi_name);

    MPI_Barrier(MPI_COMM_WORLD);
    io_timing = MPI_Wtime();

    /* read data points from file ------------------------------------------*/
    objects = mpi_read(isInFileBinary, filename, &numObjs, &numCoords,
                       MPI_COMM_WORLD);

    if (_debug) { /* print the first 4 objects' coordinates */
        int num = (numObjs < 4) ? numObjs : 4;
        for (i=0; i<num; i++) {
            char strline[1024], strdouble[16];
            sprintf(strline,"%d: objects[%d]= ",rank,i);
            for (j=0; j<numCoords; j++) {
                sprintf(strdouble,"%10f",objects[i][j]);
                strcat(strline, strdouble);
            }
            strcat(strline, "\n");
            printf("%s",strline);
        }
    }

    timing            = MPI_Wtime();
    io_timing         = timing - io_timing;
    initialization_timing = timing;

    /* allocate a 2D space for clusters[] (coordinates of cluster centers)
       this array should be the same across all processes                  */
    clusters    = (double**) malloc(numClusters *             sizeof(double*));
    assert(clusters != NULL);
    clusters[0] = (double*)  malloc(numClusters * numCoords * sizeof(double));
    assert(clusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        clusters[i] = clusters[i-1] + numCoords;

    MPI_Allreduce(&numObjs, &totalNumObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* pick first numClusters elements in feature[] as initial cluster centers*/
    /* if (rank == 0) { */
    /*     for (i=0; i<numClusters; i++) */
    /*         for (j=0; j<numCoords; j++) */
    /*             clusters[i][j] = objects[i][j]; */
    /* } */

    /* instead use mpi_scalablekmeanspp */
    set_seed(1111+rank,12345+rank);
    /* can specify either exact number for L (lmult) or constant to multiply by K (Lfactor) */
    /* by default L is listed as 0, if unspecified, default to 2*k */
    /* if both are specified, use biggest, if neither, use 2.  */
    if(lmult < 1. && Lfactor < .001)    lmult = 2.0*numClusters;
    if(Lfactor*numClusters > lmult) lmult = Lfactor * numClusters;
    if(lmult * rounds < 1.1*numClusters) lmult = 1.5 * numClusters / (1.0*rounds);
    if(_debug && rank == 0) printf("L = %f \n", lmult);
    int root = 0;
    mpi_scalablekmeanspp(objects, numObjs, numCoords, numClusters, rounds, lmult , clusters, root, MPI_COMM_WORLD);


    
    MPI_Bcast(clusters[0], numClusters*numCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    timing            = MPI_Wtime();
    initialization_timing = timing - initialization_timing;
    clustering_timing = timing;

    /* membership: the cluster id for each data object */
    membership = (int*) malloc(numObjs * sizeof(int));
    assert(membership != NULL);

    /* start the core computation -------------------------------------------*/
    mpi_kmeans(objects, numCoords, numObjs, numClusters, threshold, membership,
               clusters, MPI_COMM_WORLD);

    free(objects[0]);
    free(objects);

    timing            = MPI_Wtime();
    clustering_timing = timing - clustering_timing;

    /* output: the coordinates of the cluster centres ----------------------*/
    mpi_write(isOutFileBinary, filename, numClusters, numObjs, numCoords,
              clusters, membership, totalNumObjs, MPI_COMM_WORLD);

    free(membership);
    free(clusters[0]);
    free(clusters);

    /*---- output performance numbers ---------------------------------------*/
    if (is_output_timing) {
      double max_io_timing, max_clustering_timing, max_initialization_timing;

        io_timing += MPI_Wtime() - timing;

        /* get the max timing measured among all processes */
        MPI_Reduce(&io_timing, &max_io_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&initialization_timing, &max_initialization_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&clustering_timing, &max_clustering_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("\nPerforming **** Simple Kmeans  (MPI) ****\n");
            printf("Num of processes = %d\n", nproc);
            printf("Input file:        %s\n", filename);
            printf("numObjs          = %d\n", totalNumObjs);
            printf("numCoords        = %d\n", numCoords);
            printf("numClusters      = %d\n", numClusters);
            printf("threshold        = %.4f\n", threshold);
            printf("I/O time           = %10.4f sec\n", max_io_timing);
	    printf("Initialize time           = %10.4f sec\n", max_initialization_timing);
            printf("Computation timing = %10.4f sec\n", max_clustering_timing);
        }
    }

    MPI_Finalize();
    return(0);
}

