# vim:set ts=8 sw=8 sts=0 noet:

#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#    File:         Makefile                                                  */
#    Description:  Makefile for programs running a simple k-means clustering */
#                  algorithm                                                 */
#                                                                            */
#    Author:  Wei-keng Liao                                                  */
#             ECE Department Northwestern University                         */
#             email: wkliao@ece.northwestern.edu                             */
#    Copyright, 2005, Wei-keng Liao                                          */
#                                                                            */
#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

# Copyright (c) 2005 Wei-keng Liao
# Copyright (c) 2017 Geoffrey Z Thompson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# ------------------------------------------------------------------------------

.KEEP_STATE:

all: seq mpi

DFLAGS      =
OPTFLAGS    = -O -NDEBUG
OPTFLAGS    = -g -pg
INCFLAGS    = -I. -I/usr/include/R/
CFLAGS      = $(OPTFLAGS) $(DFLAGS) $(INCFLAGS) -DBLOCK_SHARED_MEM_OPTIMIZATION=1
LDFLAGS     = $(OPTFLAGS)
LIBS        = -lm -lRmath

CC          = gcc
MPICC       = mpicc

.c.o:
	$(CC) $(CFLAGS) -c $<

H_FILES     = kmeans-mpi.h array.h util.h

#------   MPI version -----------------------------------------
MPI_SRC     = mpi_main.c   \
              mpi_kmeans.c \
              mpi_io.c     \
	      file_io.c    \
	      wkmeans.c    \
	      mpi_scalablekmeanspp.c

MPI_OBJ     = $(MPI_SRC:%.c=%.o)


wkmeans.o: wkmeans.c $(H_FILES)
	$(MPICC) $(CFLAGS) -c $*.c  

mpi_scalablekmeanspp.o: mpi_scalablekmeanspp.c $(H_FILES)
	$(MPICC) $(CFLAGS) -c $*.c  

mpi_main.o: mpi_main.c $(H_FILES)
	$(MPICC) $(CFLAGS) -c $*.c -lm -lRmath

mpi_kmeans.o: mpi_kmeans.c $(H_FILES)
	$(MPICC) $(CFLAGS) -c $*.c 

mpi_io.o: mpi_io.c $(H_FILES)
	$(MPICC) $(CFLAGS) -c $*.c

mpi: mpi_main
mpi_main: $(MPI_OBJ) $(H_FILES)
	$(MPICC) $(LDFLAGS) -o mpi_main $(MPI_OBJ) $(LIBS)

#------   sequential version -----------------------------------------
SEQ_SRC     = seq_main.c   \
              seq_kmeans.c \
	      file_io.c    \
	      wtime.c

SEQ_OBJ     = $(SEQ_SRC:%.c=%.o)

$(SEQ_OBJ): $(H_FILES)

seq: seq_main
seq_main: $(SEQ_OBJ) $(H_FILES)
	$(CC) $(LDFLAGS) -o seq_main $(SEQ_OBJ) $(LIBS)

#---------------------------------------------------------------------
clean:
	rm -rf *.o omp_main seq_main mpi_main cuda_main \
	       core* .make.state gmon.out     \
               *.cluster_centres *.membership \
               Image_data/*.cluster_centres   \
               Image_data/*.membership        \
               profiles/
