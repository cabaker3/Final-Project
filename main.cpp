#include <iostream>
#include <omp.h>
#include <fstream>
#define OMPI_SKIP_MPICXX  /* Don't use OpenMPI's C++ bindings (they are deprecated) */
#include <mpi.h>
#include "mainJobA.cpp"
#include "helperJobB.cpp"

int main(int argc, char *argv[]) {
	//software backbone/infrastructure that allows mainJobA to carry on while helperJobB works asynchronously to update
	//the computationally costly ingredient needed by mainJobA
	
	//mainJobA runs on one node, helperJobB runs on another node; updates via MPI
	
	int N = atoi(argv[1]);
	
	omp_set_num_threads(N);

	double start = omp_get_wtime();
	
	double stop = omp_get_wtime();
	
	double time = (stop - start) * 1000;
  
}
