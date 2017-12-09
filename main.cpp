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
	
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " N\n";
		return -1;
	}
	
	const auto num_threads = atoi(argv[1]);
	omp_set_num_threads(num_threads);
	
	//read in file for mainJobA
	//read in file for helperJobB

	double start = omp_get_wtime();
	
	double stop = omp_get_wtime();
	
	double time = (stop - start) * 1000;
	
	std::cout << time << endl;
  
}
