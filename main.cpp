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
	int L;
	float g, dx, dt, IM;
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " N\n";
		return -1;
	}
	
	const auto num_threads = atoi(argv[1]);
	omp_set_num_threads(num_threads);
	
	if (argc != 3) {
		L = 1;
	}else{
		L = atoi(argv[2]);
	}
	
	if (argc != 4) {
		g = 1.4;
	}else{
		g = atoi(argv[3]);
	}
	
	if (argc != 5) {
		dx = 0.01;
	}else{
		a = atoi(argv[4]);
	}
	
	if (argc != 6) {
		dt = 1/2000;
	}else{
		dt = atoi(argv[5]);
	}
	
	if (argc != 7) {
		IM = L/dx;
	}else{
		L = atoi(argv[6]);
	}
	
	errPlots = 1;

	double start = omp_get_wtime();
	
	//create array to store result
	
	#pragma omp parallel shared (nt,n_x,a,L,tmax,errPlots)
		//mainJobA(nt,n_x,a,L,tmax,errPlots); //save to array
	
	double stop = omp_get_wtime();
	
	double time = (stop - start) * 1000;
	
	std::cout << time << endl;
	//print out result
  
}
