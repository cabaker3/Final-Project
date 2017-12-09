#include <iostream>
#include <omp.h>
#include <fstream>
#include "stopwatch.hpp"

int main(int argc, char *argv[]) {
	
	int N = atoi(argv[1]);
	
	omp_set_num_threads(N);

	double start = omp_get_wtime();
	
	double stop = omp_get_wtime();
	
	double time = (stop - start) * 1000;
  
}
