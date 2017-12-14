#include "omp.h"
#include <iostream>
#include <fstream>
#include <math.h>

//run on the GPU
//goal: update entries in J
  //might be fast enough to update all entries in J before start of a new time step
  //update might be slow to the point where only some percent of J entries are update
    //then should keep updating the remaining stale entries of J over the next time step

float helperJobB(float U[nt][n_x],int nt, int n_x, float r, float r2){
  for(int m = 2; m <= nt; m++){
    for(int i = 2; i <= (n_x - 1); i++){
      U[i][m] = r*U[i-1][m-1] + r2*U[i][m-1] + r*U[i+1][m-1]; 
    }
  }
  
  return U;
}
