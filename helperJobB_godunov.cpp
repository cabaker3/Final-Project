#include "omp.h"
#include <iostream>
#include <fstream>
#include <math.h>

//run on the GPU
//goal: update entries in J
  //might be fast enough to update all entries in J before start of a new time step
  //update might be slow to the point where only some percent of J entries are update
    //then should keep updating the remaining stale entries of J over the next time step

float helperJobB(float alpha, float E[][], float Qold[][], float Qnew[][], float F[][], float Qn1[][]){
  //create 2d array F and Qn1???
  
  for(int j = 1; j <= 3; j++){
    for(int i = 1; i <= IM; i++)
      F[j][i] = 0.5 * (E[j][i] + E[j][i+1]) - 0.5 * abs(alpha[i][]) * (Qold[j][i+1] - Qnew[j][i]);
  }
  
  for(int j = 1; j <= 3; j++){
    for(int i = 2; i <= IM; i++)
      Qn1[j][i] = Qold[j][i] - (dt/dx)*(F[j][i-1]);
  }
  
  return F, Qn1;
}
