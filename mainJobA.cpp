#include <iostream>
#include "omp.h"
#include <fstream>
#include "helperJobB.cpp"
#include <math.h>

//run on the host using OpenMP
//time integration of a deforming elastic body
  //done in FEA using numerical integration
  //implicit integration needs a Jacobian J
float mainJobA(int nt,int n_x,float a, int L, float tmax, int error_plots){
  float dx = L/(n_x - 1);
  float dt = tmax/(nt - 1);
  float r = alpha * dt / pow(dx,2);
  float r2 = 1 - 2*r;
  
  //set up iterations
  
  //create arrays to save data for export
  //x = ... L,n_x
  float **x = new float*[L];
  for(int i = 0; i < L; ++i) {
    x[i] = new float[n_x];
  }
  //t = ... tmax,n_t
  float **t = new float*[tmax];
  for(int i = 0; i < tmax; ++i) {
    t[i] = new float[nt];
  }
  //U = ... n_x,n_t
  float **U = new float*[n_x];
  for(int i = 0; i < n_x; ++i) {
    U[i] = new float[nt];
  }
  
  //call helperJobB
  helperJobB();
}
/*int max_itr = 10;

for (int itr =0; itr < max_itr ; itr ++)
{
  // compute the update for the assigned x[i]
  // call helperJobB
#pragma omp barrier
{
  // write x[i] in global memory
}
#pragma omp barrier
}*/
