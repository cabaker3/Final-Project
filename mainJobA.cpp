#include <iostream>
#include "omp.h"
#include <fstream>
#include "helperJobB.cpp"
#include <math.h>

//run on the host using OpenMP
//time integration of a deforming elastic body
  //done in FEA using numerical integration
  //implicit integration needs a Jacobian J
float mainJobA(n_t,n_x,a,L,tmax,error_plots){
  float dx = L/(n_x - 1);
  float dt = tmax/(n_t - 1);
  float r = alpha * dt / pow(dx,2);
  float r2 = 1 - 2*r;
  
  //create arrays to save data for export
  //x = ... L,n_x
  //t = ... tmax,n_t
  //U = ... n_x,n_t
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
