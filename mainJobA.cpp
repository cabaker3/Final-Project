#include <iostream>
#include "omp.h"
#include <fstream>

//run on the host using OpenMP
//time integration of a deforming elastic body
  //done in FEA using numerical integration
  //implicit integration needs a Jacobian J

int max_itr = 10;

for (int itr =0; itr < max_itr ; itr ++)
{
  // compute the update for the assigned x[i]
#pragma omp barrier
{
  // write x[i] in global memory
}
#pragma omp barrier
}
