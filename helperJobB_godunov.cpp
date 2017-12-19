#include "omp.h"
#include <iostream>
#include <math.h>

//run on the GPU

float *helperJobB_godunov(float alpha[101][1], float E[3][101], float Qold[3][101], float Qnew[3][101], float F[3][100], float IM){
  
  for(int j = 1; j <= 3; j++){
    for(int i = 1; i <= IM; i++)
      F[j][i] = 0.5 * (E[j][i] + E[j][i+1]) - 0.5 * abs(alpha[i][]) * (Qold[j][i+1] - Qnew[j][i]);
  }
  
  return F;
}
