#include<stdio.h>
#include<stdlib.h>
#include <cuda.h>
#include <iostream>
#include <math.h>

__global__ void helperJobB(float *alpha, float *E, float *Qold, float *Qnew, float *F, float IM){
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  
  for(int j = 1; j <= 3; j++){
    for(int i = 1; i <= IM; i++)
      F[j][i] = 0.5 * (E[j][i] + E[j][i+1]) - 0.5 * abs(alpha[i][1]) * (Qold[j][i+1] - Qnew[j][i]);
  }
  
}
