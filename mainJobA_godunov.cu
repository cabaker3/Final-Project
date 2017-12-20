#include <iostream>
#include "omp.h"
#include "helperJobB_kernel.cu"
#include <math.h>
#include <utility>
#include <cuda.h>

using namespace std;

//run on the host using OpenMP

void mainJobA_godunov(int L, float g, float dx, float dt, float IM){
  //constant
  int i = 1;
  const int row1 = 101;
  const int row2 = 3;
  const int col1 = 1;
  const int col2 = 101;
  const int col3 = 81;
  dim3 dimBlock(3*100, 3*100);
  dim3 dimGrid(1, 1);
  
  //create 2D array Qi
  float **Qi = new float*[row2];
  for(int i = 0; i < row2; ++i) {
    Qi[i] = new float[col2];
  }
  
  //Lax Initial Conditions
  for(int x = 0; x <= 1; x+=dx){
    if(x <= 0.5){
      Qi[1][i] = 0.445;
      Qi[2][i] = 0.311;
      Qi[3][i] = 8.928;
    }else{
      Qi[1][i] = 0.5;  
      Qi[2][i] = 0;
      Qi[3][i] = 1.4275;
    }
    i += 1; //change to for loop
  }
  
  //create array Qold = Qi
  float **Qold = new float*[row2];
  float **Qnew = new float*[row2];
  for(int i = 0; i < row2; ++i) {
    Qold[i] = new float[col2];
    Qnew[i] = new float[col2];
  }
  memcpy(Qold, Qi, sizeof(Qold));
  //create array Qnew = Qold
  memcpy(Qnew, Qold, sizeof(Qnew));
  
  //Initial Flow Properties
  
  
  //create array rhoi, ui, eti, pi, a, E, eigen
  float **rhoi = new float*[row1];
  float **ui = new float*[row1];
  float **eti = new float*[row1];
  float **pi = new float*[row1];
  float **a = new float*[row1];
  float **alpha = new float*[row1];
  for(int i = 0; i < row1; ++i) {
    rhoi[i] = new float[col1];
    ui[i] = new float[col1];
    eti[i] = new float[col1];
    pi[i] = new float[col1];
    a[i] = new float[col1];
    alpha[i] = new float[col1];
  }
  
  float **E = new float*[row2];
  float **eigen = new float*[row2];
  for(int i = 0; i < row2; ++i) {
    E[i] = new float[col2];
    eigen[i] = new float[col2];
  }
  
  #pragma omp for
  for(int i = 1; i < IM+1; i++){
    //Density
    rhoi[i][1] = Qi[1][i];
    
    //Velocity
    ui[i][1] = Qi[2][i] / rhoi[i][1];
    
    //Total Energy
    eti[i][1] = Qi[3][i] / rhoi[i][1];
    
    //Pressure, from the equation of state
    pi[i][1] = (g-1) * (rhoi[i][1]) * eti[i][1] - 0.5 * rhoi[i][1] * pow(ui[i][1],2));
      
    //Speed of Sound
    a[i][1] = sqrt(g*pi[i][1]/rhoi[i][1]);
    
    //Intial E Matrix
    E[1][i] = rhoi[i][1]*ui[i][1];
    E[2][i] = rhoi[i][1]*pow(ui[i][1],2) + pi[i][1];
    E[3][i] = eti[i][1]*rhoi[i][1]*ui[i][1]+pi[i][1]*ui[i][1]};
    
    //Eigenvalues
    eigen[1][i] = ui[i][1];
    eigen[2][i] = ui[i][1] + a[i][1]; 
    eigen[3][i] = ui[i][1]-a[i][1];
  }
  
  #pragma omp for
  for(int i = 0; i < row1; ++i) {
    alpha[i][1] = max(abs(eigen));
  }
  
  delete[] rhoi;
  delete[] ui;
  delete[] eti;
  delete[] pi;
  delete[] a;
  delete[] E;
  delete[] eigen;
  
  int k = 1;
  
  //initialize rho,u,et,p,a,E,eigen, F,Qn1
  float **rho = new float*[row1];
  float **u = new float*[row1];
  for(int i = 0; i < row1; ++i) {
    rho[i] = new float[col1];
    u[i] = new float[col1];
    et[i] = new float[col1];
    p[i] = new float[col1];
    a[i] = new float[col1];
  }
 
  float **E = new float*[row2];
  float **eigen = new float*[row2];
  float **Qn1 = new float*[row2];
  for(int i = 0; i < row2; ++i) {
    E[i] = new float[col2];
    eigen[i] = new float[col2];
    Qn1[i] = new float[col2];
  }

  float **rhom = new float*[row1];
  float **um = new float*[row1];
  float **etm = new float*[row1];
  float **pm = new float*[row1];
  for(int i = 0; i < row1; ++i) {
    rhom[i] = new float[col3];
    um[i] = new float[col3];
    etm[i] = new float[col3];
    pm[i] = new float[col3];
  }

  float **F = new float*[row2];
  for(int i = 0; i < row2; ++i) {
    F[i] = new float[100];
  }
  
  #pragma omp for
  for(int t = 0; t <= 0.16; t+=dt){
    //call helperJobB
    helperJobB_godunov<<<dimGrid, dimBlock>>>(alpha,E,Qold,Qnew,F,IM); //Flux
    
    #pragma omp barrier
    for(int j = 1; j <= 3; j++){
      for(int i = 2; i <= IM; i++){
        Qn1[j][i] = Qold[j][i] - (dt/dx) * (F[j][i] - F[j][i-1]); 
      }
    }
    
    for(int x = 0; x <= 3; x++){
    Qn1[x][1] = Qi[x][1];
    Qn1[x][101] = Qi[x][101];
    }
    
    memcpy(Qnew, Qn1, sizeof(Qnew));
    memcpy(Qold, Qnew, sizeof(Qold));
    
    #pragma omp for
    for(int i = 1; i < IM+1; i++){
    //Density
    rho[i][1] = Qnew[1][i];
    
    //Velocity
    u[i] = Qnew[2][i] / rho[i][1];
    
    //Total Energy
    et[i][1] = Qnew[3][i] / rho[i][1];
    
    //Pressure, from the equation of state
    p[i][1] = (g-1) * (rho[i][1] * et[i][1] - 0.5 * rho[i][1] * pow(u[i][1],2));
      
    //Speed of Sound
    a[i][1] = sqrt(g*p[i][1]/rho[i][1]);
    
    //Intial E Matrix
    E[1][i] = rho[i][1]*u[i][1];
    E[2][i] = rho[i][1]*pow(ui[i][1],2) + p[i][1];
    E[3][i] = et[i][1]*rho[i][1]*u[i][1]+p[i][1]*u[i][1];
    
    //Eigenvalues
    eigen[1][i] = u[i][1];
    eigen[2][i] = u[i][1] + a[i][1]; 
    eigen[3][i] = u[i][1]-a[i][1];
  }
    //Alpha
    alpha = max(abs(eigen));
    
    for(int x = 0; x <= 101; x++){
      um[x][k] =u[x][1];
      rhom[x][k] = rho[x][1];
      pm[x][k] = p[x][1];
      etm[x][k] = et[x][1];
    }
    
    k += 1;
  }
  
  delete[] rho;
  delete[] u;
  delete[] et;
  delete[] p;
  delete[] a;
  delete[] E;
  delete[] eigen;
  delete[] alpha;
  delete[] F;
  delete[] Qn1;
  delete[] Qi;
  delete[] Qnew;
  delete[] Qold;
  
  // display properties matrices 
  cout << rhom << "\n" << um << "\n" << etm << "\n" << pm << "\n";
  
  delete[] rhom;
  delete[] um;
  delete[] etm;
  delete[] pm;
}
