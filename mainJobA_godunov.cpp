#include <iostream>
#include "omp.h"
#include "helperJobB_godunov.cpp"
#include <math.h>
#include <utility>
#include <type_traits>
#include <typeinfo>
#include <cxxabi.h>

using namespace std;

//run on the host using OpenMP

void mainJobA_godunov(int L, float g, float dx, float dt, float IM){
  //constant
  int i = 1;
  int row1 = 101;
  int row2 = 3;
  int col1 = 1;
  int col2 = 101;
  int col3 = 81;
  
  //create 2D array Qi
  auto Qi = new float [row2][col2];
  
  //Lax Initial Conditions
  for(int x = 0; x <= 1; x+=dx){
    if(x <= 0.5){
      Qi[][i] = {0.445, 0.311, 8.928};
    }else{
      Qi[][i] = {0.5, 0, 1.4275};  
    }
    i += 1; //change to for loop
  }
  
  //create array Qold = Qi
  auto Qold = new float [row2][col2];
  memcpy(Qold, Qi, sizeof(Qold));
  //create array Qnew = Qold
  auto Qnew = new float [row2][col2];
  memcpy(Qnew, Qold, sizeof(Qnew));
  
  //Initial Flow Properties
  
  
  //create array rhoi, ui, eti, pi, a, E, eigen
  auto rhoi = new float [row1][col1];
  auto ui = new float [row1][col1];
  auto eti = new float [row1][col1];
  auto pi = new float [row1][col1];
  auto a = new float [row1][col1];
  auto E = new float [row2][col2];
  auto eigen = new float [row2][col2];
  
  auto alpha = new float [row1][col1];
  
  #pragma omp for
  for(int i = 1; i < IM+1; i++){
    //Density
    rhoi[i][] = Qi[1][i];
    
    //Velocity
    ui[i] = Qi[2][i] / rhoi[i][];
    
    //Total Energy
    eti[i][] = Qi[3][i] / rhoi[i][];
    
    //Pressure, from the equation of state
    pi[i][] = (g-1) * (rhoi[i][]) * eti[i][] - 0.5 * rhoi[i][] * pow(ui[i][],2));
      
    //Speed of Sound
    a[i][] = sqrt(g*pi[i][]/rhoi[i][]);
    
    //Intial E Matrix
    E[][i] = {rhoi[i][]*ui[i][], rhoi[i][]*pow(ui[i][],2) + pi[i][], eti[i][]*rhoi[i][]*ui[i][]+pi[i][]*ui[i][]};
    
    //Eigenvalues
    eigen[][i] = {ui[i][], ui[i][] + a[i][], ui[i][]-a[i][]};
  }
  
  float alpha = max(abs(eigen));
  
  delete[] rhoi;
  delete[] ui;
  delete[] eti;
  delete[] pi;
  delete[] a;
  delete[] E;
  delete[] eigen;
  
  int k = 1;
  
  //initialize rho,u,et,p,a,E,eigen, F,Qn1
  auto rho = new float [row1][col1];
  auto u = new float [row1][col1];
  auto et = new float [row1][col1];
  auto p = new float [row1][col1];
  auto a = new float [row1][col1];
  auto E = new float [row2][col2];
  auto eigen = new float [row2][col2];
  auto F = new float [row2][100];
  auto Qn1 = new float [row2][col2];
  auto rhom = new float [row1][col3];
  auto um = new float [row1][col3];
  auto etm = new float [row1][col3];
  auto pm = new float [row1][col3];
  
  #pragma omp for
  for(int t = 0; t <= 0.16; t+=dt){
    //call helperJobB
    memcpy(F,helperJobB_godunov(alpha,E,Qold,Qnew,F,IM),sizeof(F)); //Flux
    
    #pragma omp for
    for(int j = 1; j <= 3; j++){
      for(int i = 2; i <= IM; i++){
        Qn1[j][i] = Qold[j][i] - (dt/dx) * (F[j][i] - F[j][i-1]); 
      }
    }
    
    Qn1[][1] = Qi[][1];
    Qn1[][101] = Qi[][101];
    Qnew = Qn1;
    Qold = Qnew;
    
    #pragma omp for
    for(int i = 1; i < IM+1; i++){
    //Density
    rho[i][] = Qnew[1][i];
    
    //Velocity
    u[i] = Qnew[2][i] / rho[i][];
    
    //Total Energy
    et[i][] = Qnew[3][i] / rho[i][];
    
    //Pressure, from the equation of state
    p[i][] = (g-1) * (rho[i][] * et[i][] - 0.5 * rho[i][] * pow(u[i][],2));
      
    //Speed of Sound
    a[i][] = sqrt(g*p[i][]/rho[i][]);
    
    //Intial E Matrix
    E[][i] = {rho[i][]*u[i][], rho[i][]*pow(ui[i][],2) + p[i][], et[i][]*rho[i][]*u[i][]+p[i][]*u[i][]};
    
    //Eigenvalues
    eigen[][i] = {u[i][], u[i][] + a[i][], u[i][]-a[i][]};
  }
    //Alpha
    alpha = max(abs(eigen));
    
    um[][k] =u[][1];
    rhom[][k] = rho[][1];
    pm[][k] = p[][1];
    etm[][k] = et[][1];
    
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
  
  //cout
  cout << rhom << "\n" << um << "\n" << etm << "\n" << pm << "\n";
  
  delete[] rhom;
  delete[] um;
  delete[] etm;
  delete[] pm;
  //return ;
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
