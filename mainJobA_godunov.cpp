#include <iostream>
#include "omp.h"
#include "helperJobB.cpp"
#include <math.h>

//run on the host using OpenMP
//time integration of a deforming elastic body
  //done in FEA using numerical integration
  //implicit integration needs a Jacobian J
float mainJobA(int L, float g, dx, dt, IM){
  //constant
  int i = 1;
  
  //create 2D array Qi
  auto Qi = new float [][];
  
  //Lax Initial Conditions
  for(int x = 0; x <= 1; x+=dx){
    if(x <= 0.5){
      Qi[1][i] = {0.445; 0.311; 8.928};
    }else{
      Qi[1][i] = {0.5; 0; 1.4275};  
    }
    i += 1; //change to for loop
  }
  
  //create array Qold = Qi
  //create array Qnew = Qold
  
  //Initial Flow Properties
  //create array rhoi, ui, eti, pi, a, E, eigen
  for(int i = 1; i < IM+1; i++){
    //Density
    rhoi[i][] = Qi[1][i];
    
    //Velocity
    ui[i] = Qi[2][i] / rhoi[i][];
    
    //Total Energy
    eti[i][] = Qi[3][i] / rhoi[i][];
    
    //Pressure, from the equation of state
    pi[i][] = ();
      
    //Speed of Sound
    a[i][] = sqrt(g*pi[i][]/rhoi[i][]);
    
    //Intial E Matrix
    E[][i] = [rhoi[i][]*ui[i][], rhoi[i][]*pow(ui[i][],2) + pi[i][], eti[i][]*rhoi[i][]*ui[i][]+pi[i][]*ui[i][]];
    
    //Eigenvalues
    eigen[][i] = [ui[i][], ui[i][] + a[i][], ui[i][]-a[i][]];
  }
  
  float alpha = max(abs(eigen));
  
  int k = 1;
  
  //initialize rho,u,et,p
  
  for(int t = 0; t <= 0.16; t+=dt){
    //call helperJobB
    helperJobB(); //Flux & Qnew
    
    Qn1[][1] = Qi[][1];
    Qn1[][101] = Qi[][101];
    Qnew = Qn1;
    Qold = Qnew;
    
    for(int i = 1; i < IM+1; i++){
    //Density
    rho[i][] = Qnew[1][i];
    
    //Velocity
    u[i] = Qnew[2][i] / rho[i][];
    
    //Total Energy
    et[i][] = Qnew[3][i] / rho[i][];
    
    //Pressure, from the equation of state
    p[i][] = ();
      
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
