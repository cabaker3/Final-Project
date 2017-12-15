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
  
  //Lax Initial Conditions
  for(int x = 0; x <= 1; x+=dx){
    if(x <= 0.5){
      //Qi[i] = [0.445; 0.311; 8.928];
    }else{
      //Qi[i] = [0.5; 0; 1.4275];  
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
  
  for(int t = 0; t <= 0.16; t+=dt){
    //call helperJobB
    helperJobB(); //Flux & Qnew
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
