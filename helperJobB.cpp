#include "omp.h"
#include <iostream>
#include <fstream>

//run on the GPU
//goal: update entries in J
  //might be fast enough to update all entries in J before start of a new time step
  //update might be slow to the point where only some percent of J entries are update
    //then should keep updating the remaining stale entries of J over the next time step