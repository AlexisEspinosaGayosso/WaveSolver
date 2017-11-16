#include "allvars.h"
#include "proto.h"

int get_neighbours (int nx, int ny)
{
  tn = -1;
  bn = -1;
  ln = -1;
  rn = -1;
  
  if (nx == 1 && ny == 1)
    decomposition = NO_DEC;
  else
    if (nx > 1 && ny == 1)
      decomposition = X_DEC;
    else
      if (ny > 1 && nx == 1)
        decomposition = Y_DEC;
      else
        if (nx > 1 && ny > 1)
          decomposition = XY_DEC;
  
  if (decomposition == NO_DEC)
  {
    // All neighbour values are already -1;
  }
  
  if (decomposition == X_DEC)
  {
    if (ThisTask > 0)
      ln = ThisTask - 1;
    
    if (ThisTask < (NTasks-1))
      rn = ThisTask + 1;
  }
    
  if (decomposition == Y_DEC)
  {
    if (ThisTask > 0)
      tn = ThisTask - 1;
    
    if (ThisTask < (NTasks-1))
      bn = ThisTask + 1;
  }
    
  if (decomposition == XY_DEC)
  {
    row = ThisTask / px;
    col = ThisTask % px;
    
    // 
    // Set Left - Right neighbours    
    // 
    if (col > 0)
      ln = ThisTask - 1;
      
    if (col < px - 1)
      rn = ThisTask + 1;
      
    // 
    // Set Top - Bottom neighbours
    // 
    if (row > 0)
      tn = ThisTask - NTasks / py;
    
    if (row < py - 1)
      bn = ThisTask + NTasks / py;    
  }
  
  printf ("ThisTask %d  ln %d, rn %d, tn %d, bn %d\n", ThisTask, ln, rn, tn, bn );
  
  return 1;
}
