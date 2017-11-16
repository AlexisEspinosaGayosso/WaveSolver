#include "allvars.h"
#include "proto.h"

int communicate (int t0)
{
  MPI_Request rqst;
  MPI_Status  status;
  
  int i, j, k, l, m;
  
  int offset;
  
  //
  // Information to be sent, is copied to a 1-D array to
  // keep the number of messages to the minimum.
  // 
  // The message received is then copied to the corresponding
  // 3-D position.
  //
  
  //
  // Left Right Send
  //
  if (ln > -1)
  {
    offset = 0;
    
    for (i = 4; i < 8; i++)
      for (j = 4; j < Ny + 4; j++)
        for (k = 4; k < Nz + 4; k++)
          utosend_l[offset++] = u[t0][i][j][k];

    MPI_Isend (utosend_l, nx_to_sr, MPIREAL, ln, ThisTask, MPI_COMM_WORLD, &rqst);
  }

  if (rn > -1)
  {
    offset = 0;
    
    for (i = Nx; i < Nx + 4; i++)
      for (j = 4; j < Ny + 4; j++)
        for (k = 4; k < Nz + 4; k++)
          utosend_r[offset++] = u[t0][i][j][k];

    MPI_Isend (utosend_r, nx_to_sr, MPIREAL, rn, ThisTask, MPI_COMM_WORLD, &rqst);    
  }
  
  //
  // Top-Bottom Send
  //
  if (tn > -1)
  {
    offset = 0;
    
    for (i = 4; i < Nx + 4; i++)
      for (j = 4; j < 8; j++)
        for (k = 4; k < Nz + 4; k++)
          utosend_t[offset++] = u[t0][i][j][k];

    MPI_Isend (utosend_t, ny_to_sr, MPIREAL, tn, ThisTask, MPI_COMM_WORLD, &rqst);
 }
  
  if (bn > -1)
  {
    offset = 0;
    
    for (i = 4; i < Nx + 4; i++)
      for (j = Ny; j < Ny + 4; j++)
        for (k = 4; k < Nz + 4; k++)
          utosend_b[offset++] = u[t0][i][j][k];

    MPI_Isend (utosend_b, ny_to_sr, MPIREAL, bn, ThisTask, MPI_COMM_WORLD, &rqst);
  }
  
//   printf ("Task %d - Done Sending\n", ThisTask);
  MPI_Barrier (MPI_COMM_WORLD);  
  
  //
  // Left-Right Recieve
  //
  if (ln > -1)
  {
    MPI_Recv (utorecv_l, nx_to_sr, MPIREAL, ln, ln, MPI_COMM_WORLD, &status);
    
    offset = 0;
    
    for (i = 0; i < 4; i++)
      for (j = 4; j < Ny + 4; j++)
        for (k = 4; k < Nz + 4; k++)
          u[t0][i][j][k] = utorecv_l[offset++];
  }

  if (rn > -1)
  {
    MPI_Recv (utorecv_r, nx_to_sr, MPIREAL, rn, rn, MPI_COMM_WORLD, &status);
    
    offset = 0;
    
    for (i = Nx + 4; i < Nx + 8; i++)
      for (j = 4; j < Ny + 4; j++)
        for (k = 4; k < Nz + 4; k++)
          u[t0][i][j][k] = utorecv_r[offset++];
  }
  
  //
  // Top-Bottom Receive
  //
  if (tn > -1)
  {
    MPI_Recv (utorecv_t, ny_to_sr, MPIREAL, tn, tn, MPI_COMM_WORLD, &status);
    
    offset = 0;
    
    for (i = 4; i < Nx + 4; i++)
      for (j = 0; j < 4; j++)
        for (k = 4; k < Nz + 4; k++)
          u[t0][i][j][k] = utorecv_t[offset++];
  }
  
  if (bn > -1)
  {
    MPI_Recv (utorecv_b, ny_to_sr, MPIREAL, bn, bn, MPI_COMM_WORLD, &status);
    
    offset = 0;
    
    for (i = 4; i < Nx + 4; i++)
      for (j = Ny + 4; j < Ny + 8; j++)
        for (k = 4; k < Nz + 4; k++)
          u[t0][i][j][k] = utorecv_b[offset++];
  } 
  
  MPI_Barrier (MPI_COMM_WORLD);
  
  return 1; 
}