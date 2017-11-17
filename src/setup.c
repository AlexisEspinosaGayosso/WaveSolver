#include "allvars.h"
#include "proto.h"

int setup (void)
{
    // 
    // Variables used for loops
    // 
    int i, j, k, l;
    int x, y, z, t;
  
    // 
    // Length of each dimension
    // 
    Lx  = 3000;            // 3 km
    Ly  = 3000;            // 3 km
    Lz  = 3000;            // 3 km
    
    //
    // Time-step length
    //
    dt  = 0.001;           // in seconds
    
    //
    // Size of cells in simulation
    //
    dx  = 12.5;            // in meters
    dy  = 12.5;            // in meters
    dz  = 12.5;            // in meters
    
    //
    // Factor that multiplies dt and dx
    // to make cell size and time-step
    // smaller or bigger by  fac times
    //
//     fac = 1.0;
    
    dt *= fac;
    dx *= fac;
    dy *= fac;
    dz *= fac;
 
    //
    // Number of Timesteps to run simulation
    //
//     T = 1000;
    
    //
    // Delta t ^2 and Delta x^2 computed
    // in order to reduce number of 
    // operations in the run
    //    
    dt2 = dt * dt; 
    dx2 = dx * dx;

    
    //
    // Number of cells per dimension
    //
    Nx = (int) Lx / dx / px;
    Ny = (int) Ly / dy / py;
    Nz = (int) Lz / dz;

    min_cell[0] = 0;
    min_cell[1] = 0;
    min_cell[2] = 0;

    max_cell[0] = Nx;// - 1;
    max_cell[1] = Ny;// - 1;
    max_cell[2] = Nz;// - 1;
    
    
    if (decomposition == NO_DEC)
    {
      // All neighbour values are already -1;
    }
    
    if (decomposition == X_DEC)
    {
      min_cell[0] += ThisTask * Nx;
      max_cell[0] += ThisTask * Nx;
    }
      
    if (decomposition == Y_DEC)
    {
      min_cell[1] += ThisTask * Ny;
      max_cell[1] += ThisTask * Ny;
    }

    if (decomposition == XY_DEC)
    {
      min_cell[0] += Nx * (ThisTask % px);
      max_cell[0] += Nx * (ThisTask % px);
      
      min_cell[1] += Ny * (ThisTask / px);
      max_cell[1] += Ny * (ThisTask / px);
    }
    
    //
    // Source position
    //
    Source_coord[0] = 1500;
    Source_coord[1] = 1500;
    Source_coord[2] = 1500;
    
    Source_cells[0] = (int) Source_coord[0] / dx;
    Source_cells[1] = (int) Source_coord[1] / dy;
    Source_cells[2] = (int) Source_coord[2] / dz;
    
    int tmp = -1;
    
    if (Source_cells[0] < max_cell[0] && Source_cells[0] >= min_cell[0])
      if (Source_cells[1] < max_cell[1] && Source_cells[1] >= min_cell[1])
        tmp = ThisTask;
    
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Allreduce (&tmp, &SourceTask, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);       
    
    if (SourceTask == ThisTask)
    {
      Source_cells[0] = Source_cells[0] % (Nx);
      Source_cells[1] = Source_cells[1] % (Ny);
      Source_cells[2] = Source_cells[2] % (Nz);
    }
    else
    {
      Source_cells[0] = -4;
      Source_cells[1] = -4;
      Source_cells[2] = -4;      
    }
    
    nx_to_sr = Ny * Nz * 4;
    ny_to_sr = Nx * Nz * 4;
    
    if (ln > -1)
    {
      utosend_l = (REAL * ) malloc (nx_to_sr * sizeof(REAL)); 
      utorecv_l = (REAL * ) malloc (nx_to_sr * sizeof(REAL)); 
    }
    
    if (rn > -1)
    {
      utosend_r = (REAL * ) malloc (nx_to_sr * sizeof(REAL)); 
      utorecv_r = (REAL * ) malloc (nx_to_sr * sizeof(REAL)); 
    }
    
    if (tn > -1)
    {
      utosend_t = (REAL * ) malloc (ny_to_sr * sizeof(REAL));
      utorecv_t = (REAL * ) malloc (ny_to_sr * sizeof(REAL));
    }
  
    if (bn > -1)
    {
      utosend_b = (REAL * ) malloc (ny_to_sr * sizeof(REAL));
      utorecv_b = (REAL * ) malloc (ny_to_sr * sizeof(REAL));
    } 
    
    //
    // Pression U 3D array memory allocation
    // 
/*AEG quita    u = (REAL ****) malloc (3 * sizeof(REAL ***));
    for (i = 0; i < 3; i++)
    {
      u[i] = (REAL ***) malloc ((Nx+8)*sizeof(REAL **));
      for (j = 0; j < Nx+8; j++)
      {
	     u[i][j] = (REAL **) malloc ((Ny+8)*sizeof(REAL *));
        for (k = 0; k < Ny+8; k++)
	       u[i][j][k] = (REAL *) malloc ((Nz+8)*sizeof(REAL));
      }
    }
*/
   reserve_u(); //Execute the dynamic reservation of CONTIGUOUS memory
    //
    // Velocity V 3D array memory allocation
    // 
    v = (REAL ***) malloc (Nx*sizeof(REAL **));
    for (j = 0; j < Nx; j++)
    {
      v[j] = (REAL **) malloc (Ny*sizeof(REAL *));
      for (k = 0; k < Ny; k++)
	     v[j][k] = (REAL *) malloc (Nz*sizeof(REAL));
    }
    
    //
    // Initialization of pressure to zero in all grid
    //
    for (t = 0; t < 3; t++)
      for (x = 0; x < (Nx + 8); x++)
	     for (y = 0; y < (Ny + 8); y++)
   	    for (z = 0; z < (Nz + 8); z++)
            u[t][x][y][z] = 0;
    
          
    //
    // Initialization of Velocity 
    // 
    for (x = 0; x < Nx; x++)
      for (y = 0; y < Ny; y++)
        for (z = 0; z < Nz; z++)
          if ( z <= (Nz/6))      //--- Set velocity to 1.5 km/s water layer
            v[x][y][z] = 1500;   
          else
            v[x][y][z] = 2000;   //--- Set velocity to 2.0 km/s sediment layer
	    
    
    //
    // Initialization of velocity in Homogeneous medium
    //
//     for (x = 0; x < Nx; x++)
//       for (y = 0; y < Ny; y++)
//         for (z = 0; z < Nz; z++)
//           v[x][y][z] = 3000;
    
    return 1;  
}
