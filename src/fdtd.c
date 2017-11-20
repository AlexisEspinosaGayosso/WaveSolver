#include "allvars.h"
#include "proto.h"


 /*
  * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  *                                                                                 *
  *                     Wave Propagation in 3D using FDTD                           *
  *                                                                                 *
  * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  *   
  * 
  * 1/\Delta t^2 \sum_{s=l_time}^{r_time} D_t u(x,y,z,t + s\Delta t) =
  * 
  *     v(x,y,z)^2 [ 1/\Delta x^2 \sum_{g=l_space}^{r_space} D_g u(x + g\Delta x,y,z,t) +
  * 
  *                  1/\Delta y^2 \sum_{h=l_space}^{r_space} D_h u(x,y + h\Delta y,z,t) +
  * 
  *                  1/\Delta z^2 \sum_{k=l_space}^{r_space} D_k u(x,y,z + k\Delta z,t) ]
  * 
  *    
  * 
  * Solving for   u(x,y,z,t + r_time \Delta t)
  * 
  *    u(x,y,z,t + r_time \Delta t) =
  * 
  *            v(x,y,z)^2 \Delta t^2 / ( D_t u(x,y,z,t + l_time\Delta t + D_t u(x,y,z,t + l_time\Delta t ) \times
  *             
  *                [ 1/\Delta x^2 \sum_{g=l_space}^{r_space} D_g u(x + g\Delta x,y,z,t) +
  * 
  *                  1/\Delta y^2 \sum_{h=l_space}^{r_space} D_h u(x,y + h\Delta y,z,t) +
  * 
  *                  1/\Delta z^2 \sum_{k=l_space}^{r_space} D_k u(x,y,z + k\Delta z,t) ]
  * 
  * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *
  * 
  *   wave_prop_3D_cube ()
  *
  *     For general case where dx != dy != dz
  *        
  *         sum_x = 0;    sum_y = 0;    sum_z = 0;
  *         for (l=0, i=-4; l < 9; l++,i++)
  *         {
  *              sum_x += Dk[l] * u[x+i][y][z];
  *              sum_y += Dk[l] * u[x][y+i][z];
  *              sum_z += Dk[l] * u[x][y][z+i];
  *         }
  *         u[x][y][z] = v2 * dt2 * ((sum_x / dx2) + (sum_y / dy2) + (sum_z / dz2))
  * 
  * 
  */

int wave_prop_3D_cube (void)
{
    
    // 
    // Variables used for loops
    // 
    int i, j, k, l;
    int x, y, z, t;
    
    //
    // File to which output will be written
    //
    if (write_output == 1)
    {
      char name[100];
#ifdef MPI
      sprintf(name,"output_%d", ThisTask);
#else
      sprintf(name,"output");
#endif
      f = fopen (name,"w");
    }
    
    //
    // Factors for spatial differentiation
    // 
    REAL Dk [9] = {-1.0/560, 8.0/315, -1.0/5, 8.0/5, -205.0/72, 8.0/5, -1.0/5, 8.0/315, -1.0/560};
    
    //
    // Factor for time differentiation
    //
    REAL Dt [3] = {1.0, -2.0, 1.0};
    
    //
    // Integers used for evolution of
    // simulation
    //
    int t_l, t_0, t_r;
    
    //
    // Variables for source term
    //
    REAL v2;                                 // square of velocity
    REAL sqrt2pi = sqrt(2.0 * acos(-1.0));   // square root (2 * pi )
    REAL sigma = 0.011;                      // Sigma value for source term
    REAL sigma2 = sigma * sigma;             // sigma^2 to reduce # of operations
    REAL sigma3 = sigma2 * sigma;            // sigma^3 ro reduce # of operations 
    REAL c1 = sqrt2pi * sigma3;              // constant sqrt(2*pi)*sigma^3
    REAL t0 = 0.1;
    REAL t_t02;                              // ( t - t_0 )^2
    REAL ft;                                 // source term
    
    //
    // Variables for Summatory of spatial
    // and time Finite-differentiation
    //
    REAL sum_p, sum_t;
    
    //
    // x,y,z position of Source in 
    // Grid coordinates
    //    
    int Srcx = Source_cells[0] + 4;
    int Srcy = Source_cells[1] + 4;
    int Srcz = Source_cells[2] + 4;
    
    //
    // NOTE: This implementation only works when dx = dy = dz
    //       and for Nx = Ny = Nz = N

    //    
    // Time Loop
    //
    
    int block_size=9; //Nx+4; //HAE Testing the loop blocking
    for (t = 0, t_l = 0, t_0 = 1, t_r = 2; t < T; t++) // Time Loop
    {    
      //for (x = 4; x < Nx + 4; x++)  // Loop over x coordinate
      for (int xx = 4; xx < Nx + 4; xx+=block_size)  // Loop over x coordinate
      {
	      //for (y = 4; y < Ny + 4; y++)  // Loop over y coordinate
         for (int yy = 4; yy < Ny + 4; yy+=block_size)  // Loop over y coordinate
      	{
	         //for (z = 4; z < Nz + 4; z++)  // Loop over z coordinate
            for (int zz = 4; zz < Nz + 4; zz+=block_size)  // Loop over z coordinate
	         {
               for (x=xx; x<xx+block_size; x++){for (y=yy; y<yy+block_size; y++){for (z=zz; z<zz+block_size; z++){
	            v2 = v[x-4][y-4][z-4] * v[x-4][y-4][z-4];  // v^2 of current cell
	    
	            sum_p = 0;                                 // set Summatory to 0
	    
	            ft = 0;                                    // set source to 0

               //
               // If source is located at current cell
               // add source term
               //
               if ((x == Srcx) && (y == Srcy) && (z == Srcz)) 
               {
                 t_t02 =  (t+1)*dt - t0;
                 t_t02 *= t_t02;
                 ft = (1 - t_t02/sigma2) / c1 * exp(-t_t02/(2*sigma2));
               }
               
               //
               // Spatial Differentiation
               //
               /*AEG_optimizando for (l=0, i=-4; l < 9; l++,i++)
                 sum_p += Dk[l] * (u[t_0][x+i][y][z] \
                                 + u[t_0][x][y+i][z] \
                                 + u[t_0][x][y][z+i]);*/

               for (l=0, i=-4; l < 9; l++,i++)
                 sum_p += Dk[l] * ( 0 
                                 + u[t_0][x+i][y][z] 
                                 + u[t_0][x][y+i][z] 
                                 + u[t_0][x][y][z+i] 
                                  );
               //sum_p += 0
                   //     +  Dk[0] * u[t_0][x-4][y][z] + Dk[1] * u[t_0][x-3][y][z] + Dk[2] * u[t_0][x-2][y][z] + Dk[3] * u[t_0][x-1][y][z] + Dk[4] * u[t_0][x][y][z] 
                   //     + Dk[5] * u[t_0][x+1][y][z] + Dk[6] * u[t_0][x+2][y][z] + Dk[7] * u[t_0][x+3][y][z] + Dk[8] * u[t_0][x+4][y][z] 
                   //     + Dk[0] * u[t_0][x][y-4][z] + Dk[1] * u[t_0][x][y-3][z] + Dk[2] * u[t_0][x][y-2][z] + Dk[3] * u[t_0][x][y-1][z] + Dk[4] * u[t_0][x][y][z] 
                   //     + Dk[5] * u[t_0][x][y+1][z] + Dk[6] * u[t_0][x][y+2][z] + Dk[7] * u[t_0][x][y+3][z] + Dk[8] * u[t_0][x][y+4][z] 
                   //     + Dk[0] * u[t_0][x][y][z-4] + Dk[1] * u[t_0][x][y][z-3] + Dk[2] * u[t_0][x][y][z-2] + Dk[3] * u[t_0][x][y][z-1] + Dk[4] * u[t_0][x][y][z] 
                   //     + Dk[5] * u[t_0][x][y][z+1] + Dk[6] * u[t_0][x][y][z+2] + Dk[7] * u[t_0][x][y][z+3] + Dk[8] * u[t_0][x][y][z+4] 
                   //     ;

               /*AEG_optimizando sum_p += \
                          Dk[0] * u[t_0][x-4][y][z] + Dk[1] * u[t_0][x-3][y][z] + Dk[2] * u[t_0][x-2][y][z] + Dk[3] * u[t_0][x-1][y][z] + Dk[4] * u[t_0][x][y][z] \
                        + Dk[5] * u[t_0][x+1][y][z] + Dk[6] * u[t_0][x+2][y][z] + Dk[7] * u[t_0][x+3][y][z] + Dk[8] * u[t_0][x+4][y][z] \
                        + Dk[0] * u[t_0][x][y-4][z] + Dk[1] * u[t_0][x][y-3][z] + Dk[2] * u[t_0][x][y-2][z] + Dk[3] * u[t_0][x][y-1][z] + Dk[4] * u[t_0][x][y][z] \
                        + Dk[5] * u[t_0][x][y+1][z] + Dk[6] * u[t_0][x][y+2][z] + Dk[7] * u[t_0][x][y+3][z] + Dk[8] * u[t_0][x][y+4][z] \
                        + Dk[0] * u[t_0][x][y][z-4] + Dk[1] * u[t_0][x][y][z-3] + Dk[2] * u[t_0][x][y][z-2] + Dk[3] * u[t_0][x][y][z-1] + Dk[4] * u[t_0][x][y][z] \
                        + Dk[5] * u[t_0][x][y][z+1] + Dk[6] * u[t_0][x][y][z+2] + Dk[7] * u[t_0][x][y][z+3] + Dk[8] * u[t_0][x][y][z+4] \
                        ;
                */

               //
               // Time Differentiation
               //  
               sum_t = Dt[0]*u[t_l][x][y][z]  +  Dt[1]*u[t_0][x][y][z];
               
               //
               // Compute u(x,y,z,t+Dt)
               //
               u[t_r][x][y][z] = v2 * dt2 * (ft + sum_p/dx2) - sum_t;
               }}}
            }
         }
      }
      if ((t+1)%10 == 0)            // Output written every 10 time-steps
      {
       if (write_output == 1)
#ifdef OPENMP
         if(omp_get_thread_num() == 0) // Only master thread prints results
#endif         
         {
	        for (x = 4; x < Nx + 4; x++)
	        {
	           for (y = 4; y < Ny + 4; y++)
	           fprintf (f, "%6.5f  ", u[t_0][x][y][Nz/2 + 4]);
	           fprintf (f, "\n");
	        }
         }

	    printf ("time = %g\n", t*dt);
      }
      t_l = (t_l + 1) % 3;
      t_0 = (t_0 + 1) % 3;
      t_r = (t_r + 1) % 3;
      

#ifdef MPI
      MPI_Barrier (MPI_COMM_WORLD);
      communicate (t_0);
#endif
  }


  if (write_output == 1)
    fclose(f);

  //
  // Free U 3-3D array
  //   
/*AEG_textBook  my_free(u[0][0][0]);
  my_free(u[0][0]);
  my_free(u[0]);
  my_free(u);
*/
  printf("u array freed\n");
  //
  // Free V 3D array
  // 
/*AEG_quita  for (j = 0; j < Nx; j++)
  {
    for (k = 0; k < Ny; k++)
    {
      free(v[j][k]);
    }
    free(v[j]);
  }
  free(v);
*/

/*AEG_textBook  my_free(v[0][0]);
  my_free(v[0]);
  my_free(v);
*/
  printf("v array freed\n");
  
  return 1;
}
