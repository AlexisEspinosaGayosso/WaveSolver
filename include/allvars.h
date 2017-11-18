#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <omp.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>

#define NAME_LENGTH 256

#ifdef SINGLE

#define REAL    float
#define MPIREAL MPI_REAL

#else

#define REAL    double
#define MPIREAL MPI_DOUBLE

#endif
    //
    // Flag to enable - disable output
    //
    int write_output;
    
    // 
    // Length of each dimension
    // 
    REAL Lx;
    REAL Ly;
    REAL Lz;
    
    //
    // Time-step length
    //
    REAL dt;
    
    //
    // Size of cells in simulation
    //
    REAL dx;
    REAL dy;
    REAL dz;
    
    //
    // Factor that multiplies dt and dx
    // to make cell size and time-step
    // smaller or bigger by  fac times
    //
    REAL fac;
    
    //
    // Number of Timesteps to run simulation
    //
    int T;
    
    //
    //  Although these are redundant for Nx = Ny = Nz,
    //  its better to have these in case we wanted to
    //  expand the code for Nx != Ny != Nz
    //
    int Nx;
    int Ny;
    int Nz;
    
    //
    // Delta t ^2 and Delta x^2 computed
    // in order to reduce number of 
    // operations in the run
    //    
    REAL dt2;
    REAL dx2;
    
    //
    // Pression U 3D array memory allocation
    // 
    REAL **** u;
    
    //
    // Velocity V 3D array
    // 
    REAL *** v;
            
    //
    // File to which output will be written
    //
    FILE * f;

    //
    // Position Source
    //
    REAL Source_coord[3];
    int  Source_cells[3];
    
    int  min_cell[3];
    int  max_cell[3];
    
    int  SourceTask;
    
#ifdef MPI
  int    ThisTask;
  int    NTasks;
  
  
//
// Define types decompositions and 
// boundaries to make fast choice
// when communicating with neighbours
//

#define NO_DEC   0
#define  X_DEC   1
#define  Y_DEC   2
#define XY_DEC   3
  int row;
  int col;
  
  int decomposition; 
  
  int px;
  int py;
  
  int tn;
  int bn;
  int ln;
  int rn;
  
  int nx_to_sr;
  int ny_to_sr;
  
  REAL * utosend_l;
  REAL * utorecv_l;

  REAL * utosend_r;
  REAL * utorecv_r;

  REAL * utosend_t;
  REAL * utorecv_t;

  REAL * utosend_b;
  REAL * utorecv_b;

  MPI_Datatype UTOSENDLR; //HAE to be used for the send left to right and viceversa
  MPI_Datatype UTOSENDTB; //HAE to be used for the send up to bottom and viceversa


/*
 * 0 - No Decomposition
 *     px = 1
 *     py = 1
 * ________
 * |      |
 * |      |
 * |______|
 *
 * 
 */ 


/*
 * 1 - X Decomposition
 *     px > 1
 *     py = 1
 * ______________________
 * |      |      |      |
 * |  HL  |  HC  |  HR  |
 * |______|______|______|
 *
 */  
// #define   HL   10
// #define   HC   11
// #define   HR   12

/*  
 * 2 - Y Decomposition
 *     px = 1
 *     py > 1
 * ________  
 * |      |  
 * |  VT  |  
 * |______|  
 * |      |  
 * |  VC  |  
 * |______|  
 * |      |  
 * |  VB  |  
 * |______|
 *  
 */
// #define   VT   20
// #define   VC   21
// #define   VB   22

/*
 * 3 - XY Decomposition
 *     px > 1
 *     py > 1
 * ______________________
 * |      |      |      |
 * |  TL  |  TC  |  TR  |
 * |______|______|______|
 * |      |      |      |
 * |  CL  |  CC  |  CR  |
 * |______|______|______|
 * |      |      |      |
 * |  BL  |  BC  |  BR  |
 * |______|______|______|
 *
 */  
// #define   TL   0  
// #define   TC   1  
// #define   TR   2  
// #define   CL   3  
// #define   CC   4  
// #define   CR   5  
// #define   BL   6  
// #define   BC   7  
// #define   TR   8 
  
#endif




#endif
