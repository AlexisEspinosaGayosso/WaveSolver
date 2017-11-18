#include "allvars.h"
#include "proto.h"

int main (int argc, char ** argv)
{
  
#ifdef MPI  
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &NTasks);
  MPI_Comm_rank (MPI_COMM_WORLD, &ThisTask);
  
  px = atoi (argv[4]);
  py = atoi (argv[5]);
  
  if (px*py != NTasks || px < 1 || py < 1)
  {
    printf ("Error px*py is different to Ntasks\n");
    printf ("Exiting\n");    
    MPI_Finalize();
    exit(0);
  }
#endif 

  
  //
  // Get factor to multiply Dx and Dt
  //
  fac = atof(argv[1]);
  
  //
  // Number of time steps to run the simulation
  //
  T   = atoi(argv[2]);
  
  //
  // Output writing
  //
  write_output = atoi (argv[3]);

  //
  // Decomposition
  // 
  get_neighbours (px, py);
  printf ("Task %d - Done neighbours\n", ThisTask);

  
  //
  // Setup Simulation
  //
  setup ();
  printf ("Task %d - Done setup\n", ThisTask);
  
  
  //
  // Communication Test
  //
//   comm_test ();
//   printf ("Task %d - Done communication\n", ThisTask);


  //
  // Start clock for wall-time 
  // timing
  //
  struct timeval start, finish;
  gettimeofday (&start, NULL);
  printf ("size of real  %ld\n", sizeof(REAL));

  clock_t starttime = clock();
  
  //
  // MainLoop
  // 
  printf ("Calling loop\n");
  wave_prop_3D_cube ();

  //
  // Measure wall-time 
  //  
  gettimeofday (&finish, NULL);
  double total_time  = (finish.tv_sec - start.tv_sec);
  total_time += (finish.tv_usec - start.tv_usec) / 1000000.0;
  printf ("Total time %lf seconds\n",total_time);

  clock_t finaltime = clock();
  double thetime = ((double)finaltime-starttime) / CLOCKS_PER_SEC;
  printf ("The time %lf seconds\n",thetime);


  // 
  // Terminate
  //  
#ifdef MPI
  MPI_Type_free(&UTOSENDLR);
  MPI_Type_free(&UTOSENDTB);
  MPI_Finalize ();
#endif

  
  return 0;
}
