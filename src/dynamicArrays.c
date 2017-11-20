#include "allvars.h"
#include "proto.h"
 
void *my_malloc ( char *expr, size_t size ) {
    void *result = malloc( size );
    printf( "Malloc(%s) is size %lu, returning %p\n", expr, (unsigned long)size, result );
    return result;
}
void my_free( void *ptr ) {
    printf( "Free(%p)\n", ptr );
    free( ptr );
}
#define MY_MALLOC(x)    my_malloc( #x, x )
#define MY_FREE(x)      my_free(x)
 
//
//HAE Velocity V 3D array memory CONTIGUOUS allocation
//
int reserve_v (void)
{
   int max_x=Nx;
   int max_y=Ny;
   int max_z=Nz;
                    
                v = (REAL ***) malloc( max_x * sizeof *v );
   REAL   **all_y = (REAL **)  malloc( max_x * max_y * sizeof *all_y );
   REAL    *all_z = (REAL *)   malloc( max_x * max_y * max_z * sizeof *all_z );

   int i, j;

   for ( i = 0 ; i < max_x ; i++, all_y += max_y ) {
      v[i] = all_y;
      for ( j = 0 ; j < max_y ; j++, all_z += max_z ) {
          v[i][j] = all_z;
      }
   }


   return 1;
}


//
//HAE Pressure U t+3D (4D) array CONTIGUOUS memory allocation
//

int reserve_u (void)
{
   int max_t=3;
   int max_x=Nx+8;
   int max_y=Ny+8;
   int max_z=Nz+8;
    
                u = (REAL ****)malloc( max_t * sizeof *u );
   REAL  ***all_x = (REAL ***) malloc( max_t * max_x * sizeof *all_x );
   REAL   **all_y = (REAL **)  malloc( max_t * max_x * max_y * sizeof *all_y );
   REAL    *all_z = (REAL *)   malloc( max_t * max_x * max_y * max_z * sizeof *all_z );

   int t, i, j;

   for ( t = 0 ; t < max_t ; t++, all_x += max_x ) {
     u[t] = all_x;
     for ( i = 0 ; i < max_x ; i++, all_y += max_y ) {
         u[t][i] = all_y;
         for ( j = 0 ; j < max_y ; j++, all_z += max_z ) {
             u[t][i][j] = all_z;
         }
     }
   }
   
   return 1;
}
