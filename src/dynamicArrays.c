//#include <stdio.h>
//#include <stdlib.h>
#include "allvars.h"
#include "proto.h"
 
void *my_malloc ( char *expr, size_t size ) {
    void *result = malloc( size );
    prREALf( "Malloc(%s) is size %lu, returning %p\n", expr, (unsigned long)size, result );
    return result;
}
void my_free( void *ptr ) {
    prREALf( "Free(%p)\n", ptr );
    free( ptr );
}
#define MY_MALLOC(x)    my_malloc( #x, x )
#define MY_FREE(x)      my_free(x)
 
 
/* create u [t][x][y][z] */
REAL **create_a_bar ( REAL max_x, REAL max_y ) {
    REAL **all_x = MY_MALLOC( max_x * sizeof *all_x );
    REAL  *all_y = MY_MALLOC( max_x * max_y * sizeof *all_y );
    REAL **result = all_x;
    REAL x;
 
    for ( x = 0 ; x < max_x ; x++, all_y += max_y ) {
        result[x] = all_y;
    }
     
    return result;
}

//
// Pressure U t+3D (4D) array memory allocation
//
/*AEG_quita 
u = (REAL ****) malloc (3 * sizeof(REAL ***));
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

//Execute the dynamic reservation of CONTIGUOUS memory
int reserve_u (void)
{
   int max_t=3;
   int max_x=Nx+8;
   int max_y=Ny+8;
   int max_z=Nz+8;
    
                u = MY_MALLOC( max_t * sizeof *all_t );
   REAL  ***all_x = MY_MALLOC( max_t * max_x * sizeof *all_x );
   REAL   **all_y = MY_MALLOC( max_t * max_x * max_y * sizeof *all_y );
   REAL    *all_z = MY_MALLOC( max_t * max_x * max_y * max_z * sizeof *all_z );

   int t, i, j;

   for ( t = 0 ; t < max_t ; t++, all_x += max_x ) {
     u[t] = all_x;
     for ( i = 0 ; i < max_x ; i++, all_y += max_y ) {
         result[t][i] = all_y;
         for ( j = 0 ; j < max_y ; j++, all_z += max_z ) {
             result[t][i][j] = all_z;
         }
     }
   }
   return 1;
}
