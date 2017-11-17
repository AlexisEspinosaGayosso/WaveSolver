#ifndef PROTO_H_
#define PROTO_H_

int wave_prop_3D_cube (void);

int get_neighbours (int nx, int ny);

int setup (void);

void *my_malloc (char *expr, size_t size);

void my_free( void *ptr);

int communicate (int t0);

int comm_test (void);

int wave_prop_3D_cube_v2 (void);

#endif
