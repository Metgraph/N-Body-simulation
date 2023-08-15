#ifndef __BH_H_
#define __BH_H_

#include <cuda_runtime_api.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    int firstfree;
    int ff_mutex;
    int root;
    double max;
    double4 *node_pos;
    int *ents;
    int *parent;
    int *children;
    int *mutex;
    int sz;
} Octtree;


__host__
void cuda_check_error(cudaError_t error, const char *message);

// Tree
__global__
void set_tree(Octtree *tree, int ents_sz, double4 *pos, int *ents, int *parent, int *children);
__host__
void init_tree(Octtree **tree, int ents_sz, double4 **pos, int **ents, int **parent, int **children, int **mutex);
__host__
void free_tree(Octtree *tree);
__global__
void print_max_tree(Octtree *tree);
__device__
int init_node(Octtree *tree);
__host__
void reset_mutex(int *mutex, int ents_sz);


/*
__host__ void get_entities(char filename[], uint *n_ents, double **positions, double **velocities);
__host__ void count_entities_file(char *filename, uint *n);
__global__ void acceleration(double *positions, double *acc, uint ents_sz, int step);
__host__ void safe_malloc_double(double **pointer, int quantity);
__global__ void update_positions(double *positions, double *velocities, uint ents_sz, double dt, int step);
__global__ void update_velocities(double *acc, double *vel, uint ents_sz, double dt);
__host__ void run(double *positions, double *velocities, uint ents_sz, int n_steps, double dt, const char *output);
__global__ void get_bounding_box_fi(double *positions, int ents_sz, double *g_odata, int step);

__device__ void print_tree_rec(Octtree *tree, int id, char *space, uint depth);
__global__ void cuda_print(Octtree *tree, int ents_sz, char *space);
__host__ void print_tree(Octtree *tree, int ents_sz);
__global__ void get_bounding_box_fi(double *positions, int ents_sz, double *g_odata, int step);
__global__ void get_bounding_box(double *g_idata, int ents_sz, double *g_odata);
__global__ void set_max(double *max, Octtree *tree);
__host__ void bounding_box(double *positions, double *max_in, double *max_out, int ents_sz, int step, Octtree *tree);
__device__ int get_indx_loc(double4 *pos, double3 *center, double *border_size);
__device__ void set_lock(Octtree *tree, int index_mutex);
__device__ void unset_lock(Octtree *tree, int index_mutex);
__global__ void add_ent(Octtree *tree, double4 *ents, int ents_sz);
*/







#endif

