#include <stdio.h>
#include "barnes-hut.cuh"

#define PRE_ALLOC_SIZE 5

__global__
void set_tree(Octtree *tree, int ents_sz, double4 *pos, int *ents, int *parent, int *children, int *mutex) {
    tree->firstfree = ents_sz;
    tree->ff_mutex = 0;
    tree->root = ents_sz;
    tree->node_pos = pos;
    tree->ents = ents;
    tree->parent = parent;
    tree->children = children;
    tree->mutex = mutex;
    tree->sz = ents_sz * PRE_ALLOC_SIZE;

    init_node(tree);
}

__host__
void init_tree(Octtree **tree, int ents_sz, double4 **pos, int **ents, int **parent, int **children, int **mutex) {
    cudaError_t error;

    error = cudaMalloc((void **)tree, sizeof(Octtree));
    cuda_check_error(error, "Tree malloc\n");

    error = cudaMalloc((void **)pos, ents_sz * sizeof(double4) * PRE_ALLOC_SIZE); // Posizione dei nodi
    cuda_check_error(error, "Tree positions malloc\n");

    error = cudaMalloc((void **)ents, ents_sz * sizeof(int) * PRE_ALLOC_SIZE);
    cuda_check_error(error, "Tree ents malloc\n");

    error = cudaMalloc((void **)parent, ents_sz * sizeof(int) * PRE_ALLOC_SIZE);
    cuda_check_error(error, "Tree parent malloc\n");

    error = cudaMalloc((void **)children, ents_sz * sizeof(int) * 8 * PRE_ALLOC_SIZE);
    cuda_check_error(error, "Tree children malloc\n");

    error = cudaMalloc((void **)mutex, ents_sz * sizeof(int) * 8 * PRE_ALLOC_SIZE);
    cuda_check_error(error, "Tree mutex malloc\n");

    error = cudaMemset(*mutex, 0, ents_sz * sizeof(int) * 8 * PRE_ALLOC_SIZE);
    cuda_check_error(error, "Tree mutex memset\n");

    set_tree<<<1, 1>>>(*tree, ents_sz, *pos, *ents, *parent, *children, *mutex);
    cudaDeviceSynchronize();
}

__host__
void reset_mutex(int *mutex, int ents_sz) {
    cudaError_t error;
    error = cudaMemset(mutex, 0, ents_sz * sizeof(int) * 8 * PRE_ALLOC_SIZE);
    cuda_check_error(error, "Reset mutex memset\n");
}

__host__
void free_tree(Octtree **tree, int ents_sz, double4 **pos, int **ents, int **parent, int **children, int **mutex) {
    cudaFree(*tree);
    cudaFree(*pos);
    cudaFree(*ents);
    cudaFree(*parent);
    cudaFree(*children);
    cudaFree(*mutex);
}

__global__
void print_max_tree(Octtree *tree) {
    printf("Max into tree: %fl\n", tree->max);
    printf("Original max into tree: %fl\n", tree->max / 2);
}

__device__
int init_node(Octtree *tree) {
    int new_node;
    int blocked = 1;
    while (blocked) {
        if (atomicCAS(&tree->ff_mutex, 0, 1) == 0){
            __threadfence();

            if (tree->firstfree >= tree->sz) {
                printf("ERROR: The space into the nodes array is finished. Restart with a bigger size.\n");
            }

            new_node = tree->firstfree++;

            atomicExch(&tree->ff_mutex, 0);
            __threadfence();
            blocked = 0;
        }
    }

    tree->node_pos[new_node].x = 0.0;
    tree->node_pos[new_node].y = 0;
    tree->node_pos[new_node].z = 0;
    tree->node_pos[new_node].w = 0; // Mass
    tree->ents[new_node] = 0;
    tree->parent[new_node] = -1;

    for (uint i = 0; i < 8; i++) {
        tree->children[new_node * 8 + i] = -1;
    }

    return new_node;
}
