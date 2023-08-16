// https://en.wikipedia.org/wiki/N-body_simulation
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "barnes-hut.cuh"

typedef unsigned int uint;

#define BILLION 1000000000.0

//__constant__ double BIG_G = 6.67e-11; // Nella read-only cache e condivisa dai thread del blocco
__constant__ double BIG_G = 1.0; // Nella read-only cache e condivisa dai thread del blocco

__host__ void get_entities(char filename[], uint *n_ents, double **positions, double **velocities);
__host__ void count_entities_file(char *filename, uint *n);
__global__ void acceleration(double *positions, double *acc, uint ents_sz, int step);
__host__ void safe_malloc_double(double **pointer, int quantity);
__global__ void update_positions(double *positions, double *velocities, uint ents_sz, double dt, int step);
__global__ void update_velocities(double *acc, double *vel, uint ents_sz, double dt);
__host__ void run(double *positions, double *velocities, uint ents_sz, int n_steps, double dt, const char *output);
__global__ void get_bounding_box_fi(double *positions, int ents_sz, double *g_odata, int step);

uint grid_s;
uint block_s;


__host__
void cuda_check_error(cudaError_t error, const char *message) {
    if (error != cudaSuccess) {
        fprintf(stderr, "Error: %s\nLine: %s", cudaGetErrorString(error), message);
        exit(EXIT_FAILURE);
    }
}

__host__
void print_tree_rec(double4 *node_pos, int *children, int *ents, int id, char *space, uint depth, double max, int child_i) {
    //Octnode *node = &tree->nodes[id];
    // how much divide
    uint temp = 1 << depth;
    double border = (max) / (double)temp;
    printf("%sid: %d, (x:%lf, y:%lf, z:%lf), mass:%lf, border: %lf\n", space,
           id, node_pos[id].x, node_pos[id].y, node_pos[id].z, node_pos[id].w, border);
    if (ents[id] > 1) {

        int i;
        for (i = depth * 4; i < depth * 4 + 4; i++) {
            space[i] = ' ';
        }
        space[i] = '\0';
        for (int i = 0; i < 8; i++) {
            if (children[id * 8 + i] > -1) {
                print_tree_rec(node_pos, children, ents, children[id*8+i], space, depth + 1, max, i);
            }
        }
        space[depth * 4] = '\0';
    }
}

// used for debug
__host__
void print_tree(Octtree *d_tree, int ents_sz, double4 *d_node_pos, int *d_children, int *d_ents) {
    uint sz_space = 4 * 40;
    char *space = (char *)malloc(sz_space * sizeof(char));
    space[0] = '\0';


    double4 *node_pos = (double4 *)malloc(ents_sz * sizeof(double4) * 5);
    int *children = (int *) malloc(ents_sz * sizeof(int) * 8 * 5);
    int *ents = (int *) malloc(ents_sz * sizeof(int) * 5);
    Octtree *tree = (Octtree *) malloc(sizeof(Octtree));

    cudaMemcpy(node_pos, d_node_pos, ents_sz * 5 * sizeof(double4), cudaMemcpyDeviceToHost);
    cudaMemcpy(children, d_children, ents_sz * 5 * sizeof(int) * 8, cudaMemcpyDeviceToHost);
    cudaMemcpy(ents, d_ents, ents_sz * 5 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(tree, d_tree, sizeof(Octtree), cudaMemcpyDeviceToHost);

    print_tree_rec(node_pos, children, ents, ents_sz, space, 0, tree->max, -1);

    free(node_pos);
    free(children);
    free(ents);
    free(tree);
    free(space);
}

__global__
void print_local(double *val, int ents_sz){
    double max = val[0];
    int pos = 0;
    for (int i=0; i<ents_sz; i++) {
        //printf("Val pos %d: %lf\n", i, val[i]);
        if (val[i] > max){
            max = val[i];
            pos=i;
        }
    }
    printf("Find max: %lf, %lf, POS: %d\n", max, max*2, pos);
}


__global__
void get_bounding_box_fi(double *positions, int ents_sz, double *g_odata, int step)
{
    extern __shared__ double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint tid = threadIdx.x;
    double4 *g_idata;
    int myId;
    double max;

    // Move pointer to current time step
    positions = positions + (ents_sz * 4 * step);
    g_idata = (double4 *) positions;

    myId = blockIdx.x * blockDim.x + threadIdx.x;
    max = 0.0;

    if (myId < ents_sz) {
        double4 body = g_idata[myId];
        max = fabs(body.x) > max ? fabs(body.x) : max;
        max = fabs(body.y) > max ? fabs(body.y) : max;
        max = fabs(body.z) > max ? fabs(body.z) : max;
        // Save max into shared memory
        sdata[tid] = max;
    } else {
        // If thread is outer than values set local max to 0
        sdata[tid] = 0;
    }
    //printf("thread %d max %lf\n", myId, sdata[tid]);
    __syncthreads();

    if (myId < ents_sz) {
        // Do reduction in shared memory
        for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
        {
            if (s == 3)
                s = 4;
            if (tid < s && sdata[tid + s] > sdata[tid]) {
                sdata[tid] = sdata[tid + s];
            }
            __syncthreads();
        }
        // Write result for this block to global mem
        if (tid == 0) {
            g_odata[blockIdx.x] = sdata[0];
        }
    }
}

__global__
void get_bounding_box(double *g_idata, int ents_sz, double *g_odata)
{
    extern __shared__ double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint tid = threadIdx.x;
    int myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        sdata[tid] = g_idata[myId];
        //printf("Tread %d insert %lf\n", myId, sdata[tid]);
    } else {
        sdata[tid] = 0;
    }
    __syncthreads();

    // Do reduction in shared memory
    if (myId < ents_sz) {
        for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
        {
            // Se il blocco ha una dimensione che è multiplo di 3
            // vengono saltate alcuni elementi, rimettendo la stride
            // a 4 ci sarà una piccola sovrapposizione però si esaminano
            // tutti gli elementi
            if (s == 3)
                s = 4;
            if (tid < s && sdata[tid + s] > sdata[tid])  {
                sdata[tid] = sdata[tid + s];
            }
            __syncthreads();
        }
        // Write result for this block to global mem
        if (tid == 0) {
            g_odata[blockIdx.x] = sdata[0];
        }
    }
}

__global__ void set_branch_values(Octtree *tree)
{
    int sz_threads = gridDim.x * blockDim.x;
    // used as cache
    int children_cache[8];
    double center[3];
    int child_indx, last_node, ents;
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int cache_sz, null_counter;
    // new_mass is needed because we need old and new mass value at the same time
    double mass, new_mass;
    last_node = tree->firstfree - 1;
    // TODO check if it's better to give the first not allocated node to the first free thread
    for (int indx = last_node - id; indx >= tree->root; indx -= sz_threads)
    {
        // setting value on a new node, initialize variables
        cache_sz = 0;
        new_mass = 0;
        mass = 0;
        ents = 0;
        center[0] = 0;
        center[1] = 0;
        center[2] = 0;
        null_counter = 0;
        for (int i = 0; i < 8; i++)
        {
            // get the child index
            child_indx = tree->children[indx * 8 + i];
            if (child_indx != -1)
            {
                // move at the begin not null branches
                // TODO maybe an if on null_counter can be more efficient
                // even better save new children array and copy once
                tree->children[indx * 8 + i - null_counter] = child_indx;
                // if child is not calculated yet
                if (tree->node_pos[child_indx].w == 0)
                {
                    children_cache[cache_sz] = child_indx;
                    cache_sz++;
                }
                else
                {
                    new_mass += tree->node_pos[child_indx].w;
                    ents += 1;
                    // TODO optimized, just divide at the end
                    center[0] = (tree->node_pos[child_indx * 3].x * tree->node_pos[child_indx].w / new_mass) + (center[0] * mass / new_mass);
                    center[1] = (tree->node_pos[child_indx * 3 + 1].y * tree->node_pos[child_indx].w / new_mass) + (center[1] * mass / new_mass);
                    center[2] = (tree->node_pos[child_indx * 3 + 2].z * tree->node_pos[child_indx].w / new_mass) + (center[2] * mass / new_mass);
                    mass = new_mass;
                }
            }
            else
            {
                null_counter++;
            }
        }
        // set the null branches at the end
        for (int i = 8 - null_counter; i < 8; i++)
        {
            tree->children[indx * 8 + i] = -1;
        }
        // TODO resolve divergence
        do
        {
            int completed = 0;
            // more easy to iterate
            for (int i = 0; i < cache_sz; i++)
            {
                child_indx = children_cache[i];
                if (tree->node_pos[child_indx].w > 0)
                {
                    new_mass += tree->node_pos[child_indx].w;
                    ents += 1;
                    // TODO optimized, just divide at the end
                    center[0] = (tree->node_pos[child_indx * 3].x * tree->node_pos[child_indx].w / new_mass) + (center[0] * mass / new_mass);
                    center[1] = (tree->node_pos[child_indx * 3 + 1].y * tree->node_pos[child_indx].w / new_mass) + (center[1] * mass / new_mass);
                    center[2] = (tree->node_pos[child_indx * 3 + 2].z * tree->node_pos[child_indx].w / new_mass) + (center[2] * mass / new_mass);
                    mass = new_mass;
                    completed++;
                }
                else
                {
                    // compact not completed children at the beginning of the cache
                    children_cache[i - completed] = children_cache[i];
                }
            }
            cache_sz -= completed;

            if (cache_sz == 0)
            {
                //copy_arrs3(&tree->center[indx * 3], center, 0);
                //tree->ents[indx] = ents;
                tree->node_pos[indx].x = center[0];
                tree->node_pos[indx].y = center[1];
                tree->node_pos[indx].z = center[2];

                tree->node_pos[indx].w = mass;
                __threadfence();
            }
        } while (cache_sz > 0);
    }
}

__global__
void set_max(double *max, Octtree *tree) {
    tree->max = *max * 2;
}

__host__
void bounding_box(double *positions, double *max_in, double *max_out, int ents_sz, int step, Octtree *tree) {

    dim3 grid(grid_s, 1, 1);
    dim3 block(block_s, 1, 1);
    int shsz = block.x; // Shared memory size

    // Prima iterazione per estrarre i max dalle xyz di ogni corpo
    get_bounding_box_fi<<<grid, block, shsz * sizeof(double)>>>(positions, ents_sz, max_in, step);
    cudaDeviceSynchronize();

    // Trova il massimo tra i massimi
    double *temp;
    int ents = ents_sz;
    int grid_size = ceil((float)ents / (float)block.x);
    grid.x = grid_size;

    // Ad ogni iterazione per trovare il massimo rimangono un numero di elementi
    // pari a quello del numero dei blocchi.
    while(ents > 1) {
        get_bounding_box<<<grid_size, block, shsz * sizeof(double)>>>(max_in, ents, max_out);
        cudaDeviceSynchronize();
        // Prendo le misure per la prossima iterazione
        ents = grid_size;
        grid_size = ceil((float) grid.x / (float) block.x);

        grid.x = grid_size;

        temp = max_in;
        max_in = max_out;
        max_out = temp;
    }
    set_max<<<1, 1>>>(max_in, tree);
    cudaDeviceSynchronize();
}

__device__
int get_indx_loc(double4 *pos, double3 *center, double *border_size) {
    int indx;
    int x, y, z;
    // used to calculate the new center, it's the new border divided by 2,
    // equals
    double bord_div4 = *border_size / 4;
    z = pos->z >= center->z;
    y = pos->y >= center->y;
    x = pos->x >= center->x;
    indx = z * 4 + y * 2 + x;
    // used to calculate new center
    center->x += x ? bord_div4 : -(bord_div4);
    center->y += y ? bord_div4 : -(bord_div4);
    center->z += z ? bord_div4 : -(bord_div4);
    *border_size /= 2;

    return indx;
}

__global__
void add_ent(Octtree *tree, double4 *ents, int ents_sz) {
    // allocated is used as a boolean
    int allocated, node_indx, body_pos;
    //Octnode *node;
    double border_size;
    double3 volume_center, prev_vol;


    int myId = threadIdx.x + blockIdx.x * blockDim.x;
    if (myId < ents_sz) {
        double4 ent = ents[myId];
        //printf("Thread %d body: x: %lf, y: %lf, z: %lf, mass: %lf\n", myId, ent.x, ent.y, ent.z, ent.w);
        // set init value
        allocated = 0;
        // keep last visited node index
        node_indx = tree->root;
        //node = &tree->nodes[node_indx];
        border_size = tree->max;

        // set center of whole volume
        volume_center = {0, 0, 0};

        while (!allocated) {
            // center and border_size are updated to the next branch value
            prev_vol = volume_center;
            body_pos = get_indx_loc(&ent, &volume_center, &border_size);
            int i = node_indx * 8 + body_pos;

            // Se si usa c&s nel modo classico con in while va in deadlock
            // Con l'esecuzione parallela delle istruzioni dei thread in un warp
            // nessuno riesce a passare oltre il while.
            if (atomicCAS(&tree->mutex[i], 0, 1) == 0) {
                __threadfence();
                //set_lock(tree, i);
                //omp_set_lock(&node->writelocks[body_pos]);
                //printf("Thread %d gain lock for node: %d child: %d, border size: %lf\n", myId, node_indx, body_pos, border_size);

                if (tree->children[i] == -1) {
                    atomicExch(&tree->children[i], myId); // Scrivo sul nodo chi è il figlio
                    atomicAdd(&tree->ents[node_indx], 1); // Aggiorno la quantità di entità

                    tree->node_pos[myId] = ent; // Copio i dati del corpo sul nodo libero
                    //tree->node_pos[myId].x = -myId;
                    //tree->node_pos[myId].y = -myId;
                    //tree->node_pos[myId].z = body_pos;

                    tree->ents[myId] = 1; // Aggiorno il numero di entità sulla foglia
                    tree->parent[myId] = node_indx; // Sistemo il padre della foglia
                    //printf("Thread %d body POSIZIONATO SU LIBERO %d child %d: x: %lf, y: %lf, z: %lf, mass: %lf\n", myId, node_indx, body_pos, tree->node_pos[myId].x, tree->node_pos[myId].y, tree->node_pos[myId].z, tree->node_pos[myId].w);

                    /*
                    tree->nodes[myId].center = ent->pos;
                    tree->nodes[myId].mass = ent->mass;
                    tree->nodes[myId].ents = 1;
                    tree->nodes[myId].parent = node_indx;
                    */

                    allocated = 1;
                    //omp_unset_lock(&node->writelocks[body_pos]);
                    atomicExch(&tree->mutex[i], 0);
                    __threadfence();
                } else {
                    // if the location is occupied by a body-leaf
                    // [leaf, leaf, leaf, ..., root, branch, branch, ...]
                    if (tree->children[i] < tree->root) {
                        //printf("Tread %d trovato occupato da %d\n", myId, tree->children[i]);
                        // other is the other leaf
                        double3 other_center = volume_center;
                        double other_border = border_size;
                        //int other = node->children[body_pos];
                        int other_id = tree->children[i];
                        int other_indx = body_pos;
                        int prev_node_indx;

                        // Add new body to count in the current node
                        atomicAdd(&tree->ents[node_indx], 1);
                        ////printf("Thread %d atomic add node %d ents %d\n", myId, node_indx, tree->ents[node_indx]);

                        // When the leaves will be in different position exit the loop
                        while (body_pos == other_indx) {

                            // take first free location and set the parent of the new
                            // branch
                            //printf("Thread %d body pos: %d, other body pos: %d\n", myId, body_pos, other_indx);
                            int free;

                            free = init_node(tree);
                            //printf("Thread %d nodo %d inizializzato - parent %d\n", myId, free, node_indx);
                            // Imposto il padre della nuova foglia
                            tree->parent[free] = node_indx;

                            // set the new node as child
                            atomicExch(&tree->children[node_indx * 8 + body_pos], free);
                            atomicExch(&tree->ents[free], 2);

                            // get leaves position in the new branch
                            int old_body_pos = body_pos;
                            body_pos = get_indx_loc(&ent, &volume_center, &border_size);

                            // the center of the leaf is the position of the entity
                            // associated
                            other_indx = get_indx_loc(&tree->node_pos[other_id], &other_center, &other_border);
                            // use the new branch as the current one
                            prev_node_indx = node_indx;
                            node_indx = free;
                            //omp_set_lock(&tree->nodes[node_indx].writelocks[body_pos]);
                            //set_lock(tree, (node_indx * 8 + body_pos));

                            // Il nodo è nuovo, nessuno può averlo lockato finche non rilascio
                            // il lock al nodo padre
                            atomicCAS(&tree->mutex[node_indx*8+body_pos], 0, 1);
                            if (other_indx != body_pos) {
                                //omp_set_lock( &tree->nodes[node_indx].writelocks[other_indx]);
                                //set_lock(tree, (node_indx * 8 + other_indx));
                                atomicCAS(&tree->mutex[node_indx*8+other_indx], 0, 1);
                            }
                            //omp_unset_lock(&node->writelocks[old_body_pos]);
                            //unset_lock(tree, (prev_node_indx * 8 + old_body_pos));
                            atomicExch(&tree->mutex[prev_node_indx*8+old_body_pos], 0);
                            __threadfence();
                            //node = &tree->nodes[node_indx];
                        }

                        // set new parent in the leaves values
                        tree->parent[other_id] = node_indx;
                        tree->parent[myId] = node_indx;

                        tree->ents[myId] = 1;
                        tree->node_pos[myId] = ent;
                        //tree->node_pos[myId].x = myId;
                        //tree->node_pos[myId].y = myId;
                        //tree->node_pos[myId].z = body_pos;

                        // set the leaves as branch children
                        atomicExch(&tree->children[node_indx * 8 + body_pos], myId);
                        atomicExch(&tree->children[node_indx * 8 + other_indx], other_id);

                        //omp_unset_lock(&node->writelocks[body_pos]);
                        //omp_unset_lock(&node->writelocks[other_indx]);
                        //unset_lock(tree, (node_indx * 8 + body_pos));
                        //unset_lock(tree, (node_indx * 8 + other_indx));


                        ////printf("Thread %d body POSIZIONATO SU FONDO: x: %lf, y: %lf, z: %lf, mass: %lf\n", myId, tree->node_pos[myId].x, tree->node_pos[myId].y, tree->node_pos[myId].z, tree->node_pos[myId].w);
                        //printf("Thread %d riposizionati body: my=%d - other=%d\n", myId, myId, other_id);
                        atomicExch(&tree->mutex[node_indx*8+body_pos], 0);
                        atomicExch(&tree->mutex[node_indx*8+other_indx], 0);
                        __threadfence();

                        allocated = 1;
                    } else { // Descend into the tree
                        // The current node will have one more body in its subtree
                        atomicAdd(&tree->ents[node_indx], 1);
                        //omp_unset_lock(&node->writelocks[body_pos]);
                        //unset_lock(tree, i);
                        atomicExch(&tree->mutex[i], 0);
                        __threadfence();
                        // cross the branch
                        //printf("Thread %d lascia %d attraversa verso %d\n", myId, i, body_pos);
                        node_indx = tree->children[node_indx * 8 + body_pos];
                        //node = &tree->nodes[node_indx];
                    }
                }
            } else {
                //printf("Thread %d lock is closed node: %d child: %d, borders size: %lf\n", myId, node_indx, body_pos, border_size);
                // Se trovo lockato ripristino i valori a ingresso del loop
                border_size *= 2;
                volume_center = prev_vol;
            }
        }
    }
}

__global__
void center_of_mass(Octtree *tree, int ents_sz, int block) {
    int myId, myNode;
    //int sz_threads = gridDim.x * blockDim.x;

    // The calculation of the centers of mass is done starting from
    // the bottom and reaching up to the root of the tree.
    myId = blockIdx.x * blockDim.x + threadIdx.x;
    myNode = (tree->firstfree - 1) - myId;
    int done = 0;
    //for (myNode = (tree->firstfree -1) - myId; myNode >= tree->root; myNode -= sz_threads){}

    if (myNode >= ents_sz){
        double new_mass, l_mass;
        double4 child;
        double3 l_center;
        int j;
        int counter = 0;

        //int last_child = -1;

        l_center.x = tree->node_pos[myNode].x;
        l_center.y = tree->node_pos[myNode].y;
        l_center.z = tree->node_pos[myNode].z;
        l_mass = tree->node_pos[myNode].w;

        j = tree->children[myNode * 8 + counter];
        while (!done) {
            //Octnode *my_node = &tree->nodes[n];
            if (counter < 8) {
                if (j == -1) {
                    counter++;
                    j = tree->children[myNode * 8 + counter];
                } else if (tree->node_pos[j].w != 0) {
                    child = tree->node_pos[j];

                    new_mass = l_mass + child.w;

                    l_center.x = (child.x * child.w / new_mass) + (l_center.x * l_mass / new_mass);
                    l_center.y = (child.y * child.w / new_mass) + (l_center.y * l_mass / new_mass);
                    l_center.z = (child.z * child.w / new_mass) + (l_center.z * l_mass / new_mass);
                    l_mass = new_mass;

                    printf("Thread %d node %d calculated child %d node %d\n", myId, myNode, counter, j);

                    counter++;
                    j = tree->children[myNode * 8 + counter];
                }
            } // if counter < 8

            if (counter >= 8) {
                tree->node_pos[myNode].x = l_center.x;
                tree->node_pos[myNode].y = l_center.y;
                tree->node_pos[myNode].z = l_center.z;
                tree->node_pos[myNode].w = l_mass;

                __threadfence();
                printf("Thread %d node %d done counter %d\n", myId, myNode, counter);
                done = 1;
            }
        } // While done
    } // if mynode>ents_sz
}

__global__
void sort_tree_leaves(Octtree *tree, int *sorted_nodes, int ents_sz) {
    int myId;
    int offset;
    int child_i;
    int not_done;
    int my_node_i;
    int i;
    //int child_ents;

    myId = blockIdx.x * blockDim.x + threadIdx.x;
    offset = 0;

    my_node_i = myId + ents_sz;

    if (my_node_i < tree->firstfree){
        if (my_node_i == ents_sz){ // Sono la root
            tree->mutex[my_node_i] = 1;
            //printf("Mutex root: %d\n", tree->mutex[ents_sz]);
        }

        not_done = 1;
        //my_parent = (my_node_i == ents_sz) ? 0 : tree->parent[myId + ents_sz];
        //if (myId == 0) printf("My parent: %d - ents_sz: %d\n", my_parent, ents_sz);
        //my_node_i = myId + ents_sz;
        while (not_done) {
            if (tree->mutex[my_node_i] == 1) {
                //printf("Thread %d pre check positions node: %d, parent: %d\n", myId, my_node_i, my_parent);
                offset = (my_node_i == ents_sz) ? 0 : tree->mutex[my_node_i+tree->firstfree];
                i = -1;
                for (i = 0; i < 8; i++) {
                    child_i = tree->children[my_node_i * 8 + i];
                    if (child_i != -1) {
                        sorted_nodes[offset] = child_i;

                        tree->mutex[child_i+tree->firstfree] = offset;
                        tree->mutex[child_i] = 1;

                        offset += tree->ents[child_i];
                        __threadfence();

                        // Restore mutex array
                        tree->mutex[my_node_i] = 0;
                        tree->mutex[my_node_i+tree->firstfree] = 0;
                        //printf("Thread %d unlock node: %d\n", myId, child_i);
                    }
                    if (i == 7) {
                        //printf("Thread %d finish write positions\n", myId);
                        not_done = 0;
                    }
                } // for over children
            } // if parent == 1
        } // while not_done
    }
}

__global__
void print_sorted(int *sorted_nodes, int ents_sz) {
    for (int i = 0; i < ents_sz; i++) {
        printf("%d ", sorted_nodes[i]);
    }
    printf("\n");
}

__global__
void ci_sono_tutti_i_numeri(int *sorted_nodes, int ents_sz) {
    for (int i = 0; i < ents_sz; i++) {
        for (int j = 0; j < ents_sz; j++) {
            if (i == sorted_nodes[j]){
                break;
            }
        }
        printf("Missing: %d\n", i); // qua dovrebbe arrivare solo se completa il giro
    }
}

__global__
void find_dups(int *sorted_nodes, int ents_sz) {
    for (int i = 0; i < ents_sz; i++) {
        for (int j = 0; j < ents_sz; j++) {
            if (sorted_nodes[i] == sorted_nodes[j] && i != j){
                printf("Duplicate of %d at pos %d found at %d\n", sorted_nodes[i], i, sorted_nodes[j]);
            }
        }
    }

}

__global__
void mini_print(Octtree *tree, int ents_sz) {
    printf("Mini print\n");
    double4 *node_pos = tree->node_pos;
    for (int id = 0; id < ents_sz; id++) {
        printf("id: %d, (x:%lf, y:%lf, z:%lf), mass:%lf, parent: %d\n",
           id, node_pos[id].x, node_pos[id].y, node_pos[id].z, node_pos[id].w, tree->parent[id]);
    }
}


int main(int argc, char *argv[]) {
    if (argc < 8) {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename grid_size block_size\n", argv[0]);
        exit(1);
    }

    grid_s = atoi(argv[6]);
    block_s = atoi(argv[7]);

    uint n_ents;
    float start, end, dt;
    int n_steps;
    struct timespec s, e;

    double *positions;
    double *velocities;

    get_entities(argv[1], &n_ents, &positions, &velocities);

    start = strtof(argv[2], NULL);
    end = strtof(argv[3], NULL);
    dt = strtof(argv[4], NULL);

    n_steps = (end - start) / dt;

    printf("Start: %f, end: %f, delta time: %f, time steps: %d, ents: %d\n", start, end, dt, n_steps, n_ents);

    clock_gettime(CLOCK_REALTIME, &s);
    run(positions, velocities, n_ents, n_steps, dt, argv[5]);
    clock_gettime(CLOCK_REALTIME, &e);

    printf("Completed. Output file: %s\n", argv[5]);

    double time_spent =
        (e.tv_sec - s.tv_sec) + (e.tv_nsec - s.tv_nsec) / BILLION;

    printf("Elapsed wall time: %f s\n", time_spent);


    free(positions);
    free(velocities);
}

__host__
void calculate_block_dim() {
    //
}

__host__
void do_all_mallocs_and_copy(int ents_sz, int n_steps, double **d_positions, double **d_velocities, double **d_accelerations, double **d_loc_max_in, double **d_loc_max_out, double *h_positions, double *h_velocities, int **d_sorted_nodes) {

    cudaError_t error;

    // Malloc
    error = cudaMalloc((void **)d_positions, ents_sz * 4 * sizeof(double) * n_steps + 1);
    cuda_check_error(error, "Devce positions malloc\n");

    error = cudaMalloc((void **)d_velocities, ents_sz * 3 * sizeof(double));
    cuda_check_error(error, "Devce velocities malloc\n");

    error = cudaMalloc((void **)d_accelerations, ents_sz * 3 * sizeof(double));
    cuda_check_error(error, "Devce accelerations malloc\n");

    error = cudaMalloc((void **)d_loc_max_in, ents_sz* sizeof(double));
    cuda_check_error(error, "Devce local max malloc\n");

    error = cudaMalloc((void **)d_loc_max_out, ents_sz * sizeof(double));
    cuda_check_error(error, "Device local max malloc\n");

    error = cudaMalloc((void **)d_sorted_nodes, ents_sz * sizeof(int));
    cuda_check_error(error, "Device sorted nodes malloc\n");


    // Copy
    error = cudaMemcpy(*d_positions, h_positions, ents_sz * 4 * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error, "Host do Device copy positions\n");

    error = cudaMemcpy(*d_velocities, h_velocities, ents_sz * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error, "Host to Device copy velocities\n");

}

__host__
void sort_tree(Octtree *tree, int *sorted_nodes, int *mutex, int ents_sz, int first_free) {
    int size = ceil((double)first_free/(double)block_s);
    dim3 grid(size, 1, 1);
    dim3 block(block_s, 1, 1);

    printf("Sorting tree grid: %d, block size: %d, first free: %d\n", size, block_s, first_free);
    sort_tree_leaves<<<grid, block>>>(tree, sorted_nodes, ents_sz);
    cudaDeviceSynchronize();
    printf("Sorting finished\n");

    reset_mutex(mutex, ents_sz);
}



__host__
void run(double *h_positions, double *h_velocities, uint ents_sz, int n_steps, double dt, const char *output) {
    FILE *fpt;
    cudaError_t error;
    double *d_positions;
    double *d_velocities;
    double *d_accelerations;
    double *d_loc_max_in, *d_loc_max_out;

    // Tree data
    Octtree *d_tree, *h_tree;
    double4 *d_node_pos;
    int *d_ents, *d_parent, *d_children, *d_mutex, *d_sorted_nodes;

    dim3 grid(grid_s, 1, 1);
    dim3 block(block_s, 1, 1);
    int shsz = block.x; // Shared memory size

    if (block_s * grid_s < ents_sz) {
        fprintf(stderr, "Insufficent threads!\nTotal cuda threads requested: %d - Total bodies: %d\n", block_s * grid_s, ents_sz);
        exit(1);
    }

    do_all_mallocs_and_copy(ents_sz, n_steps, &d_positions, &d_velocities, &d_accelerations, &d_loc_max_in, &d_loc_max_out, h_positions, h_velocities, &d_sorted_nodes);

    init_tree(&d_tree, ents_sz, &d_node_pos, &d_ents, &d_parent, &d_children, &d_mutex);

    bounding_box(d_positions, d_loc_max_in, d_loc_max_out, ents_sz, 0, d_tree);

    add_ent<<<grid, block>>>(d_tree, (double4 *)d_positions, ents_sz);
    cudaDeviceSynchronize();

    //print_tree(d_tree, ents_sz, d_node_pos, d_children, d_ents);
    //cudaDeviceSynchronize();


    h_tree = (Octtree *)malloc(sizeof(Octtree));

    error = cudaMemcpy(h_tree, d_tree, sizeof(Octtree), cudaMemcpyDeviceToHost);
    cuda_check_error(error, "Tree to Device copy velocities\n");

    int g = ((double)(h_tree->firstfree) / (double)block_s) + 1;
    dim3 grid2(g, 1, 1);

    printf("Center of mass\n");
    printf("Centering tree grid: %d, block size: %d, first free: %d\n", g, block_s, h_tree->firstfree);
    center_of_mass<<<g, block>>>(d_tree, ents_sz, block_s);
    //set_branch_values<<<grid, block>>>(d_tree);
    cudaDeviceSynchronize();

    //print_tree(d_tree, ents_sz, d_node_pos, d_children, d_ents);
    //cudaDeviceSynchronize();

    //print_mutex<<<1,1>>>(d_tree, ents_sz);
    //cudaDeviceSynchronize();

    sort_tree(d_tree, d_sorted_nodes, d_mutex, ents_sz, h_tree->firstfree);
    printf("Sort finish");

    print_sorted<<<1, 1>>>(d_sorted_nodes, ents_sz);
    cudaDeviceSynchronize();


    ci_sono_tutti_i_numeri<<<1, 1>>>(d_sorted_nodes, ents_sz);
    cudaDeviceSynchronize();


    find_dups<<<1, 1>>>(d_sorted_nodes, ents_sz);
    cudaDeviceSynchronize();
    // TODO: calcolo accelerazione
    exit(0);

    acceleration<<<grid, block, shsz * sizeof(double) * 4>>>(d_positions, d_accelerations, ents_sz, 0);

    for (int t = 1; t <= n_steps; t++) {

        update_velocities<<<grid, block>>>(d_accelerations, d_velocities, ents_sz, dt);
        cudaDeviceSynchronize();

        update_positions<<<grid, block>>>(d_positions, d_velocities, ents_sz, dt, t);
        cudaDeviceSynchronize();

        acceleration<<<grid, block, shsz * sizeof(double) * 4>>>(d_positions, d_accelerations, ents_sz, t);
        cudaDeviceSynchronize();

        update_velocities<<<grid, block>>>(d_accelerations, d_velocities, ents_sz, dt);
        cudaDeviceSynchronize();
    }

    double *h_pos;
    safe_malloc_double(&h_pos, ents_sz * 4 * n_steps);
    error = cudaMemcpy(h_pos, d_positions, (size_t)(ents_sz * 4 * sizeof(double) * n_steps), cudaMemcpyDeviceToHost);
    cuda_check_error(error, "Device to Host copy positions\n");
    double4 *h_poss = (double4 *) h_pos;

    fpt = fopen(output, "w");
    for (int i = 0; i < ents_sz * n_steps; i++)
        fprintf(fpt, "%d,%lf,%lf,%lf,%lf\n", i%ents_sz, h_poss[i].x, h_poss[i].y, h_poss[i].z, h_poss[i].w);


    free(h_pos);
    cudaFree(d_positions);
    cudaFree(d_velocities);
    cudaFree(d_accelerations);
    cudaFree(d_loc_max_in);
    cudaFree(d_loc_max_out);
    cudaFree(d_sorted_nodes);

    fclose(fpt);
}

/**
 * Estimate the number of bodies by counting the lines of the file
 */
__host__
void count_entities_file(char *filename, uint *n){
    FILE *file;
    char c;

    file = fopen(filename, "r");
    if (file == NULL){
        fprintf(stderr, "Error opening file '%s'\n", filename);
        exit(1);
    }

    *n = 0;
    for (c = getc(file); c != EOF; c = getc(file))
        if (c == '\n')
            (*n)++;
    fclose(file);
    if (n == 0){
        fprintf(stderr, "No bodies found into file. Closing\n");
        exit(1);
    } else {
        // For prevent files that do not have a newline character at the end
        (*n)++;
    }
}

/**
 * Allocate array for doubles and check for errors
 */
__host__
void safe_malloc_double(double **pointer, int quantity) {
    *pointer = (double *) malloc(quantity * sizeof(double));
    if (*pointer == NULL){
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }
}

/**
 * Read file and generate an arrays for masses, positions and velocities for the bodies
 */
__host__
void get_entities(char filename[], uint *n_ents, double **positions, double **velocities) {
    FILE *file;
    int i;

    count_entities_file(filename, n_ents);

    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error while opening file '%s'\n", filename);
        exit(1);
    }

    safe_malloc_double(positions, *n_ents * 4);
    safe_malloc_double(velocities, *n_ents * 3);
    i = 0;

    while ((fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
                        &((*positions)[i*4  ] ),
                        &((*positions)[i*4+1] ),
                        &((*positions)[i*4+2] ),
                        &((*velocities)[i*3  ]),
                        &((*velocities)[i*3+1]),
                        &((*velocities)[i*3+2]),
                        &((*positions)[i*4+3]))) == 7) {
        i++;
    }
    // check if while ended because the end of the file has been reached
    if (fgetc(file) != EOF) {
        fprintf(stderr, "Error while reading file '%s': file is not well formed\n", filename);
        fclose(file);
        exit(1);
    }
    // Update n_ents with the correct number of scanned lines
    *n_ents = i;
    fclose(file);
}

__global__
void acceleration(double *positions, double *acc, uint ents_sz, int step) {
    extern __shared__ double4 shPositions[];
    double4 *gPositions;
    double3 *gAcc;
    double3 r_vector, l_acc;
    double inv_r3;
    int myId;
    int i;
    int iteration;
    int offset;
    double3 myBodyPos;

    positions = positions + (ents_sz * 4 * step);
    gPositions = (double4 *) positions;
    //gPositions += ents_sz * step;
    gAcc = (double3 *) acc;
    l_acc = {0.0, 0.0, 0.0};
    myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        // Mi tengo le mie posizioni per non accedere alla memoria globale
        myBodyPos = {gPositions[myId].x, gPositions[myId].y, gPositions[myId].z};
    }

    for (i = 0, iteration = 0; i < ents_sz; i += blockDim.x, iteration++) {
        offset = blockDim.x * iteration; // Calcola offset
        if (offset + threadIdx.x < ents_sz){
            // Ogni thread carica un corpo nella shared memory
            shPositions[threadIdx.x] = gPositions[offset + threadIdx.x];
        }
        __syncthreads(); // Mi assicuro che tutti i thread carichino un corpo
        for (int j = 0; j < blockDim.x && offset+j < ents_sz; j++){
            if (threadIdx.x != j) {

                r_vector.x = shPositions[j].x - myBodyPos.x;
                r_vector.y = shPositions[j].y - myBodyPos.y;
                r_vector.z = shPositions[j].z - myBodyPos.z;

                inv_r3 = r_vector.x * r_vector.x + r_vector.y * r_vector.y +
                            r_vector.z * r_vector.z + 0.01;
                inv_r3 = pow(inv_r3, -1.5);

                l_acc.x += BIG_G * r_vector.x * inv_r3 * shPositions[j].w;
                l_acc.y += BIG_G * r_vector.y * inv_r3 * shPositions[j].w;
                l_acc.z += BIG_G * r_vector.z * inv_r3 * shPositions[j].w;
            }
        }
    }
    __syncthreads();

    if (myId < ents_sz) {
        gAcc[myId].x = l_acc.x;
        gAcc[myId].y = l_acc.y;
        gAcc[myId].z = l_acc.z;
    }
}

__global__
void update_velocities(double *acc, double *vel, uint ents_sz, double dt) {
    int myId;

    myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        vel[myId*3    ] += acc[myId*3    ] * dt / 2.0;
        vel[myId*3 + 1] += acc[myId*3 + 1] * dt / 2.0;
        vel[myId*3 + 2] += acc[myId*3 + 2] * dt / 2.0;
    }
}

__global__
void update_positions(double *pos, double *velocities, uint ents_sz, double dt, int step) {
    double4 *gPos;
    double3 *gVelocities;
    int myId;

    gPos = (double4 *) pos;
    gVelocities = (double3 *) velocities;

    myId = blockIdx.x * blockDim.x + threadIdx.x;
    int new_pos = step * ents_sz + myId;
    int old_pos = (step - 1) * ents_sz + myId;

    if (myId < ents_sz) {
        gPos[new_pos].x = gPos[old_pos].x + gVelocities[myId].x * dt;
        gPos[new_pos].y = gPos[old_pos].y + gVelocities[myId].y * dt;
        gPos[new_pos].z = gPos[old_pos].z + gVelocities[myId].z * dt;
        gPos[new_pos].w = gPos[old_pos].w;
    }
}

