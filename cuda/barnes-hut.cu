#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#define LOCKED -2
#define PRINT_LOOP

//compute-sanitizer --tool memcheck
// https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html
typedef unsigned int uint;

typedef struct
{
    double x;
    double y;
    double z;
} RVec3;

typedef struct
{
    double *pos; // allocation size must be three times the size of mass
    double *vel; // allocation size must be three times the size of mass
    double *mass;
} Entities;

typedef struct
{
    int firstfree; // first location free
    int root;
    double max;
    // Octnodes
    uint *ents;
    double *mass;
    double *center; // allocation size must be three times the size of mass
    int *parent;
    int *children; // must have 8 slots for each node

} Octtree;

typedef struct {
    int id;
    double3 center;
    double mass;
    double border;
} Stacknode;

#define RESULTS

// const double BIG_G = 6.67e-11;
const double BIG_G = 1;
const double THETA = 0.5; // 1;
const uint FULL_MASK=0xFFFFFFFF;
// TODO check from version

// a lot of optimization like use shift instead of multiplication will be made by compiler

uint get_entities(char filename[], Entities *ents)
{
    // Entity e_buff;
    double pos_buff[3];
    double vel_buff[3];
    double mass_buff;
    uint ret_size;
    uint size;
    double *pos_ret;
    double *vel_ret;
    double *mass_ret;

    FILE *f = fopen(filename, "r");

    // Check if file has been open correctly, if not return NULL
    if (!f)
    {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }

    
    pos_ret = (double *)malloc(3 * sizeof(double));
    vel_ret = (double *)malloc(3 * sizeof(RVec3));
    mass_ret = (double *)malloc(1 * sizeof(double));

    size = 0;
    ret_size = 1;
    // fscanf return the number of input items successfully matched and assigned
    while (
        (fscanf(f, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &pos_buff[0],
                         &pos_buff[1], &pos_buff[2], &vel_buff[0], &vel_buff[1],
                         &vel_buff[2], &mass_buff)) == 7)
    {
        size++;
        if (ret_size < size)
        {
            
            ret_size *= 2;
            // ret = (Entity *)realloc((void *)ret, ret_size * sizeof(Entity));
            pos_ret = (double *)realloc(pos_ret, ret_size * 3 * sizeof(double));
            vel_ret = (double *)realloc(vel_ret, ret_size * 3 * sizeof(double));
            mass_ret = (double *)realloc(mass_ret, ret_size * sizeof(double));
        }
        // Save value in first free location
        // ret[size - 1] = e_buff;
        pos_ret[size * 3 - 3] = pos_buff[0];
        pos_ret[size * 3 - 2] = pos_buff[1];
        pos_ret[size * 3 - 1] = pos_buff[2];
        vel_ret[size * 3 - 3] = vel_buff[0];
        vel_ret[size * 3 - 2] = vel_buff[1];
        vel_ret[size * 3 - 1] = vel_buff[2];
        mass_ret[size - 1] = mass_buff;
    }

    // check if while ended because the end of the file has been reached
    if (fgetc(f) != EOF)
    {
        fprintf(stderr, "Error reading file '%s': file is not well formed\n",
                filename);
        fclose(f);
        return 0;
    }

    // *ents = ret;
    ents->pos = pos_ret;
    ents->vel = vel_ret;
    ents->mass = mass_ret;
    fclose(f);
    return size;
}

__device__ double border_tree(Octtree *tree) { return tree->max * 2; }

// copy array of 3 elements, usually used for position array
__device__ void copy_arrs3(double *dest, double *src, uint indx)
{
    dest[indx] = src[indx];
    dest[indx + 1] = src[indx + 1];
    dest[indx + 2] = src[indx + 2];
}

__device__ void init_node(Octtree *tree, int indx, int depth)
{
    tree->center[indx * 3] = 0;
    tree->center[indx * 3 + 1] = 0;
    tree->center[indx * 3 + 2] = 0;
    tree->mass[indx] = 0;
    tree->ents[indx] = 0;
    tree->parent[indx] = -1;
    // tree->depth[indx] = depth;
    for (uint i = 0; i < 8; i++)
    {
        tree->children[indx * 8 + i] = -1;
    }
}

// pos and center must point to the first axis of position
__device__ int get_indx_loc(double *pos, double *center, double *border_size)
{
    int indx;
    int x, y, z;
    // used to calculate the new center, it's the new border divided by 2,
    // equals
    double bord_div4 = *border_size / 4;
    z = pos[2] >= center[2];
    y = pos[1] >= center[1];
    x = pos[0] >= center[0];
    indx = z * 4 + y * 2 + x;
    // used to calculate new center
    center[0] +=
        x ? bord_div4 : -(bord_div4); // double(x)*2*border_size - border_size
    center[1] += y ? bord_div4 : -(bord_div4);
    center[2] += z ? bord_div4 : -(bord_div4);
    *border_size /= 2;

    return indx;
}

// KERNEL 2
//  add a entity in the tree
//  it's create all the needed branch
//  the leaf of the tree are biunivocaly (e' inglese?) associated to an entity
__global__ void add_ent(Octtree *tree, Entities *ent, uint ents_sz){
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    if(id< ents_sz){
        // printf("body %d will be added in the tree\n", id);
        int allocated, node_indx, body_pos, child_val, child_indx, root;
        double border_size;
        double volume_center[3]={0,0,0};
        double pos_ent[3];
        allocated = 0;
        root = tree->root;
        node_indx = root;
        border_size=tree->max;
        copy_arrs3(pos_ent, &ent->pos[id * 3], 0);
        copy_arrs3(&tree->center[id*3], pos_ent, 0);
        tree->mass[id] = ent->mass[id];
        tree->ents[id] = 1;

        while(!allocated){
            body_pos= get_indx_loc(pos_ent, volume_center, &border_size);
            child_indx = node_indx * 8 + body_pos;
            do{
                child_val = tree->children[child_indx];
                if(child_val!=LOCKED){
                    if(child_val == atomicCAS(&tree->children[child_indx], child_val, LOCKED)){
                        atomicAdd(&tree->ents[node_indx], 1);
                        if(child_val == -1){
                            // copy_arrs3(&tree->center[id*3], pos_ent, 0);
                            tree->children[child_indx] = id;
                            tree->parent[id] = node_indx;
                            allocated=1;
                            __threadfence();
                        }else{
                            if(child_val < root){
                                int other = child_val;
                                double other_center[3];
                                copy_arrs3(other_center, volume_center, 0);
                                double other_pos_vol[3];
                                copy_arrs3(other_pos_vol, &tree->center[other*3],0);
                                double other_border=border_size;
                                int other_pos = body_pos;

                                while(body_pos == other_pos){
                                    int free = atomicAdd(&tree->firstfree, 1);
                                    init_node(tree, free, 0);
                                    tree->parent[free]= node_indx;
                                    tree->ents[free]=2;
                                    
                                    body_pos=get_indx_loc(pos_ent, volume_center, &border_size);
                                    other_pos=get_indx_loc(other_pos_vol, other_center, &other_border);
                                    node_indx=free;
                                    tree->children[node_indx*8 + body_pos]=LOCKED;
                                    if(body_pos!=other_pos){
                                        tree->children[node_indx*8 + other_pos]=LOCKED;
                                    }
                                    tree->children[child_indx]=node_indx;
                                    child_indx=node_indx * 8 + body_pos;
                                    __threadfence();
                                }

                                tree->parent[other]=node_indx;
                                tree->parent[id]=node_indx;
                                tree->children[child_indx]=id;
                                tree->children[node_indx*8+other_pos]=other;
                                __threadfence();
                                allocated = 1;

                            }else{
                                node_indx=child_val;
                                tree->children[child_indx] = child_val;
                                __threadfence();
                            }
                        }
                    }else{
                        child_val=LOCKED;
                    }
                }
            }while(child_val==LOCKED);
        }
        // printf("body %d added in the tree\n", id);

    }
}

// KERNEL 1
__global__ void get_bounding_box(double *g_idata, int ents_sz, double *g_odata)
{
    extern __shared__ double sdata[];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint tid = threadIdx.x;
    uint space = (blockDim.x);
    uint i = (blockIdx.x * (blockDim.x * 2) + threadIdx.x);
    double v1, v2;
    // if thread are more than values give them a 0 value
    // the first if has no divergence, becuase threads have same ents_sz value
    if (i < ents_sz)
    {
        v1 = fabs(g_idata[i])*2;
        if (i + space < ents_sz)
        {
            v2 = fabs(g_idata[i + space])*2;
            sdata[tid] = v1 > v2 ? v1 : v2;
        }
        else
        {
            sdata[tid] = v1;
        }
    }
    else
    {
        sdata[tid] = 0;
    }
    __syncthreads();
    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s && sdata[tid + s] > sdata[tid])
        {
            sdata[tid] = sdata[tid + s];
        }
        __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0)
    {

        g_odata[blockIdx.x] = sdata[0];
    }
}

// KERNEL 3
//  if mass==0 node is not ready
__global__ void center_of_mass(Octtree *tree)
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
                // even better save new children array and copy once
                tree->children[indx * 8 + i - null_counter] = child_indx;
                // if child is not calculated yet
                if (tree->mass[child_indx] == 0)
                {
                    children_cache[cache_sz] = child_indx;
                    cache_sz++;
                }
                else
                {
                    new_mass += tree->mass[child_indx];
                    // ents += 1;
                    ents += tree->ents[child_indx];
                    center[0] = (tree->center[child_indx * 3] * tree->mass[child_indx] / new_mass) + (center[0] * mass / new_mass);
                    center[1] = (tree->center[child_indx * 3 + 1] * tree->mass[child_indx] / new_mass) + (center[1] * mass / new_mass);
                    center[2] = (tree->center[child_indx * 3 + 2] * tree->mass[child_indx] / new_mass) + (center[2] * mass / new_mass);
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
        
        do
        {
            int completed = 0;
            // more easy to iterate
            for (int i = 0; i < cache_sz; i++)
            {
                child_indx = children_cache[i];
                if (tree->mass[child_indx] > 0)
                {
                    new_mass += tree->mass[child_indx];
                    // ents += 1;
                    ents += tree->ents[child_indx];
                    center[0] = (tree->center[child_indx * 3] * tree->mass[child_indx] / new_mass) + (center[0] * mass / new_mass);
                    center[1] = (tree->center[child_indx * 3 + 1] * tree->mass[child_indx] / new_mass) + (center[1] * mass / new_mass);
                    center[2] = (tree->center[child_indx * 3 + 2] * tree->mass[child_indx] / new_mass) + (center[2] * mass / new_mass);
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
                copy_arrs3(&tree->center[indx * 3], center, 0);
                // tree->ents[indx] = ents;
                tree->mass[indx] = mass;
                __threadfence();
            }
        } while (cache_sz > 0);
    }
}

// KERNEL 4
//  here tree->parent will be used to store the "s", the position where to write values (tree->parent is not needed anymore)
__global__ void sort_ents(Octtree *tree, int *sorted)
{
    //all __syncthreads() have been added to increase performance
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int node_id;
    if(id<tree->root){
        node_id=tree->root;
        int count_ents=0;
        while(tree->ents[node_id]>1){
            int i=0;
            int child_ents=0;
            int child;
            do
            {
                count_ents+=child_ents;
                child=tree->children[node_id*8+i];
                i++;
                child_ents = tree->ents[child];
            // } while (count_ents+child_ents<=id);
            } while (id>=count_ents+child_ents);
            __syncwarp();
            node_id=child;
        }
        // __syncthreads();
        sorted[id]=node_id;
    }
}

// il paper mette il calcolo delle accelerazioni nel kernel 5 e l'aggiornamento e velocita' nel kernel 6
__device__ void calculate_acceleration(double3 *pos_ent, double mass_ent, double3 *pos_node, double mass_node,
                                       double3 *acc)
{
    double3 r_vector;
    r_vector.x = pos_node->x - pos_ent->x;
    r_vector.y = pos_node->y - pos_ent->y;
    r_vector.z = pos_node->z - pos_ent->z;

    double inv_r3 = r_vector.x * r_vector.x + r_vector.y * r_vector.y +
                    r_vector.z * r_vector.z + 0.01;
    inv_r3 = pow(inv_r3, -1.5);

    acc->x += BIG_G * r_vector.x * inv_r3 * mass_node;
    acc->y += BIG_G * r_vector.y * inv_r3 * mass_node;
    acc->z += BIG_G * r_vector.z * inv_r3 * mass_node;
}

__device__ double get_distance(double3 *r1, double3 *r2)
{
    return sqrt((r1->x - r2->x) * (r1->x - r2->x) + (r1->y - r2->y) * (r1->y - r2->y) + (r1->z - r2->z) * (r1->z- r2->z));
}

__device__ int get_next_node(int curr_node, int parent, int children[], int last_node)
{
    if (parent == last_node)
    {
        return children[0];
    }
    for (int i = 0; i < 7; i++)
    {
        if (children[i] == last_node)
        {
            if (children[i + 1] > -1)
            {
                return children[i + 1];
            }
            else
            {
                return curr_node;
            }
            break;
        }
    }
    return curr_node;
}

__global__ void acceleration_w_stack(Octtree *tree, Entities *ents,
                                     int ents_sz, int* sorted_nodes, double *acc_buff, size_t total_mem) {
    
    extern __shared__ Stacknode stacks[];
    // change 32 with warp size
    int entId, laneId, myWarpId, stackSz, totalWarps, myId,
        lastIndx;    // laneId is the thread id in the warp
    double distance, my_mass; // actBorder is the actual border
    double3 myAcc;
    Stacknode *myStack;
    double3 my_pos;
    myId = threadIdx.x + blockIdx.x * blockDim.x;
    laneId = threadIdx.x % 32;
    myWarpId = threadIdx.x / 32;
    totalWarps = (blockDim.x - 1) / 32 + 1;
    stackSz = total_mem / sizeof(Stacknode)/ totalWarps;
    lastIndx = 0;
    // my_pos_mass = tree->node_pos[myId];
    if(myId<ents_sz){
        entId = sorted_nodes[myId];
        my_mass = ents->mass[entId];
        my_pos.x = ents->pos[entId*3];
        my_pos.y = ents->pos[entId*3+1];
        my_pos.z = ents->pos[entId*3+2];
    }else{
        my_pos = {0,0,0};
        my_mass = 1;
        entId=0;
    }
    myAcc.x = 0;
    myAcc.y = 0;
    myAcc.z = 0;
    // __shared__ int lastIndx=0;

    myStack = &stacks[stackSz * myWarpId];
    if (laneId == 0) {
        double3 center= {tree->center[tree->root*3], tree->center[tree->root*3+1], tree->center[tree->root*3+2]};
        myStack[lastIndx] = {tree->root, center, tree->mass[tree->root], tree->max};
    }
    __syncwarp();

    while (lastIndx >= 0) {
        distance = get_distance(&myStack[lastIndx].center, &my_pos);
        double border = myStack[lastIndx].border;
        int nodeId = myStack[lastIndx].id;

        // if one or more thread need to continue the visit, all the threads in
        // the warp will continue in this way we get more accurate result with
        // the same execution time
        if (__all_sync(FULL_MASK, border / distance < THETA) ||
            tree->ents[nodeId] == 1) {
            if (nodeId != entId) {
                calculate_acceleration(&my_pos, my_mass,
                                       &myStack[lastIndx].center, myStack[lastIndx].mass, &myAcc);
            }
            if (laneId == 0) {
                lastIndx--;
            }
            lastIndx = __shfl_sync(FULL_MASK, lastIndx, 0);
        } else {
            if (laneId < 8) {
                // decrement lastIndx because we remove the node from stack
                lastIndx--;
                int indx = tree->children[nodeId * 8 + laneId];
                int mask = __ballot_sync(0x000000FF, indx != -1);
                
                int pushPos = __popc(~(FULL_MASK << (laneId+1)) & mask);
                if (indx != -1) {
                    double3 center={tree->center[indx*3],tree->center[indx*3+1],tree->center[indx*3+2]};
                    myStack[lastIndx + pushPos] = {indx, center, tree->mass[indx],
                                                   border / 2};
                }
                // last thread has the number of new node inside pushPos
                if (laneId == 7) {
                    lastIndx += pushPos;
                }
            }
            __syncwarp();
            lastIndx = __shfl_sync(FULL_MASK, lastIndx, 7);
        }
        
    }
    if(myId<ents_sz){
        // printf("acc %d. x: %lf, y: %lf, z:%lf\n", myId, myAcc.x, myAcc.y, myAcc.z);
        acc_buff[entId*3]=myAcc.x;
        acc_buff[entId*3+1]=myAcc.y;
        acc_buff[entId*3+2]=myAcc.z;

    }
}

void print_values(double *pos, double *vel, int ents_sz, FILE *fpt)
{
    for (int entity_idx = 0; entity_idx < ents_sz; entity_idx++)
    {
        fprintf(fpt, "%u,%lf,%lf,%lf,%lf,%lf,%lf \n", entity_idx, pos[entity_idx * 3],
                pos[entity_idx * 3 + 1], pos[entity_idx * 3 + 2], vel[entity_idx * 3], vel[entity_idx * 3 + 1],
                vel[entity_idx * 3 + 2]);
    }
}

// tot_threads is the minimum number of threads to be reach
// if tot_threads is 0 choose the maximum number of thread
void get_opt_grid(cudaDeviceProp *prop, uint tot_threads, uint regs_sz, uint *blocks_sz, uint *threads_sz)
{
    uint temp_block, temp_thread;
    temp_thread = prop->maxThreadsPerMultiProcessor / prop->maxBlocksPerMultiProcessor; // calculate value
    if (temp_thread * regs_sz > prop->regsPerMultiprocessor / prop->maxBlocksPerMultiProcessor)
    {
        temp_thread = prop->regsPerMultiprocessor / prop->maxBlocksPerMultiProcessor / regs_sz;
    }
    if (tot_threads == 0)
    {
        temp_block = prop->maxBlocksPerMultiProcessor * prop->multiProcessorCount;
    }
    else
    {
        temp_block = (tot_threads - 1) / temp_thread + 1;
    }

    *blocks_sz = temp_block;
    *threads_sz = temp_thread;
}

void check_error(cudaError_t err, char *str=NULL)
{
    if (err != cudaSuccess)
    {
        if(str!=NULL){
            printf("ERROR: %s\nline: %s\n", cudaGetErrorString(err), str);
        }else{
            printf("ERROR: %s\n", cudaGetErrorString(err));

        }
        exit(-1);
    }
}

__global__ void print_node(Octtree *tree, int i){
    printf("id: %d, ", i);
        printf("pos: [x: %lf, y: %lf, z:%lf], ", tree->center[i*3], tree->center[i*3+1], tree->center[i*3+2]);
        printf("children: [");
        for(int j=i*8; j<i*8+8; j++){
            printf("%d, ", tree->children[j]);
        }
        printf("], parent: %d\n", tree->parent[i]);
}

__global__ void set_tree(Octtree *tree)
{
    int root = tree->root;
    switch (threadIdx.y)
    {
    case 0:
        if (threadIdx.x == 0)
        {
            tree->firstfree = root + 1;
        }
        break;
    case 1:
        if (threadIdx.x < 3)
        {
            tree->center[root * 3 + threadIdx.x] = 0;
        }
        break;
    case 2:
        if (threadIdx.x == 0)
        {
            tree->mass[root] = 0;
        }
        break;

    case 3:
        if (threadIdx.x == 0)
        {
            tree->ents[root] = 0;
        }
        break;
    case 4:
        if (threadIdx.x == 0)
        {
            tree->parent[root] = -1;
        }
        break;

    case 5:
        if (threadIdx.x < 8)
        {
            tree->children[root * 8 + threadIdx.x] = -1;
        }
        break;
    case 6:
        if (threadIdx.x == 0)
        {
            // tree->depth[root] = 0;
        }
    default:
        break;
    }
}

__global__ void ci_sono_tutti_i_numeri(int *sorted_nodes, int ents_sz)
{
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    __shared__ int cache[1024];

    int cache_id = id % 1024;

    int miss = id<ents_sz ? 1 : 0;
    for (int i = 0; i < ents_sz; i += blockDim.x)
    {
        if (i + cache_id < ents_sz)
        {
            cache[cache_id] = sorted_nodes[i + cache_id];
        }
        __syncthreads();

        int iteration_end = ents_sz > i + blockDim.x ? blockDim.x : ents_sz - i;
        for (int j = 0; j < iteration_end; j++)
        {
            miss = miss && id != cache[j];
        }
        if (__syncthreads_and(!miss))
        {
            return;
        }
    }

    if (miss)
    {
        printf("Missing: %d\n", id); // qua dovrebbe arrivare solo se completa il giro
    }
}

__global__ void find_dups(int *sorted_nodes, int ents_sz)
{
    int id = threadIdx.x + blockDim.x * blockIdx.x;
    if (ents_sz > id)
    {
        int my_val=sorted_nodes[id];
        __shared__ int cache[1024];
        int cache_id = id % 1024;
        for (int i = 0; i < ents_sz; i+=blockDim.x)
        {
            if (i + cache_id < ents_sz)
            {
                cache[cache_id] = sorted_nodes[i + cache_id];
            }
            __syncthreads();
            int iteration_end = ents_sz > i + blockDim.x ? blockDim.x : ents_sz - i;
            for (int j = 0; j < iteration_end; j++)
            {
                if (my_val == cache[j] && id != i+j)
                {
                    printf("Duplicate of %d at pos %d found at %d\n", my_val, id, i+j);
                }
            }
            __syncthreads();
        }
    }
}

__global__ void update_velocities(double *acc, Entities *ents, uint ents_sz,
                                  double dt) {
    int myId;

    myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        ents->vel[myId * 3] += acc[myId * 3] * dt / 2.0;
        ents->vel[myId * 3 + 1] += acc[myId * 3 + 1] * dt / 2.0;
        ents->vel[myId * 3 + 2] += acc[myId * 3 + 2] * dt / 2.0;
    }
}

__global__ void update_positions(Entities *ents, uint ents_sz,
                                 double dt) {
    double3 entPos;
    double3 entVelocities;
    int myId;
    myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        entPos.x=ents->pos[myId*3];
        entPos.y=ents->pos[myId*3+1];
        entPos.z=ents->pos[myId*3+2];
        entVelocities.x=ents->vel[myId*3];
        entVelocities.y=ents->vel[myId*3+1];
        entVelocities.z=ents->vel[myId*3+2];


        ents->pos[myId*3] = entPos.x + entVelocities.x * dt;
        ents->pos[myId*3+1] = entPos.y + entVelocities.y * dt;
        ents->pos[myId*3+2] = entPos.z + entVelocities.z * dt;
    }
}

void print_tree_rec(Octtree *tree, int id, char *space, uint depth) {
    // Octnode *node = &tree->nodes[id];
    // how much divide
    uint temp = 1 << depth;
    double border = (tree->max) / (double)temp;
    printf("%sid: %d, (x:%lf, y:%lf, z:%lf), border: %lf, ents: %d, mass: %lf\n", space, id,
           tree->center[id*3], tree->center[id*3+1], tree->center[id*3+2], border, tree->ents[id], tree->mass[id]);
    if (tree->ents[id] > 1) {

        int i;
        for (i = depth * 4; i < depth * 4 + 4; i++) {
            space[i] = ' ';
        }
        space[i] = '\0';
        for (int i = 0; i < 8; i++) {
            if (tree->children[id*8+i] > -1) {
                print_tree_rec(tree, tree->children[id*8+i], space, depth + 1);
            }
        }
        space[depth * 4] = '\0';
    }
}

void print_tree(Octtree *tree) {
    uint sz_space = 4 * 40;
    char *space = (char*) malloc(sz_space * sizeof(char));
    space[0] = '\0';
    print_tree_rec(tree, tree->root, space, 0);
    free(space);
}

void print_csv(FILE *f, double *d_epos, double *d_emass, Entities *buff, int ents_sz){
    cudaError_t cuda_err;
    cuda_err=cudaMemcpy(buff->pos, d_epos, sizeof(double)*3*ents_sz, cudaMemcpyDeviceToHost);
    check_error(cuda_err, (char *)"print -> copy pos");
    // cudaMemcpy(buff->vel, d_ents->vel, sizeof(double)*3*ents_sz, cudaMemcpyDeviceToHost);
    cuda_err=cudaMemcpy(buff->mass, d_emass, sizeof(double)*ents_sz, cudaMemcpyDeviceToHost);
    check_error(cuda_err, (char *)"print -> copy mass");
    for (size_t i = 0; i < ents_sz; i++)
    {
        double3 pos;
        double mass;
        pos.x=buff->pos[i*3];
        pos.y=buff->pos[i*3+1];
        pos.z=buff->pos[i*3+2];
        mass=buff->mass[i];

        fprintf(f, "%lu,%lf,%lf,%lf,%lf\n", i, pos.x, pos.y, pos.z, mass);
    }
    

}

void print_time(struct timespec *s, struct timespec *e){
    const double BILLION = 1000000000.0;
    double time_spent =
    (e->tv_sec - s->tv_sec) + (e->tv_nsec - s->tv_nsec) / BILLION;

    printf("Elapsed wall time: %f s\n", time_spent);
}


int main(int argc, char *argv[])
{
    // opt_thread is value to optimize
    // the quantity of memory needed for each node
    //                   pos                  mass             children          parent        ents_sz        depth
    const int node_mem = sizeof(double) * 3 + sizeof(double) + sizeof(int) * 8 + sizeof(int) + sizeof(uint) + sizeof(int) + sizeof(uint);
    cudaDeviceProp cuda_prop;
    cudaError_t cuda_err;
    uint n_ents, *d_tents, opt_thread, opt_block;
    dim3 block;
    int *d_tchildren, *d_tparent, *d_sorted_nodes;
    size_t free_mem, total_mem, n_steps;
    #ifdef PRINT_KERNEL_TIME
    struct timespec s, e;
    #endif
    struct timespec s_total, e_total;
    struct cudaFuncAttributes funcAttrib;
    // *_e* memory for entity data
    // *_t* memory for tree data
    // h_le* (host) memory locked for entity data copies
    double *d_epos, *d_evel, *d_emass, *d_tcenter, *d_tmass, *h_lepos, *h_level, *d_reduce1, *d_reduce2, *d_acc;
    Entities h_ents_struct, *d_ents_struct, h_ents_cpy;
    Octtree h_tree_cpy, *d_tree;
    float start, end, dt;
    // if (argc < 6 || argc > 7)
    // {
    //     fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename [cache_sz_MB]\n", argv[0]);
    //     return 1;
    // }
    if(argc != 6){
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
        return -1;
    }

    n_ents = get_entities(argv[1], &h_ents_struct);
    start = strtof(argv[2], NULL);
    end = strtof(argv[3], NULL);
    dt = strtof(argv[4], NULL);
    n_steps = (end - start) / dt;

    cudaGetDeviceProperties(&cuda_prop, 0);

    int nodes_sz = 0; 
    check_error(cudaFuncGetAttributes(&funcAttrib, add_ent), (char *)"get regs");

    cuda_err = cudaMallocHost(&h_lepos, sizeof(double) * n_ents * 3);
    check_error(cuda_err);
    cuda_err = cudaMallocHost(&h_level, sizeof(double) * n_ents * 3);
    check_error(cuda_err);
    cuda_err = cudaMemcpy(h_lepos, h_ents_struct.pos, sizeof(double) * n_ents * 3, cudaMemcpyHostToHost);
    check_error(cuda_err);
    cuda_err = cudaMemcpy(h_level, h_ents_struct.vel, sizeof(double) * n_ents * 3, cudaMemcpyHostToHost);
    check_error(cuda_err);

    cuda_err = cudaMalloc(&d_ents_struct, sizeof(Entities));
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_tree, sizeof(Octtree));
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_epos, sizeof(double) * n_ents * 3);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_evel, sizeof(double) * n_ents * 3);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_emass, sizeof(double) * n_ents);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_sorted_nodes, sizeof(int) * n_ents);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_acc, sizeof(double) * 3 * n_ents);
    check_error(cuda_err);

    cuda_err = cudaMalloc(&d_reduce1, sizeof(double) * (n_ents * 3));
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_reduce2, sizeof(double) * ((n_ents * 3 - 1) / 1024 + 1));
    check_error(cuda_err);

    cuda_err = cudaMemGetInfo(&free_mem, &total_mem);
    check_error(cuda_err);
    // nodes_sz = free_mem * 3 / 4 / node_mem;
    // add int to allocate an array for reordered
    nodes_sz = free_mem * 3 / 4 / (node_mem + sizeof(int));
    printf("FREE MEM: %lu\nLEN NODES: %d\n", free_mem, nodes_sz);
    // cuda_err = cudaMalloc(&d_sorted_nodes, sizeof(int) * nodes_sz);
    cuda_err = cudaMalloc(&d_tcenter, sizeof(double) * nodes_sz * 3);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_tmass, sizeof(double) * nodes_sz);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_tents, sizeof(uint) * nodes_sz);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_tparent, sizeof(uint) * nodes_sz);
    check_error(cuda_err);
    cuda_err = cudaMalloc(&d_tchildren, sizeof(uint) * nodes_sz * 8);
    check_error(cuda_err);

    // initialize struct to be copied in vram
    h_ents_cpy.pos = d_epos;
    h_ents_cpy.vel = d_evel;
    h_ents_cpy.mass = d_emass;
    h_tree_cpy.max = 0;
    h_tree_cpy.ents = d_tents;
    h_tree_cpy.mass = d_tmass;
    h_tree_cpy.center = d_tcenter;
    h_tree_cpy.parent = d_tparent;
    h_tree_cpy.children = d_tchildren;
    h_tree_cpy.root = n_ents;

    cuda_err = cudaMemcpy(d_tree, &h_tree_cpy, sizeof(Octtree), cudaMemcpyHostToDevice);
    check_error(cuda_err);
    cuda_err = cudaMemcpy(d_ents_struct, &h_ents_cpy, sizeof(Entities), cudaMemcpyHostToDevice);
    check_error(cuda_err);
    cuda_err = cudaMemcpy(d_epos, h_lepos, sizeof(double) * n_ents * 3, cudaMemcpyHostToDevice);
    check_error(cuda_err);
    cuda_err = cudaMemcpy(d_evel, h_level, sizeof(double) * n_ents * 3, cudaMemcpyHostToDevice);
    check_error(cuda_err);
    cuda_err = cudaMemcpy(d_emass, h_ents_struct.mass, sizeof(double) * n_ents, cudaMemcpyHostToDevice);
    check_error(cuda_err);

    //pos and mass will be used like a buffer for the prints
    free(h_ents_struct.vel);
    
    int max_threads = cuda_prop.maxThreadsPerBlock;
    printf("Initialization completed\n");
    FILE *fpt = fopen(argv[5], "w");
    
    clock_gettime(CLOCK_REALTIME, &s_total);
    double *d_reduce_in, *d_reduce_out;
    block.x = cuda_prop.warpSize;
    block.y = 7;
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &s);
    #endif
    set_tree<<<1, block>>>(d_tree);
    cuda_err = cudaDeviceSynchronize();
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &e);
    print_time(&s, &e);
    #endif
    
    check_error(cuda_err);
    d_reduce_in = d_reduce1;
    d_reduce_out = d_reduce2;
    cudaMemcpy(d_reduce_in, d_epos, sizeof(double)*n_ents*3, cudaMemcpyDeviceToDevice);
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &s);
    #endif
    for(int i=n_ents*3; i>1;i=((i-1)/(max_threads*2)+1)){
        get_bounding_box<<<(i-1)/(max_threads*2)+1, max_threads, max_threads * sizeof(double)>>>(d_reduce_in, i, d_reduce_out);
        cuda_err = cudaDeviceSynchronize();
        check_error(cuda_err);
        double *temp=d_reduce_out;
        d_reduce_out=d_reduce_in;
        d_reduce_in=temp;
    }
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &e);
    print_time(&s, &e);
    #endif
    cudaMemcpy(&d_tree->max, d_reduce_in, sizeof(double), cudaMemcpyDeviceToDevice);
    get_opt_grid(&cuda_prop, n_ents, funcAttrib.numRegs, &opt_block, &opt_thread);
    printf("opt_grid: %u %u\n", opt_block, opt_thread);
    // uses 56 registers

    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &s);
    #endif
    add_ent<<<opt_block, opt_thread>>>(d_tree, d_ents_struct, n_ents);
    cuda_err = cudaDeviceSynchronize();
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &e);
    print_time(&s, &e);
    #endif
    check_error(cuda_err);

    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &s);
    #endif
    center_of_mass<<<opt_block, opt_thread>>>(d_tree);
    cuda_err = cudaDeviceSynchronize();
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &e);
    print_time(&s, &e);
    #endif
    check_error(cuda_err);

    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &s);
    #endif
    sort_ents<<<(n_ents-1)/max_threads+1, max_threads>>>(d_tree, d_sorted_nodes);
    cuda_err = cudaDeviceSynchronize();
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &e);
    print_time(&s, &e);
    #endif
    check_error(cuda_err, (char *)"sort_ents call");


    cudaMemset(d_acc, 0, sizeof(double)*3*n_ents);
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &s);
    #endif
    
    acceleration_w_stack<<<(n_ents-1)/(cuda_prop.warpSize*2)+1,cuda_prop.warpSize*2, cuda_prop.sharedMemPerBlock>>>(d_tree, d_ents_struct, n_ents, d_sorted_nodes, d_acc, cuda_prop.sharedMemPerBlock);
    cuda_err = cudaDeviceSynchronize();
    #ifdef PRINT_KERNEL_TIME
    clock_gettime(CLOCK_REALTIME, &e);
    print_time(&s, &e);
    #endif
    check_error(cuda_err, (char *)"acceleration line");


    int grid_sz=(n_ents-1)/max_threads+1;
    printf("STARTING LOOP\n");
    printf("number of steps: %lu\n", n_steps);
    for(size_t t = 1; t <= n_steps; t++){
        #ifdef PRINT_LOOP
        printf("loop %d\n", t);
        #endif
         #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        update_velocities<<<grid_sz, max_threads>>>(d_acc, d_ents_struct,
                                        n_ents, dt);
        cuda_err=cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"update vel loop");

        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        update_positions<<<grid_sz, max_threads>>>(d_ents_struct, n_ents,
                                        dt);
        cuda_err=cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"update pos loop");
#ifdef RESULTS
        print_csv(fpt, d_epos, d_emass, &h_ents_struct, n_ents);
#endif

        block.x = cuda_prop.warpSize;
        block.y = 7;
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        set_tree<<<1, block>>>(d_tree);
        cuda_err=cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"set tree loop");

        d_reduce_in = d_reduce1;
        d_reduce_out = d_reduce2;
        
        cuda_err=cudaMemcpy(d_reduce_in, d_epos, sizeof(double)*n_ents*3, cudaMemcpyDeviceToDevice);
        check_error(cuda_err, (char *)"copy d_epos to d_reduce_in loop");
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        for(int i=n_ents*3; i>1;i=((i-1)/(max_threads*2)+1)){
            get_bounding_box<<<(i-1)/(max_threads*2)+1, max_threads, max_threads * sizeof(double)>>>(d_reduce_in, i, d_reduce_out);
            cuda_err = cudaDeviceSynchronize();
            check_error(cuda_err);
            double *temp=d_reduce_out;
            d_reduce_out=d_reduce_in;
            d_reduce_in=temp;
        }
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif

        cuda_err=cudaMemcpy(&d_tree->max, d_reduce_in, sizeof(double), cudaMemcpyDeviceToDevice);
        check_error(cuda_err, (char *)"copy d_reduce_in to d_tree.max loop");

        get_opt_grid(&cuda_prop, n_ents, 56, &opt_block, &opt_thread);
        
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        add_ent<<<opt_block, opt_thread>>>(d_tree, d_ents_struct, n_ents);
        cuda_err = cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"add ent loop");

        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        center_of_mass<<<opt_block, opt_thread>>>(d_tree);
        cuda_err = cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"center of mass loop");

        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        sort_ents<<<(n_ents-1)/max_threads+1, max_threads>>>(d_tree, d_sorted_nodes);
        cuda_err = cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"sort ents loop");
        
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        acceleration_w_stack<<<(n_ents-1)/(cuda_prop.warpSize*2)+1,cuda_prop.warpSize*2, cuda_prop.sharedMemPerBlock>>>(d_tree, d_ents_struct, n_ents, d_sorted_nodes, d_acc, cuda_prop.sharedMemPerBlock);
        cuda_err = cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"acceleration loop");

        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &s);
        #endif
        update_velocities<<<grid_sz, max_threads>>>(d_acc, d_ents_struct,
                                        n_ents, dt);
        cuda_err=cudaDeviceSynchronize();
        #ifdef PRINT_KERNEL_TIME
        clock_gettime(CLOCK_REALTIME, &e);
        print_time(&s, &e);
        #endif
        check_error(cuda_err, (char *)"update vel 2 loop");
    }
    clock_gettime(CLOCK_REALTIME, &e_total);
    // print_csv(fpt, d_epos, d_emass, &h_ents_struct, n_ents);
    print_time(&s_total, &e_total);
        

    fclose(fpt);
    cudaFreeHost(h_lepos);
    cudaFreeHost(h_level);
    cudaFree(d_acc);
    cudaFree(d_ents_struct);
    cudaFree(d_tree);
    cudaFree(d_epos);
    cudaFree(d_evel);
    cudaFree(d_emass);
    cudaFree(d_tcenter);
    cudaFree(d_tmass);
    cudaFree(d_tents);
    cudaFree(d_tparent);
    cudaFree(d_tchildren);
    cudaFree(d_sorted_nodes);
    cudaFree(d_reduce1);
    cudaFree(d_reduce2);
    free(h_ents_struct.pos);
    free(h_ents_struct.mass);
}
