#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#define LOCKED -2

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
    int *depth;

} Octtree;

const double BIG_G = 6.67e-11;
const double THETA = 0.5; // 1;
// TODO check from version

// a lot of optimization like use shift instead of multiplication will be made by compiler

// TODO put in a common file
uint get_entities(char filename[], Entities *ents)
{
    // Entity e_buff;
    double pos_buff[3];
    double vel_buff[3];
    double mass_buff;
    int status;
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

    // TODO Check for error in allocation
    pos_ret = (double *)malloc(3 * sizeof(double));
    vel_ret = (double *)malloc(3 * sizeof(RVec3));
    mass_ret = (double *)malloc(1 * sizeof(double));

    size = 0;
    ret_size = 1;
    // fscanf return the number of input items successfully matched and assigned
    while (
        (status = fscanf(f, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &pos_buff[0],
                         &pos_buff[1], &pos_buff[2], &vel_buff[0], &vel_buff[1],
                         &vel_buff[2], &mass_buff)) == 7)
    {
        size++;
        if (ret_size < size)
        {
            // TODO Check for error in allocation
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
    tree->depth[indx] = depth;
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

//KERNEL 2
// add a entity in the tree
// it's create all the needed branch
// the leaf of the tree are biunivocaly (e' inglese?) associated to an entity
__global__ void add_ent(Octtree *tree, Entities *ent, uint ents_sz)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    if (id < ents_sz)
    {

        double pos_ent[3];
        // allocated is used as a boolean
        int allocated, node_indx, ent_loc, child_indx, child_val, other, curr_depth;
        double border_size;
        // Octnode node;
        // set center of whole volume
        double volume_center[3] = {0, 0, 0};
        // set init value
        allocated = 0;

        copy_arrs3(pos_ent, &ent->pos[id * 3], 0);

        // keep last visited node index
        node_indx = tree->root;

        border_size = border_tree(tree);

        curr_depth=0;

        while (!allocated)
        {
            // center and border_size are updated to the next branch value
            ent_loc = get_indx_loc(pos_ent, volume_center, &border_size);
            // the position where entity should be placed in children array
            child_indx = node_indx * 8 + ent_loc;
            // current value in the child_indx location

            // temporany solution
            do
            {
                child_val = tree->children[child_indx];
                // atomic access to child
                if (child_val != LOCKED)
                {
                    if (child_val ==
                        atomicCAS(&tree->children[child_indx], child_val, LOCKED))
                    {
                        curr_depth++;
                        // if nothing is located, allocate the leaf
                        if (child_val == -1)
                        {
                            tree->children[child_indx] = id;
                            allocated = 1;
                            __threadfence();
                            // if there is a leaf, start to divide until the leaves will
                            // be allocated in 2 different place in children array
                        }
                        // if there is already a leaf
                        else if (child_val < tree->root)
                        {
                            allocated = 0;
                            // the leaf that was already there
                            other = child_val;
                            int other_loc;
                            double other_center[3];
                            copy_arrs3(other_center, volume_center, 0);
                            double other_border = border_size;
                            int new_node;
                            // node_indx is used to store last node created
                            while (!allocated)
                            {
                                // new node location
                                // atomic add return old value
                                new_node = atomicAdd(&tree->firstfree, 1);
                                // do{
                                //     new_node=tree->firstfree;
                                // }while(new_node==LOCKED || new_node!=atomicCAS(&tree->firstfree, new_node, LOCKED));
                                // tree->firstfree=new_node+1;

                                init_node(tree, new_node, curr_depth);
                                tree->parent[new_node] = node_indx;
                                curr_depth++;

                                // get leaf position in the new branch
                                ent_loc = get_indx_loc(pos_ent, volume_center,
                                                       &border_size);

                                // TODO load once tree->center[other * 3]
                                //  the location of the previous children
                                other_loc = get_indx_loc(&tree->center[other * 3],
                                                         other_center, &other_border);

                                allocated = other_loc != ent_loc;

                                tree->children[new_node * 8 + ent_loc] = LOCKED;
                                // not execute in the last loop
                                // TODO optimize
                                if (allocated)
                                {
                                    // lock location that will be explored
                                    tree->children[new_node * 8 + other_loc] = LOCKED;
                                }
                                // child_indx is not updated, so we can use it to
                                // get the location of new_node in children array
                                tree->children[child_indx] = new_node;
                                __threadfence();
                                // slot in children where will be created new node or will be put the body
                                child_indx = new_node * 8 + ent_loc;

                                node_indx = new_node;
                            }

                            tree->children[child_indx] = id;
                            tree->children[node_indx * 8 + other_loc] = other;
                            // already done at the end
                            // tree->parent[id] = node_indx;
                            tree->parent[other] = node_indx;
                            __threadfence();
                            // if there is a branch
                        }
                        else if (child_val >= tree->root)
                        {
                            node_indx = child_val;
                            // unlock child setting the value before lock
                            tree->children[child_indx] = child_val;
                            __threadfence();
                        }
                        else
                        {
                            // ERRORs
                        }
                    }
                    else
                    {
                        child_val = LOCKED;
                    }
                }
            } while (child_val == LOCKED);
        }
        // set the leaf value
        // set children should be not needed, but i'm not sure
        copy_arrs3(&tree->center[id * 3], pos_ent, 0);
        tree->mass[id] = ent->mass[id];
        tree->ents[id] = 1;
        tree->parent[id] = node_indx;
        tree->depth[id] = curr_depth;
    }
}

//KERNEL 1
// TODO write a more optimized function
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
    if (i < ents_sz * 3)
    {
        v1 = fabs(g_idata[i]);
        if (i + space < ents_sz)
        {
            v2 = fabs(g_idata[i + space]);
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

//KERNEL 3
// if mass==0 node is not ready
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
                if (tree->mass[child_indx] == 0)
                {
                    children_cache[cache_sz] = child_indx;
                    cache_sz++;
                }
                else
                {
                    new_mass += tree->mass[child_indx];
                    ents += 1;
                    // TODO optimized, just divide at the end
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
        // TODO resolve divergence
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
                    ents += 1;
                    // TODO optimized, just divide at the end
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
                tree->ents[indx] = ents;
                tree->mass[indx] = mass;
                __threadfence();
            }
        } while (cache_sz > 0);
    }
}

//KERNEL 4
// here tree->parent will be used to store the "s", the position where to write values (tree->parent is not needed anymore)
__global__ void order_ents(Octtree *tree, int *sorted)
{
    uint id = threadIdx.x + blockIdx.x * blockDim.x;
    uint stride = blockDim.x * gridDim.x;
    uint offset = 0;
    const int root = tree->root;

    // s is the position in the array where to write the body
    int s = 0;
    // TODO check if there is a better way to wait that values have been loaded
    if (threadIdx.x == 0)
    {
        for (int i = 0; i < 8; i++)
        {
            int node = tree->children[root * 8 + i];

            if (node >= root)
            { // not a leaf node
                tree->parent[node] = s;
                s += tree->ents[node];
            }
            else if (node >= 0)
            { // leaf node
                sorted[s] = node;
                s++;
            }
        }
        __threadfence();
    }
    __syncthreads();

    int cell = root + id + 1;
    int tree_sz = tree->firstfree;
    while ((cell + offset) < tree_sz)
    {
        s = tree->parent[cell + offset];

        if (s >= 0)
        {

            for (int i = 0; i < 8; i++)
            {
                int node = tree->children[8 * (cell + offset) + i];

                if (node >= root)
                { // not a leaf node
                    tree->parent[node] = s;
                    s += tree->ents[node];
                }
                else if (node >= 0)
                { // leaf node
                    sorted[s] = node;
                    s++;
                }
            }
            offset += stride;
            __threadfence();
        }
    }
}

//il paper mette il calcolo delle accelerazioni nel kernel 5 e l'aggiornamento e velocita' nel kernel 6
__device__ void calculate_acceleration(double *pos_node, double mass_node, double *pos_ent, double mass_ent,
                                       double acc[])
{
    // RVec3 r_vector;
    double r_vector[3], r_unit_vector[3];
    r_vector[0] = pos_ent[0] - pos_node[0];
    r_vector[1] = pos_ent[1] - pos_node[1];
    r_vector[2] = pos_ent[2] - pos_node[2];

    double r_mag = sqrt(r_vector[0] * r_vector[0] + r_vector[1] * r_vector[1] +
                        r_vector[2] * r_vector[2]);

    double acceleration = -1.0 * BIG_G * (mass_node) / (r_mag * r_mag);

    r_unit_vector[0] = r_vector[0] / r_mag;
    r_unit_vector[1] = r_vector[1] / r_mag;
    r_unit_vector[2] = r_vector[2] / r_mag;

    acc[0] += acceleration * r_unit_vector[0];
    acc[1] += acceleration * r_unit_vector[1];
    acc[2] += acceleration * r_unit_vector[2];
}

__device__ double get_distance(double *r1, double *r2)
{
    return sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]) + (r1[2] - r2[2]) * (r1[2] - r2[2]));
}

//KERNEL 5
//KERNEL 6
// TODO children must be compacted at the start of array
__global__ void get_acceleration(Octtree *tree, Entities *ents, int ents_sz, size_t dt, uint shared_sz)
{
    uint threads_sz = gridDim.x * blockDim.x * blockDim.y;
    uint tid;
    // 1024*48 is the total shared memory, remove space taken from double and from cache_node, divide for 32 (number of warps)
    int id_in_warp, curr_node, depth, my_stack, warp_id, stack_level, warps_sz;
    double pos_ent[3], mass_ent, border, acc_ent[3], ddt;
    warps_sz = blockDim.y;
    // extern __shared__ char shared_mem[];
    // double *pos_node = (double *)shared_mem;
    // //first cast to void so to shift in bytes
    // double *mass_node = (double *)(((void *)pos_node) + 3 * warps_sz * sizeof(double));
    // int *cache_node = (int *)(((void *)mass_node) + warps_sz * sizeof(double));
    // int *stack = (int *)(((void *)cache_node) + warps_sz * sizeof(int));
    __shared__ double pos_node[3 * 2];
    __shared__ double mass_node[2];
    __shared__ int cache_node[2];
    __shared__ int depth_node[2];
    __shared__ int stack[3050 / 4 + 1];

    const int cache_sz = shared_sz - 3 * warps_sz * sizeof(double) - warps_sz * sizeof(double) - warps_sz * sizeof(int);
    // TODO adapt for more block
    tid = (threadIdx.x + threadIdx.y * blockDim.x) + blockIdx.x * gridDim.x;
    id_in_warp = threadIdx.x;
    border = tree->max * 2;
    warp_id = threadIdx.y;
    // my_stack = &stack[cache_sz / warps_sz * warp_id];
    my_stack = 3050 / 4 / 2 * warps_sz +1;
    ddt = (double)dt;

    // root is the first node located after the leaves, so each id must be an id of a leaf
    for (int i = tid; i < ents_sz; i += threads_sz)
    {
        // initialize values
        pos_ent[0] = ents->pos[i * 3];
        pos_ent[1] = ents->pos[i * 3 + 1];
        pos_ent[2] = ents->pos[i * 3 + 2];
        mass_ent = ents->mass[i];
        acc_ent[0] = 0;
        acc_ent[1] = 0;
        acc_ent[2] = 0;
        stack_level = 0;

        depth = 0;
        curr_node = tree->root;

        // start iterate algorithm, when depth<0 the algorithm has ended
        while (stack_level >= 0)
        {
            for (int j = 0; j < 8; j++)
            {
                // save in shared memory the node
                // TODO load only id and break if <0
                if (id_in_warp == 0)
                {
                    cache_node[warp_id] = tree->children[curr_node * 8 + j];
                    
                }
                // be sure that every thread of block (but we are interested in thread in the same warp) can see the node
                __threadfence_block();
                if (cache_node[warp_id] >= 0)
                {
                    if (cache_node[warp_id] >= 0)
                    {
                        pos_node[warp_id * 3] = tree->center[cache_node[warp_id] * 3];
                        pos_node[warp_id * 3 + 1] = tree->center[cache_node[warp_id] * 3 + 1];
                        pos_node[warp_id * 3 + 2] = tree->center[cache_node[warp_id] * 3 + 2];
                        mass_node[warp_id] = tree->mass[cache_node[warp_id]];
                        depth_node[warp_id] = tree->depth[cache_node[warp_id]];
                    }
                    if (cache_node[warp_id] != i)
                    {
                        double distance = get_distance(pos_ent, pos_node);
                        // if the node is a leaf or all thread in warp have their entity far enough from the body
                        if (cache_node[warp_id] < tree->root || __all_sync(0xFFFFFFFF, (border / (1<<depth_node[warp_id]))/distance < THETA))
                        {
                            //probabilmente non calcola il sole
                            calculate_acceleration(&pos_node[warp_id * 3], mass_node[warp_id], pos_ent, mass_ent, acc_ent);
                        }
                        else
                        {
                            // save node in stack, when node will be taken from the stack, the programm will check its children
                            depth++;
                            if (id_in_warp == 0)
                            {
                                // my_stack[stack_level] = cache_node[warp_id];
                                stack[my_stack + stack_level] = cache_node[warp_id];
                            }
                            stack_level++;
                            __threadfence_block();
                        }
                    }
                }
                else
                {
                    depth = max(0, depth - 1);
                }
            }
            depth--;
            stack_level--;
            // curr_node = my_stack[stack_level];
            curr_node = stack[my_stack + stack_level];
        }
        // update velocities
        ents->vel[i * 3] += acc_ent[0] * ddt;
        ents->vel[i * 3 + 1] += acc_ent[1] * ddt;
        ents->vel[i * 3 + 2] += acc_ent[2] * ddt;
    }

    // maybe that can be included in previous loop
    for (int i = tid; i < ents_sz; i += threads_sz)
    {
        ents->pos[i * 3] += ents->vel[i * 3] * ddt;
        ents->pos[i * 3 + 1] += ents->vel[i * 3 + 1] * ddt;
        ents->pos[i * 3 + 2] += ents->vel[i * 3 + 2] * ddt;
    }
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
        if(threadIdx.x==0){
            tree->depth[root]=0;
        }
    default:
        break;
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

void check_error(cudaError_t err)
{
    if (err != cudaSuccess)
    {
        printf("ERROR: %s\n", cudaGetErrorString(err));
        exit(-1);
    }
}

// TODO check errors from malloc and kernels and fopen
// TODO implement cache for calculated value
int main(int argc, char *argv[])
{
    // opt_thread is value to optimize
    // the quantity of memory needed for each node
    //                   pos                  mass             children          parent        ents_sz        depth
    const int node_mem = sizeof(double) * 3 + sizeof(double) + sizeof(int) * 8 + sizeof(int) + sizeof(uint) + sizeof(int);
    cudaDeviceProp cuda_prop;
    cudaError_t cuda_err;
    uint n_ents, *d_tents, opt_thread, opt_block, cache_sz;
    dim3 block;
    int *d_tchildren, *d_tparent, *d_sorted_nodes, *d_depth;
    size_t free_mem, total_mem;
    // *_e* memory for entity data
    // *_t* memory for tree data
    // h_le* (host) memory locked for entity data copies
    double *d_epos, *d_evel, *d_emass, *d_tcenter, *d_tmass, *h_lepos, *h_level, *d_reduce1, *d_reduce2;
    Entities h_ents_struct, *d_ents_struct, h_ents_cpy;
    Octtree h_tree_cpy, *d_tree;
    size_t start, end, dt;
    cache_sz = 0;
    // if (argc < 6 || argc > 7)
    // {
    //     fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename [cache_sz_MB]\n", argv[0]);
    //     return 1;
    // }

    // n_ents = get_entities(argv[1], &h_ents_struct);
    // start = strtoul(argv[2], NULL, 10);
    // end = strtoul(argv[3], NULL, 10);
    // dt = strtoul(argv[4], NULL, 10);
    // if (argc == 7)
    // {
    //     cache_sz = strtoul(argv[6], NULL, 10);
    // }
    n_ents = get_entities("/home/prop/Documents/multicore/N-Body-simulation/tests/sun_earth.csv", &h_ents_struct);
    start = 0;
    end = 86400 * 500;
    dt = 86400;

    cudaGetDeviceProperties(&cuda_prop, 0);

    int nodes_sz = 0; // TODO get size of number of node

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
    check_error(cuda_err);
    if (n_ents > 1024 * 2)
    {
        cuda_err = cudaMalloc(&d_reduce1, sizeof(double) * ((n_ents * 3 - 1) / 1024 + 1));
        check_error(cuda_err);
    }
    else
    {
        // the final destination of result is the tree attribute max
        // need to be defined also here because first reduction is out of the loop (if loop will be executed, the program will not enter this branch)
        d_reduce1 = &d_tree->max;
    }
    if (n_ents > 1024 * 2 * 1024 * 2)
    {
        cuda_err = cudaMalloc(&d_reduce2, sizeof(double) * ((n_ents * 3 - 1) / 1024 / 1024 + 1));
        check_error(cuda_err);
    }

    cuda_err = cudaMemGetInfo(&free_mem, &total_mem);
    check_error(cuda_err);
    // nodes_sz = free_mem * 3 / 4 / node_mem;
    // add int to allocate an array for reordered
    nodes_sz = free_mem * 3 / 4 / (node_mem + sizeof(int));
    cuda_err = cudaMalloc(&d_sorted_nodes, sizeof(int) * nodes_sz);
    check_error(cuda_err);
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
    cuda_err = cudaMalloc(&d_depth, sizeof(int)*nodes_sz);

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
    h_tree_cpy.depth = d_depth;

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

    free(h_ents_struct.pos);
    free(h_ents_struct.vel);
    free(h_ents_struct.mass);
    // TODO recalculate thread and block size
    int max_threads = cuda_prop.maxThreadsPerBlock;
    printf("Initialization completed\n");
    // FILE *fpt = fopen(argv[5], "w");
    FILE *fpt = fopen("/home/prop/Documents/multicore/N-Body-simulation/tests/output/cuda_barn.csv", "w");
    for (size_t t = start; t < end; t += dt)
    {
        double *d_reduce_in, *d_reduce_out, *temp_swap;
        block.x = cuda_prop.warpSize;
        block.y = 7;
        set_tree<<<1, block>>>(d_tree);
        cuda_err = cudaDeviceSynchronize();
        check_error(cuda_err);
        get_bounding_box<<<(n_ents - 1) / (max_threads * 2) + 1, max_threads, max_threads * sizeof(double)>>>(d_epos, n_ents, d_reduce1);
        cuda_err = cudaDeviceSynchronize();
        check_error(cuda_err);
        d_reduce_in = d_reduce2;
        d_reduce_out = d_reduce1;
        // TODO check if num of blocks are correct
        for (int temp_sz = (n_ents - 1) / (max_threads * 2) + 1; temp_sz > 1; temp_sz = (temp_sz - 1) / (max_threads * 2) + 1)
        {
            if (temp_sz <= max_threads * 2)
            {
                temp_swap = d_reduce_in;
                d_reduce_in = d_reduce_out;
                d_reduce_out = temp_swap;
            }
            else
            {
                d_reduce_in = d_reduce_out;
                // the final destination of result is the tree attribute max
                d_reduce_out = &d_tree->max;
            }
            // uses 12 registers
            get_bounding_box<<<(temp_sz - 1) / (max_threads * 2) + 1, max_threads, max_threads * sizeof(double)>>>(d_reduce_in, temp_sz, d_reduce_out);
            cuda_err = cudaDeviceSynchronize();
            check_error(cuda_err);
        }
        get_opt_grid(&cuda_prop, n_ents, 56, &opt_block, &opt_thread);
        // uses 56 registers
        add_ent<<<opt_block, opt_thread>>>(d_tree, d_ents_struct, n_ents);
        // add_ent<<<1, 128>>>(d_tree, d_ents_struct, n_ents);
        cuda_err = cudaDeviceSynchronize();

        check_error(cuda_err);

        get_opt_grid(&cuda_prop, 0, 62, &opt_block, &opt_thread);
        // uses 62 registers
        set_branch_values<<<opt_block, opt_thread>>>(d_tree);
        // set_branch_values<<<1, 10>>>(d_tree);
        cuda_err = cudaDeviceSynchronize();
        check_error(cuda_err);

        // maybe the max threads could be n_ents or some fraction like n_ents/2
        get_opt_grid(&cuda_prop, 0, 16, &opt_block, &opt_thread);
        // uses 16 registers
        order_ents<<<opt_block, opt_thread>>>(d_tree, d_sorted_nodes);
        cuda_err = cudaDeviceSynchronize();
        check_error(cuda_err);

        get_opt_grid(&cuda_prop, n_ents, 6, &opt_block, &opt_thread);
        // uses 6 registers
        block.x = cuda_prop.warpSize;
        block.y = opt_thread / cuda_prop.warpSize;
        get_acceleration<<<opt_block, block, cuda_prop.sharedMemPerBlock / cuda_prop.maxBlocksPerMultiProcessor>>>(d_tree, d_ents_struct, n_ents, dt, cuda_prop.sharedMemPerBlock);
        cuda_err = cudaDeviceSynchronize();
        check_error(cuda_err);
        cudaMemcpy(h_level, d_evel, sizeof(double) * n_ents * 3, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_lepos, d_epos, sizeof(double) * n_ents * 3, cudaMemcpyDeviceToHost);
        print_values(h_lepos, h_level, n_ents, fpt);
    }

    fclose(fpt);
    cudaFreeHost(h_lepos);
    cudaFreeHost(h_level);
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
    cudaFree(d_depth);
    if (n_ents > 1024 * 2)
    {
        cudaFree(d_reduce1);
    }
    if (n_ents > 1024 * 2 * 1024 * 2)
    {
        cudaFree(d_reduce2);
    }
}
