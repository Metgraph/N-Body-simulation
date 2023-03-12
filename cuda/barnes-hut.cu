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

} Octtree;

const double BIG_G = 6.67e-11;
const double THETA = 0.5; // 1;
// TODO check from version
const int max_thread_block = 1024;

//a lot of optimization like use shift instead of multiplication will be made by compiler

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
            vel_ret = (double *)realloc(pos_ret, ret_size * 3 * sizeof(double));
            mass_ret = (double *)realloc(pos_ret, ret_size * sizeof(double));
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

__device__ void init_node(Octtree *tree, int indx)
{
    tree->center[indx * 3] = 0;
    tree->center[indx * 3 + 1] = 0;
    tree->center[indx * 3 + 2] = 0;
    tree->mass[indx] = 0;
    tree->ents[indx] = 0;
    tree->parent[indx] = -1;
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

// add a entity in the tree
// it's create all the needed branch
// the leaf of the tree are biunivocaly (e' inglese?) associated to an entity
__global__ void add_ent(Octtree *tree, Entities *ent, uint ents_sz)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    if (ents_sz < id)
    {

        double pos_ent[3];
        // allocated is used as a boolean
        int allocated, node_indx, ent_loc, child_indx, child_val, other;
        double border_size;
        // Octnode node;
        // set center of whole volume
        double volume_center[3] = {0, 0, 0};
        // set init value
        allocated = 0;

        // keep last visited node index
        node_indx = tree->root;

        border_size = border_tree(tree);

        while (!allocated)
        {
            // center and border_size are updated to the next branch value
            ent_loc = get_indx_loc(&ent->pos[id * 3], volume_center, &border_size);
            // the position where entity should be placed in children array
            child_indx = node_indx * 8 + ent_loc;
            // current value in the child_indx location
            child_val = tree->children[child_indx];
            if (child_val != LOCKED)
            {
                // only one thread can enter
                if (child_val ==
                    atomicCAS(&tree->children[child_indx], child_val, LOCKED))
                {
                    // if nothing is located, allocate the leaf
                    if (child_val == -1)
                    {
                        tree->children[child_indx] = id;
                        __threadfence();
                        // if there is a leaf, start to divide until the leaves will
                        // be allocated in 2 different place in children array
                    }
                    else if (child_val < tree->root)
                    {
                        allocated = 0;
                        // the leaf that was already there
                        other = child_val;
                        int other_loc;
                        double other_center[3];
                        copy_arrs3(other_center, volume_center, 0);
                        double other_border = border_size;
                        while (!allocated)
                        {
                            // new node location
                            // atomic add return old value
                            int new_node = atomicAdd(&tree->firstfree, 1);

                            init_node(tree, new_node);
                            tree->parent[new_node] = node_indx;

                            // get leaves position in the new branch
                            ent_loc = get_indx_loc(&ent->pos[id * 3], volume_center,
                                                   &border_size);

                            // the location of the previous children
                            other_loc = get_indx_loc(&tree->center[other * 3],
                                                     other_center, &other_border);

                            allocated = other_loc != ent_loc;
                            // not execute in the last loop
                            if (!allocated)
                            {
                                // lock location where
                                tree->children[node_indx * 8 + ent_loc] = LOCKED;
                                // child_indx is not updated, so we can use it to
                                // get the location of new_node in children array
                                tree->children[child_indx] = new_node;
                                __threadfence();
                            }
                            // slot in children where thread is working
                            child_indx = node_indx * 8 + ent_loc;

                            node_indx = new_node;

                            // use the new branch as the current one
                            // indx = get_indx_loc(&tree->center[]);
                            // set the new branch as child
                            // unlock children
                            // update first free location
                        }

                        tree->children[child_indx] = id;
                        tree->children[node_indx * 8 + other_loc] = other;
                        tree->parent[id] = node_indx;
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
            }
        }
        // set the leaf value
        copy_arrs3(tree->center, ent->pos, id * 3);
        tree->mass[id] = ent->mass[id];
        tree->ents[id] = 1;
        tree->parent[id] = node_indx;
    }
}

//TODO write a more optimized function
__global__ void get_bounding_box(double *g_idata, int ents_sz, double *g_odata)
{
    __shared__ double sdata[1024];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint tid = threadIdx.x;
    uint space = (blockDim.x);
    uint id = blockIdx.x * blockDim.x + threadIdx.x;
    uint i = (blockIdx.x * (blockDim.x * 2) + threadIdx.x);

    //if thread are more than values give them a 0 value
    //the first if has no divergence, becuase threads have same ents_sz value
    if(ents_sz%2==0){
        if(id<ents_sz/2){
            sdata[tid] = fabs(g_idata[i]) > fabs(g_idata[i + space]) ? fabs(g_idata[i]) : fabs(g_idata[i + space]);
        }else{
            sdata[tid] = 0;
        }
    }else{
        if(id<ents_sz-1/2+1){
            // the condition after the or is if it is the last element, and since ents_sz is odd, the last thread has only a value instead of 2
            sdata[tid] = fabs(g_idata[i]) > fabs(g_idata[i + space]) || id==ents_sz/2 ? fabs(g_idata[i]) : fabs(g_idata[i + space]);
        }else{
            sdata[tid] = 0;
        }
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

// if mass==0 node is not ready
__global__ void set_branch_values(Octtree *tree)
{
    int sz_threads = gridDim.x * blockDim.x;
    // used as cache
    int children_cache[8];
    double center[3];
    int child_indx, last_node, sz_nodes, ents;
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int cache_sz;
    // new_mass is needed because we need old and new mass value at the same time
    double mass, new_mass;
    last_node = tree->firstfree - 1;
    sz_nodes = last_node - tree->root;
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

        for (int i = 0; i < 8; i++)
        {
            // get the child index
            child_indx = tree->children[indx * 8 + i];
            if (child_indx != -1)
            {
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
        }
        while (cache_sz > 0)
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
        }
        copy_arrs3(&tree->center[indx * 3], center, 0);
        tree->ents[indx] = ents;
        tree->mass[indx] = mass;
        __threadfence();
    }
}

//here tree->parent will be used to store the "s", the position where to write values (tree->parent is not needed anymore)
__global__ void order_ents(Octtree *tree, int *sorted){
    uint id = threadIdx.x + blockIdx.x*blockDim.x;
	uint stride = blockDim.x*gridDim.x;
	uint offset = 0;
    const int root= tree->root;

    //s is the position in the array where to write the body
	int s = 0;
	if(threadIdx.x == 0){
		for(int i=0;i<8;i++){
			int node = tree->children[i*root*8];

			if(node >= root){  // not a leaf node
				tree->parent[node] = s;
				s += tree->ents[node];
			}
			else if(node >= 0){  // leaf node
				sorted[s] = node;
				s++;
			}
		}
        __threadfence();
	}

	int cell = root + id;
	int tree_sz = tree->firstfree;
	while((cell + offset) < tree_sz){
		s = tree->parent[cell + offset];
	
		if(s >= 0){

			for(int i=0;i<8;i++){
				int node = tree->children[8*(cell+offset) + i];

				if(node >= root){  // not a leaf node
					tree->parent[node] = s;
					s += tree->ents[node];
				}
				else if(node >= 0){  // leaf node
					sorted[s] = node;
					s++;
				}
			}
			offset += stride;
            __threadfence();
		}
	}
}


__device__ void calculate_acceleration(double *pos_node, double mass_node, double *pos_ent, double mass_ent,
                                       double acc[])
{
    // RVec3 r_vector;
    double r_vector[3], r_unit_vector[3];
    r_vector[0] = pos_node[0] - pos_ent[0];
    r_vector[1] = pos_node[1] - pos_ent[1];
    r_vector[2] = pos_node[2] - pos_ent[2];

    double r_mag = sqrt(r_vector[0] * r_vector[0] + r_vector[1] * r_vector[1] +
                        r_vector[2] * r_vector[2]);

    double acceleration = -1.0 * BIG_G * (mass_node) / (r_mag * r_mag);

    r_unit_vector[0] = r_vector[0] / r_mag;
    r_unit_vector[1] = r_vector[1] / r_mag;
    r_unit_vector[2] = r_vector[2] / r_mag;

    acc[0] = acceleration * r_unit_vector[0];
    acc[1] = acceleration * r_unit_vector[1];
    acc[2] = acceleration * r_unit_vector[2];
}

__device__ double get_distance(double *r1, double *r2)
{
    return sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]) + (r1[2] - r2[2]) * (r1[2] - r2[2]));
}

// TODO children must be compacted at the start of array
// TODO remove id_in_warp, call bidimensional kernel
__global__ void get_acceleration(Octtree *tree, Entities *ents, int ents_sz, size_t dt)
{
    uint threads_sz = gridDim.x * blockDim.x;
    uint tid, id, stack_pointer;
    // 1024*48 is the total shared memory, remove space taken from double and from cache_node, divide for 32 (number of warps)
    int id_in_warp, curr_node, depth, *my_stack, warp_id, stack_level;
    double pos_ent[3], vel_ent[3], mass_ent, border, acc_ent[3], ddt;
    __shared__ double pos_node[3 * 32], mass_node[32];

    const int cache_sz = 1024 * 48 - (4 * 32 * (8 / 4)) - 32;
    __shared__ int cache_node[32], stack[cache_sz];
    tid = threadIdx.x;
    id_in_warp = tid % 32;
    id = blockIdx.x * blockDim.x + tid;
    border = tree->max * 2;
    warp_id = tid / 32;
    my_stack = &stack[cache_sz / 32 * warp_id];
    ddt = (double)dt;

    // root is the first node located after the leaves, so each id must be an id of a leaf
    for (int i = tid; i < ents_sz; i += threads_sz)
    {
        // initialize values
        pos_ent[0] = ents->pos[i * 3];
        pos_ent[1] = ents->pos[i * 3 + 1];
        pos_ent[2] = ents->pos[i * 3 + 2];
        vel_ent[0] = ents->vel[i * 3];
        vel_ent[1] = ents->vel[i * 3 + 1];
        vel_ent[2] = ents->vel[i * 3 + 2];
        mass_ent = ents->mass[i];
        acc_ent[0] = 0;
        acc_ent[1] = 0;
        acc_ent[2] = 0;
        if (id_in_warp == 0)
        {
            stack_level = 0;
        }

        depth = 0;
        stack_pointer = 0;
        curr_node = tree->root;

        // start iterate algorithm, when depth<0 the algorithm has ended
        while (depth >= 0)
        {
            for (int j = 0; j < 8; j++)
            {
                // save in shared memory the node
                // TODO load only id and break if <0
                if (id_in_warp == 0)
                {
                    cache_node[warp_id] = tree->children[curr_node * 8 + j];
                    pos_node[warp_id * 3] = tree->center[cache_node[warp_id] * 3];
                    pos_node[warp_id * 3 + 1] = tree->center[cache_node[warp_id] * 3 + 1];
                    pos_node[warp_id * 3 + 2] = tree->center[cache_node[warp_id] * 3 + 2];
                    mass_node[warp_id] = tree->mass[cache_node[warp_id]];
                }
                // be sure that every thread of block (but we are interested in thread in the same warp) can see the node
                __threadfence_block();
                if (cache_node[warp_id] >= 0)
                {
                    double distance = get_distance(pos_ent, pos_node);
                    // if the node is a leaf or all thread in warp have their entity far enough from the body
                    if (cache_node[warp_id] < tree->root || __all_sync(0xFFFFFFFF, (border / depth) < THETA))
                    {
                        calculate_acceleration(&pos_node[warp_id * 3], mass_node[warp_id], pos_ent, mass_ent, acc_ent);
                    }
                    else
                    {
                        // save node in stack, when node will be taken from the stack, the programm will check its children
                        depth++;
                        if (id_in_warp == 0)
                        {
                            my_stack[stack_level] = cache_node[warp_id];
                            stack_level++;
                        }
                        __threadfence_block();
                    }
                }
                else
                {
                    depth = max(0, depth - 1);
                }
            }
            depth--;
            stack_level--;
            curr_node = my_stack[stack_level];
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
    default:
        break;
    }
}


// TODO check errors from malloc and kernels
int main(int argc, char *argv[])
{
    // opt_thread is value to optimize
    uint n_ents, *d_tents, opt_thread, opt_block, cache_sz;
    dim3 block, grid;
    int *d_tchildren, *d_tparent, *d_sorted_ents;
    // *_e* memory for entity data
    // *_t* memory for tree data
    // h_le* (host) memory locked for entity data copies
    double *d_epos, *d_evel, *d_emass, *d_tcenter, *d_tmass, *h_lepos, *h_level, *d_reduce1, *d_reduce2;
    Entities h_ents_struct, *d_ents_struct, h_ents_cpy;
    Octtree h_tree_cpy, *d_tree;
    size_t start, end, dt;
    cache_sz = 0;
    if (argc < 6 || argc > 7)
    {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename [cache_sz_MB]\n", argv[0]);
        return 1;
    }

    n_ents = get_entities(argv[1], &h_ents_struct);
    start = strtoul(argv[2], NULL, 10);
    end = strtoul(argv[3], NULL, 10);
    dt = strtoul(argv[4], NULL, 10);
    if (argc == 7)
    {
        cache_sz = strtoul(argv[6], NULL, 10);
    }

    int sz =         // TODO get size of number of node
        opt_thread = // calculate value
        opt_block = (n_ents - 1) / opt_thread + 1;

    // TODO order entities
    // not allocating memory for h_ents_struct
    cudaMallocHost(&h_lepos, sizeof(double) * n_ents * 3);
    cudaMallocHost(&h_level, sizeof(double) * n_ents * 3);
    cudaMemcpy(h_lepos, h_ents_struct.pos, sizeof(double) * n_ents * 3, cudaMemcpyHostToHost);
    cudaMemcpy(h_level, h_ents_struct.vel, sizeof(double) * n_ents * 3, cudaMemcpyHostToHost);

    cudaMalloc(&d_ents_struct, sizeof(Entities));
    cudaMalloc(&d_tree, sizeof(Octtree));
    cudaMalloc(&d_epos, sizeof(double) * n_ents * 3);
    cudaMalloc(&d_evel, sizeof(double) * n_ents * 3);
    cudaMalloc(&d_emass, sizeof(double) * n_ents);
    cudaMalloc(&d_tcenter, sizeof(double) * sz * 3);
    cudaMalloc(&d_tmass, sizeof(double) * sz);
    cudaMalloc(&d_tents, sizeof(uint) * sz);
    cudaMalloc(&d_tparent, sizeof(uint) * sz);
    cudaMalloc(&d_tchildren, sizeof(uint) * sz * 8);
    cudaMalloc(&d_sorted_ents, sizeof(int)*n_ents);
    if(n_ents>1024*2){
        cudaMalloc(&d_reduce1, sizeof(double)*((n_ents*3-1)/1024+1));
    }else{
        //the final destination of result is the tree attribute max
        //need to be defined also here because first reduction is out of the loop (if loop will be executed, the program will not enter this branch)
        d_reduce1=&d_tree->max;
    }
    if(n_ents>1024*2*1024*2){
        cudaMalloc(&d_reduce2, sizeof(double)*((n_ents*3-1)/1024/1024+1));

    }

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

    cudaMemcpy(d_tree, &h_tree_cpy, sizeof(Octtree), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ents_struct, &h_ents_cpy, sizeof(Entities), cudaMemcpyHostToDevice);
    cudaMemcpy(d_epos, h_lepos, sizeof(double) * n_ents * 3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_evel, h_level, sizeof(double) * n_ents * 3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_emass, h_ents_struct.mass, sizeof(double) * n_ents, cudaMemcpyHostToDevice);

    free(h_ents_struct.pos);
    free(h_ents_struct.vel);
    free(h_ents_struct.mass);
    // TODO recalculate thread and block size
    for (size_t t = start; t < end; t += dt)
    {
        double *d_reduce_in, *d_reduce_out, *temp_swap;
        block.x = 32;
        block.y = 6;
        set_tree<<<1, block>>>(d_tree);
        uint ents = n_ents;
        get_bounding_box<<<(sz-1)/(1024*2)+1,1024>>>(d_epos, sz, d_reduce1);
        cudaDeviceSynchronize();
        d_reduce_in=d_reduce2;
        d_reduce_out=d_reduce1;
        for(int temp_sz=(sz-1)/(1024*2)+1; temp_sz>1; temp_sz=(temp_sz-1)/(1024*2)+1){
            if(temp_sz<=1024*2){
                temp_swap=d_reduce_in;
                d_reduce_in=d_reduce_out;
                d_reduce_out=temp_swap;

            }else{
                d_reduce_in=d_reduce_out;
                //the final destination of result is the tree attribute max
                d_reduce_out=&d_tree->max;
            }

            get_bounding_box<<<(sz-1)/(1024*2)+1,1024>>>(d_reduce_in, sz, d_reduce_out);
            cudaDeviceSynchronize();
        }

        add_ent<<<opt_block, opt_thread>>>(d_tree, d_ents_struct, n_ents);
        cudaDeviceSynchronize();
        set_branch_values<<<opt_block, opt_thread>>>(d_tree);
        cudaDeviceSynchronize();
        order_ents<<<opt_block, opt_thread>>>(d_tree, d_sorted_ents);
        cudaDeviceSynchronize();
        get_acceleration<<<opt_block, opt_thread>>>(d_tree, d_ents_struct, n_ents, dt);
        cudaDeviceSynchronize();
        cudaMemcpy(h_level, d_evel, sizeof(double) * n_ents * 3, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_lepos, d_epos, sizeof(double) * n_ents * 3, cudaMemcpyDeviceToHost);
    }

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
    cudaFree(d_sorted_ents);
    if(n_ents>1024*2){
        cudaFree(d_reduce1);
    }
    if(n_ents>1024*2*1024*2){
        cudaFree(d_reduce2);

    }
}
