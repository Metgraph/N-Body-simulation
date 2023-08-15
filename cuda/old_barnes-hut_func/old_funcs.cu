#include <cuda_runtime_api.h>
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
    uint *shift;

} Octtree;

#define LOCKED -2

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
                        if(child_val == -1){
                            // copy_arrs3(&tree->center[id*3], pos_ent, 0);
                            tree->children[child_indx] = id;
                            atomicAdd(&tree->ents[node_indx], 1);
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
                                atomicAdd(&tree->ents[node_indx], 1);
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
//  TODO write a more optimized function
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
            tree->depth[root] = 0;
        }
    default:
        break;
    }
}

//function doesn't work, it's just an example how to use get_bounding_box kernel
void get_max(){
    // must be allocated
    double *d_reduce1, *d_reduce2;
    //must be allocated and initialized
    double *d_epos;
    int n_ents;
    Octtree *d_tree;
    cudaDeviceProp cuda_prop;
    //must be declared
    cudaError_t cuda_err;

    //START CODE
    int max_threads = cuda_prop.maxThreadsPerBlock;
    double *d_reduce_in, *d_reduce_out, *temp_swap;
    d_reduce_in = d_reduce1;
    d_reduce_out = d_reduce2;
    #ifdef STEP_PRINT
    printf("BOUNDING BOX\n");
    #endif
    cudaMemcpy(d_reduce_in, d_epos, sizeof(double)*n_ents*3, cudaMemcpyDeviceToDevice);
    for(int i=n_ents*3; i>1;i=((i-1)/2048+1)){
        get_bounding_box<<<(i-1)/2048+1, max_threads, max_threads * sizeof(double)>>>(d_reduce_in, i, d_reduce_out);
        cuda_err = cudaDeviceSynchronize();
        // check_error(cuda_err);
        double *temp=d_reduce_out;
        d_reduce_out=d_reduce_in;
        d_reduce_in=temp;
    }
    cudaMemcpy(&d_tree->max, d_reduce_in, sizeof(double), cudaMemcpyDeviceToDevice);
}