#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdlib.h>

typedef unsigned int uint;

typedef struct
{
    double x;
    double y;
    double z;
} RVec3;

typedef struct
{
    double *pos; //allocation size must be three times the size of mass
    double *vel; //allocation size must be three times the size of mass
    double *mass;
} Entities;


// use int instead of uint for indexs so -1 can be used as a sort of null value
//used to simplify code, not used in global memory
typedef struct
{
    uint ents; // entities in this section
    double mass;
    double center[3];    // mass center
    int parent;      // index of parent
    int children[8]; // indexs of children
} Octnode;



typedef struct
{
    int sz;        // number of total slot in array
    int firstfree; // first location free
    int root;
    double max;
    // Octnodes
    uint *ents;
    double *mass;
    double *center; //allocation size must be three times the size of mass
    int *parent;
    int *children; // must have 8 slots for each node

} Octtree;

const double BIG_G = 6.67e-11;
const double THETA = 0.5; // 1;
//TODO check from version
const int max_thread_block = 1024;

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
    while ((status =
                fscanf(f, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &pos_buff[0],
                       &pos_buff[1], &pos_buff[2], &vel_buff[0],
                       &vel_buff[1], &vel_buff[2], &mass_buff)) == 7)
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
        pos_ret[size *3 - 3] = pos_buff[0];
        pos_ret[size *3 - 2] = pos_buff[1];
        pos_ret[size *3 - 1] = pos_buff[2];
        vel_ret[size *3 - 3] = vel_buff[0];
        vel_ret[size *3 - 2] = vel_buff[1];
        vel_ret[size *3 - 1] = vel_buff[2];
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

__device__ double border_tree(Octtree *tree){
    return tree->max*2;
}

__device__ void init_node(Octtree *tree, int indx)
{
    tree->center[indx*3] = 0;
    tree->center[indx*3+1] = 0;
    tree->center[indx*3+2] = 0;
    tree->mass[indx]=0;
    tree->ents[indx]=0;
    tree->parent[indx]=-1;
    for (uint i = 0; i < 8; i++)
    {
        tree->children[indx*8+i] = -1;
    }
}

//pos and center must point to the first axis of position
__device__ int get_indx_loc(double *pos, double *center, double *border_size)
{
    int indx;
    int x, y, z;
    double bord_div4 = *border_size / 4;
    z = pos[2] >= center[2];
    y = pos[1] >= center[1];
    x = pos[0] >= center[0];
    indx = z * 4 + y * 2 + x;
    // used to calculate new center
    center[0] += x ? bord_div4 : -(bord_div4); // double(x)*2*border_size - border_size
    center[1] += y ? bord_div4 : -(bord_div4);
    center[2] += z ? bord_div4 : -(bord_div4);
    *border_size /= 2;

    return indx;
}

//add a entity in the tree
//it's create all the needed branch
__global__ void add_ent(Octtree *tree, Entities *ent, int id)
{
    // allocated is used as a boolean
    int allocated, node_indx, indx, child_indx;
    double border_size;
    // Octnode node;
    // set center of whole volume
    double volume_center[3]={0,0,0};
    // set init value
    allocated = 0;

    // keep last visited node index
    node_indx = tree->root;
    // node.ents = tree->ents[node_indx];
    // node.center[0] = tree->ents[node_indx*3];
    // node.center[1] = tree->ents[node_indx*3+1];
    // node.center[2] = tree->ents[node_indx*3+2];
    // for(int i=0; i<8; i++){
    //     node.children[i] = tree->ents[node_indx*8+i];
    // }
    // node.mass = tree->mass[node_indx];
    // node.parent = -1;

    border_size = border_tree(tree);


    do
    {
        // center and border_size are updated to the next branch value
        indx = get_indx_loc(&ent->pos[id*3], volume_center, &border_size);
        child_indx=node_indx*8+indx;
        allocated = tree->children[child_indx] == -1;
        if (allocated)
        {
            tree->children[child_indx] = id;
        }
        else
        {
            // if the location is occupied by a leaf (an entity)
            if (tree->children[child_indx] < tree->root)
            {
                // other is the other leaf
                double other_center[3];
                other_center[0] = volume_center[0];
                other_center[1] = volume_center[1];
                other_center[2] = volume_center[2];
                double other_border = border_size;
                int other = tree->children[child_indx];
                int other_indx;
                do
                {   
                    //we will allocate whole memory
                    // double space if tree is full
                    // if (tree->firstfree >= tree->sz)
                    // {
                        
                    //     // printf("ID: %d, other: %d, border: %lf, other_border: %lf, distance: %lf\n", id, other, border_size, other_border, distance);
                    //     double_Octtree(tree);
                    //     //update the pointer to new address
                    //     node = &tree->nodes[node_indx];
                    // }

                    // take first free location and set the parent of the new branch
                    init_node(tree, tree->firstfree);
                    tree->parent[tree->firstfree] = node_indx;
                    // set the new branch as child
                    tree->children[child_indx] = tree->firstfree;

                    // get leaves position in the new branch
                    indx = get_indx_loc(&ent->pos[id*3], volume_center, &border_size);
                    child_indx=node_indx*8+indx;
                    // the center of the leaf is the position of the entity associated
                    other_indx = get_indx_loc(&tree->center[other*3], other_center, &other_border);
                    // double distance = get_distance(&tree->nodes[other].center, &ent->pos);
                    // printf("ID: %d, other: %d, border: %lf, other_border: %lf, distance: %lf, indx: %d, other_indx: %d\n",id, other, border_size, other_border, distance, indx, other_indx);
                    // printf("centerX: %lf, centerY: %lf, centerZ: %lf\notherX: %lf, otherY: %lf, otherZ: %lf\n", volume_center.x, volume_center.y, volume_center.z, other_center.x, other_center.y, other_center.z);
                    // use the new branch as the current one
                    node_indx = tree->firstfree;
                    // update first free location
                    tree->firstfree++;

                    // if the leaves will be in different position exit the loop
                } while (indx == other_indx);

                // set new parent in the leaves values
                tree->parent[other]= node_indx;
                tree->parent[id] = node_indx;

                // set the leaves as branch children
                tree->children[child_indx] = id;
                tree->children[node_indx*8+other_indx] = other;

                allocated = 1;
            }
            else
            {
                // change current node to the child node
                node_indx = node->children[indx];
                node = &tree->nodes[node_indx];
            }
        }
    } while (!allocated);

    tree->nodes[id].center = ent->pos;
    tree->nodes[id].mass = ent->mass;
    tree->nodes[id].ents = 1;
    tree->nodes[id].parent = node_indx;
}

__device__ void init_node(Octtree *tree, int indx)
{
    tree->center[indx*3]=0;
    tree->center[indx*3+1]=0;
    tree->center[indx*3+2]=0;
    tree->mass[indx]=0;
    tree->ents[indx]=0;
    tree->ents[indx] = -1;
    for (uint i = 0; i < 8; i++)
    {
        tree->children[indx*8+i] = -1;
    }
}

// to check more than max_threads doubles check commented code below
__global__ void get_bounding_box(double *ents, int ents_sz, uint padding)
{
    uint t = threadIdx.x;
    __shared__ double partial_max[max_thread_block];
    int sz = ents_sz - blockDim.x * blockIdx.x;
    sz = sz > blockDim.x ? blockDim.x : sz;
    // copy value from global memory
    if (t < sz)
    {
        partial_max[t] = ents[(t+blockDim.x*blockIdx.x)<<padding];
    }
    else
    {
        // initialize value to 0 (the smallest possible value) to avoid problem with reduce
        partial_max[t] = 0;
    }
    // __syncthreads(); //they are syncend in for loops

    // powerpoint 21 page 42
    for (uint stride = blockDim.x / 2; stride >= 1; stride >>= 1)
    {
        __syncthreads();
        if (t < stride)
        {
            partial_max[t] = partial_max[t] > partial_max[t + stride] ? partial_max[t] : partial_max[t + stride];
        }
    }
    if (t == 0)
    {
        ents[(blockIdx.x * blockDim.x)<<padding] = partial_max[0];
    }
}
    // uint padding=0;
    // //over 3 cycles it will break
    // for(int temp_sz=sz; temp_sz>1; temp_sz=(temp_sz-1)/max_thread_block + 1){
    //     int num_blocks=(temp_sz-1)/max_thread_block + 1;
    //     old_bounding<<<num_blocks, max_thread_block>>>(d_nums, temp_sz, padding);
    //     padding+=10;
    //     cudaDeviceSynchronize();
    // }
    // //result will be at the begin of d_nums

cudaError_t init_tree(Entities *ents, int ents_sz, Octtree *tree)
{
    // calculate the minimum quantity of branch required to save ents_sz bodies
    cudaError_t err;
    int sz = (ents_sz - 2) / 3 + 1; // = round up (ents_sz-1)/3
    sz *= 2;                        // double the size to leave some space without need to reallocate
    // add the space required for the bodies
    sz += ents_sz;
    tree->firstfree = ents_sz + 1;
    tree->sz = sz;
    tree->root = ents_sz;
    // tree->nodes = malloc(sz * sizeof(Octnode));
    err=cudaMalloc(&tree->center, sz * 3 * sizeof(double));
    if(err!=cudaSuccess)
        return err;

    err=cudaMalloc(&tree->mass, sz * sizeof(double));
    if(err!=cudaSuccess)
        return err;

    err=cudaMalloc(&tree->ents, sz * sizeof(uint));
    if(err!=cudaSuccess)
        return err;

    err=cudaMalloc(&tree->children, sz * 8 * sizeof(int));
    if(err!=cudaSuccess)
        return err;

    err=cudaMalloc(&tree->parent, sz * sizeof(int));
    if(err!=cudaSuccess)
        return err;

    //TODO ERROR, tree pointers are in device memory
    // Octnode root;
    // root.center.x = 0;
    tree->center[ents_sz*3]=0;
    // root.center.y = 0;
    tree->center[ents_sz*3+1]=0;
    // root.center.z = 0;
    tree->center[ents_sz*3+2]=0;
    // root.mass = 0;
    tree->mass[ents_sz]=0;
    // root.parent = -1;
    tree->parent[ents_sz]=-1;
    // root.ents = 0;
    tree->ents[ents_sz]=0;
    for(int i=0; i<8; i++){
        tree->children[ents_sz*8+i]=-1;
    }
    // tree->nodes[ents_sz] = root;
    return cudaSuccess;
}

