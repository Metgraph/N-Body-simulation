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
    double *pos; //allocation size must be three times the size of mazz
    double *vel; //allocation size must be three times the size of mazz
    double *mass;
} Entities;

// use int instead of uint for indexs so -1 can be used as a sort of null value
/*
typedef struct
{
    uint ents; // entities in this section
    double mass;
    RVec3 center;    // mass center
    int parent;      // index of parent
    int children[8]; // indexs of children
} Octnode;

*/

typedef struct
{
    int sz;        // number of total slot in array
    int firstfree; // first location free
    int root;
    double max;
    // Octnodes
    uint *ents;
    double *mass;
    RVec3 *center;
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
