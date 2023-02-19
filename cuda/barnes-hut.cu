#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

typedef unsigned int uint;

typedef struct
{
    double x;
    double y;
    double z;
} RVec3;

typedef struct
{
    RVec3 *pos;
    RVec3 *vel;
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
    //Octnodes
    uint *ents;
    double *mass;
    RVec3 *center;
    int *parent;
    int *children; //must have 8 slots for each node

} Octtree;