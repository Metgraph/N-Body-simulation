#include "stdlib.h"
#include "stdio.h"

typedef double real;

typedef struct 
{
    real x;
    real y;
    real z;
} RVec3;

typedef struct
{
    RVec3 pos;
    RVec3 vec;
    double mass;
} Entity;

