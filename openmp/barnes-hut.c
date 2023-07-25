#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define MUTEX

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
    RVec3 pos;
    RVec3 vel;
    double mass;
} Entity;

typedef struct
{
    int ents; // how many ents are in the node
    double mass;
    RVec3 center;
    int children[8]; // must be an array of size 8
#ifdef MUTEX
    omp_lock_t writelocks[8];
#endif
    int parent;
    int depth;

} Octnode;

typedef struct
{
    int firstfree;
    int root;
    double max;
    Octnode *nodes;
    int sz;
    omp_lock_t reallocnodes;
} Octtree;

const double BIG_G = 6.67e-11;
const double THETA = 0.5; // 1;

uint get_entities(char filename[], Entity **ents)
{
    Entity e_buff;
    int status;
    uint ret_size;
    uint size;
    Entity *ret;
    FILE *f = fopen(filename, "r");

    // Check if file has been open correctly, if not return NULL
    if (!f)
    {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }

    // TODO Check for error in allocation
    ret = (Entity *)malloc(1 * sizeof(Entity));

    size = 0;
    ret_size = 1;
    // fscanf return the number of input items successfully matched and assigned
    while ((status =
                fscanf(f, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &e_buff.pos.x,
                       &e_buff.pos.y, &e_buff.pos.z, &e_buff.vel.x,
                       &e_buff.vel.y, &e_buff.vel.z, &e_buff.mass)) == 7)
    {
        size++;
        if (ret_size < size)
        {
            // TODO Check for error in allocation
            ret_size *= 2;
            ret = (Entity *)realloc((void *)ret, ret_size * sizeof(Entity));
        }
        // Save value in first free location
        ret[size - 1] = e_buff;
    }

    // check if while ended because the end of the file has been reached
    if (fgetc(f) != EOF)
    {
        fprintf(stderr, "Error reading file '%s': file is not well formed\n",
                filename);
        fclose(f);
        return 0;
    }

    *ents = ret;
    fclose(f);
    return size;
}

double border_tree(Octtree *tree) { return tree->max * 2; }

void copy_RVec3(RVec3 *dest, RVec3 *src)
{
    dest->x = src->x;
    dest->y = src->y;
    dest->z = src->z;
}

int get_indx_loc(RVec3 *pos, RVec3 *center, double *border_size)
{
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
    center->x +=
        x ? bord_div4 : -(bord_div4); // double(x)*2*border_size - border_size
    center->y += y ? bord_div4 : -(bord_div4);
    center->z += z ? bord_div4 : -(bord_div4);
    *border_size /= 2;

    return indx;
}

int simulate_compare_and_swap(int *ptr, int expected, int new_val)
{
    // int old_val;
    // #pragma omp atomic capture
    // {
    //     old_val = *ptr;
    //     if (old_val == expected) {
    //         *ptr = new_val;
    //     }
    // }
    // return old_val;
    // int s_cap;

    // here we capture the shared variable and also update it if p is larger
    int ret;
#pragma omp atomic compare capture
    {
        ret = *ptr;
        // s_cap = s;
        if (*ptr == expected)
        {
            *ptr = new_val;
        }
    }
    return ret;
}

//TODO create a version with single node allocation
int init_node(Octtree *tree, int depth)
{
    // omp_set_lock(&tree->reallocnodes);
    int new_node;
    // #pragma omp atomic capture
    #pragma omp critical
    {
        new_node = tree->firstfree++;
        if (tree->sz <= tree->firstfree)
        {
            tree->sz *= 2;
            tree->nodes = realloc(tree->nodes, tree->sz*sizeof(Octnode));
        }
    }
    // omp_unset_lock(&tree->reallocnodes);
    tree->nodes[new_node].center.x = 0;
    tree->nodes[new_node].center.y = 0;
    tree->nodes[new_node].center.z = 0;
    tree->nodes[new_node].mass = 0;
    tree->nodes[new_node].ents = 0;
    tree->nodes[new_node].parent = -1;
    tree->nodes[new_node].depth = depth;
    for (uint i = 0; i < 8; i++)
    {
        tree->nodes[new_node].children[i] = -1;
    }
    #ifdef MUTEX
    for (uint i = 0; i < 8; i++)
    {
        omp_init_lock(&tree->nodes[new_node].writelocks[i]);
    }
    #endif
    return new_node;
    //need to move here because must be sure that memory will not be moved during creation
    // omp_unset_lock(&tree->reallocnodes);
}

//CAS is not tested and for sure has some problems. i left it just in case i will want to recover it later
void add_ent(Octtree *tree, Entity *ent, int id)
{
    double border_size;
    int allocated, node_indx, ent_loc, child_indx, child_val, other, curr_depth;
    // Octnode *node;
    RVec3 volume_center = {0, 0, 0};
    RVec3 pos_ent;
    copy_RVec3(&pos_ent, &ent->pos);
    allocated = 0;
    node_indx = tree->root;
    // node = &tree->nodes[node_indx];
    border_size = border_tree(tree);
    curr_depth = 0;

    while (!allocated)
    {
        child_indx = get_indx_loc(&pos_ent, &volume_center, &border_size);
#ifdef CAS
        do
        {
            child_val = node->children[chtree->root, 
            {
                if (simulate_compare_and_swap(&node->children[child_indx], child_val, LOCKED))
                {

#endif
#ifdef MUTEX
                    omp_set_lock(&tree->nodes[node_indx].writelocks[child_indx]);
                    child_val = tree->nodes[node_indx].children[child_indx];
#endif
                    // child_indx=node->children[ent_loc];
                    curr_depth++;
                    // if nothing is located, allocate the leaf
                    if (child_val == -1)
                    {
#pragma omp atomic write // used to resolve false sharing
                        tree->nodes[node_indx].children[child_indx] = id;
                        omp_unset_lock(&tree->nodes[node_indx].writelocks[child_indx]);
                        allocated = 1;
                        // if there is already a leaf
                    }
                    else if (child_val < tree->root)
                    {
                        allocated = 0;
                        other = child_val;
                        int other_loc;
                        RVec3 other_center;
                        copy_RVec3(&other_center, &volume_center);
                        double other_border = border_size;
                        int new_node;
                        while (!allocated)
                        {
// #pragma omp atomic capture
                            // new_node = tree->firstfree++; // new_node has the old firsfree value
                            // TODO create program that manage array
                            new_node=init_node(tree, curr_depth);

                            tree->nodes[new_node].parent = node_indx;
                            curr_depth++;

                            ent_loc = get_indx_loc(&pos_ent, &volume_center, &border_size);
                            other_loc = get_indx_loc(&tree->nodes[other].center, &other_center, &other_border);

                            allocated = other_loc != ent_loc;

#ifdef CAS
                            #pragma omp atomic write
                            tree->nodes[new_node].children[ent_loc] = LOCKED;
                            if (allocated)
                            {
                                #pragma omp atomic write
                                tree->nodes[new_node].children[other_loc] = LOCKED;
                            }
                            #pragma omp atomic write
                            node->children[child_indx] = new_node;
// #pragma omp flush(node->children[child_indx], tree->nodes[new_node].children[ent_loc], tree->nodes[new_node].children[other_loc])

                            node = &tree->nodes[new_node];
#endif
#ifdef MUTEX
                            omp_set_lock(&tree->nodes[new_node].writelocks[ent_loc]);
                            if(allocated){
                                omp_set_lock(&tree->nodes[new_node].writelocks[other_loc]);
                            }
                            tree->nodes[node_indx].children[child_indx]=new_node;
                            omp_unset_lock(&tree->nodes[node_indx].writelocks[child_indx]);
#endif
                            child_indx = ent_loc;
                            node_indx = new_node;
                        }
                        #pragma omp atomic write
                        tree->nodes[node_indx].children[child_indx] = id;
                        omp_unset_lock(&tree->nodes[node_indx].writelocks[child_indx]);
                        #pragma omp atomic write
                        tree->nodes[node_indx].children[other_loc] = other;
                        omp_unset_lock(&tree->nodes[node_indx].writelocks[other_loc]);
                        #pragma omp atomic write
                        tree->nodes[other].parent = node_indx;
// #pragma omp flush(tree->nodes[other].parent, node->children[other_loc], node->children[child_indx])
                    }
                    else if (child_val >= tree->root)
                    {
                        #ifdef CAS
                        #pragma omp atomic write
                        node->children[child_indx] = child_val;

                        #endif
                        #ifdef MUTEX
                            omp_unset_lock(&tree->nodes[node_indx].writelocks[child_indx]);
                        #endif
// #pragma omp flush(node->children[child_indx])
                        node_indx = child_val;
                        // node = &tree->nodes[node_indx];
                    }
                    else
                    {
                        // ERROR
                    }
#ifdef CAS
                }
                else
                {
                    child_val = LOCKED;
                }
            }
        } while (child_val == LOCKED);
#endif
    }
    copy_RVec3(&tree->nodes[id].center, &pos_ent);
    tree->nodes[id].mass = ent->mass;
    tree->nodes[id].ents = 1;
    tree->nodes[id].parent = node_indx;
    tree->nodes[id].depth = curr_depth;
}

void add_ents(Octtree *tree, Entity *ents, uint ents_sz)
{
    #pragma omp parallel for
    for (int i = 0; i < ents_sz; i++)
    {
        add_ent(tree, &ents[i], i);
    }
}

void init_tree(Octtree *tree, int n_ents){
    tree->firstfree=n_ents;
    omp_init_lock(&tree->reallocnodes);
    tree->sz=tree->firstfree*2;
    tree->nodes=malloc(sizeof(Octnode)*tree->sz);
    tree->max=0;
    tree->root=n_ents;

    // Octnode *root = tree->nodes[tree->root];
    // root->center={0,0,0};
}

void default_max(double *max)
{
    *max = 0;
}

void update_max(double *max, RVec3 *val)
{
    // update max
    *max = fabs(val->x) > *max ? fabs(val->x) : *max;
    *max = fabs(val->y) > *max ? fabs(val->y) : *max;
    *max = fabs(val->z) > *max ? fabs(val->z) : *max;
}

//TODO parallelize
void get_bounding_box(Entity ents[], int ents_sz, double *max)
{
    default_max(max);
    for (int i = 0; i < ents_sz; i++)
    {
        update_max(max, &ents[i].pos);
    }
}

void run(Octtree *tree, Entity *ents, uint ents_sz, size_t start, size_t end, size_t dt){
    // for (size_t t = start; t < end; t+=dt)
    // {
        init_node(tree, 0);
        get_bounding_box(ents, ents_sz, &tree->max);
        add_ents(tree, ents, ents_sz);
    // }
    
}

void destroy_tree(Octtree *tree){
    #pragma omp parallel for
    for(int i=0; i<tree->sz; i++){
        for(int j=0; j<8; j++){
            omp_destroy_lock(&tree->nodes[i].writelocks[j]);

        }
    }
    free(tree->nodes);
}

int main(int argc, char *argv[])
{
    Entity *ents;
    uint n_ents;
    size_t start, end, dt;
    // if (argc < 6 || argc > 7)
    // {
    //     fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename [cache_sz_MB]\n", argv[0]);
    //     return 1;
    // }

    // n_ents = get_entities(argv[1], &h_ents_struct);
    // n_ents= get_entities(argv[1], &ents);
    // start = strtoul(argv[2], NULL, 10);
    // end = strtoul(argv[3], NULL, 10);
    // dt = strtoul(argv[4], NULL, 10);

    n_ents =  get_entities("./tests/100_bodies.csv", &ents);
    start=0;
    end=100;
    dt=1;

    Octtree tree;
    init_tree(&tree, n_ents);
    run(&tree, ents, n_ents, start, end, dt);
    free(ents);
    destroy_tree(&tree);

}