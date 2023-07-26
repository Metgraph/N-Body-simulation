// https://lewiscoleblog.com/barnes-hut
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define BILLION 1000000000.0

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

// Check L1 chache line size for correct padding size
// On linux: $ getconf LEVEL1_DCACHE_LINESIZE
typedef struct {
    double max;
    double pad[7];
} max_v;

// TODO put in a common file
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
    //useless, but not yet tested with its removal
    omp_lock_t reallocnodes;
} Octtree;

//const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;
const double THETA = 0.5; // Theta = 0: senza approssimazione
int thread_count;


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


void print_tree_rec(Octtree *tree, int id, char *space, uint depth){
    Octnode *node=&tree->nodes[id];
    //how much divide
    uint temp=1<<depth;
    double border=(tree->max)/(double)temp;
    printf("%sid: %d, (x:%lf, y:%lf, z:%lf), border: %lf\n", space, id, node->center.x, node->center.y, node->center.z, border);
    if(node->ents>1){


        int i;
        for(i=depth*4; i<depth*4+4; i++){
            space[i]=' ';

        }
        space[i]='\0';
        for(int i=0; i<8; i++){
            if(node->children[i]>-1){
                print_tree_rec(tree, node->children[i], space, depth+1);
            }
        }
        space[depth*4]='\0';
    }
}

//used for debug
void print_tree(Octtree *tree){
    uint sz_space=4*40;
    char *space=malloc(sz_space*sizeof(char));
    space[0]='\0';
    print_tree_rec(tree, tree->root, space, 0);
    free(space);
}

void print_tree_indented(Octtree *tree, Octnode *node, int tabs, int pos) {
    printf("Body position %d -> ", pos);
    if (node->ents == 1){
        printf("Total ents: %d, parent: %d\n", node->ents, node->parent);
        return;
    }
    printf("Total ents: %d, parent: %d\n", node->ents, node->parent);
    for (int i=0; i<8; i++){
        if (node->children[i] == -1) {
            for (int j=0; j<tabs; j++)
                printf("\t|");
            printf("Child %d is empty\n", i);
        } else {
            for (int j=0; j<tabs; j++)
                printf("\t|");
            printf("Child %d: ", i);
            print_tree_indented(tree, &tree->nodes[node->children[i]], tabs+1, node->children[i]);
        }
    }
}

double get_distance(RVec3 *r1, RVec3 *r2)
{
    return sqrt(pow(r1->x - r2->x, 2) + pow(r1->y - r2->y, 2) + pow(r1->z - r2->z, 2));
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

//now init_node manages firstfree variable and memory reallocation
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
//if one variable between CAS and MUTEX is defined, the other one must be undefined

void double_Octtree(Octtree *tree)
{
    tree->sz *= 2;
    tree->nodes=realloc(tree->nodes, tree->sz*sizeof(Octnode));
}
void add_ent(Octtree *tree, Entity *ent, int id) {
    // allocated is used as a boolean
    int allocated, node_indx, body_pos;
    Octnode *node;
    double border_size;
    RVec3 volume_center;
    // set init value
    allocated = 0;
    // keep last visited node index
    node_indx = tree->root;
    node = &tree->nodes[node_indx];
    border_size = tree->max;

    // set center of whole volume
    volume_center.x = 0;
    volume_center.y = 0;
    volume_center.z = 0;

    while (!allocated) {
        // center and border_size are updated to the next branch value
        body_pos = get_indx_loc(&ent->pos, &volume_center, &border_size);
        omp_set_lock(&node->writelocks[body_pos]);
        if (node->children[body_pos] == -1) {
            #pragma omp atomic write
            tree->nodes[node_indx].children[body_pos] = id;
            #pragma omp atomic update
            tree->nodes[node_indx].ents++;
            allocated = 1;
            omp_unset_lock(&node->writelocks[body_pos]);
        } else {
            // if the location is occupied by a body-leaf
            // [leaf, leaf, leaf, ..., root, branch, branch, ...]
            if (node->children[body_pos] < tree->root) {
                // other is the other leaf
                RVec3 other_center = volume_center;
                double other_border = border_size;
                int other = node->children[body_pos];
                int other_indx = body_pos;

                // Add new body to count in the current node
                #pragma omp atomic update
                tree->nodes[node_indx].ents++;

                // When the leaves will be in different position exit the loop
                while (body_pos == other_indx) {
                    // double up the space if tree is full
                    if (tree->firstfree >= tree->sz) {
                        printf(
                            "WARNING: Sto allargando la memoria dell'albero\n");
                        double_Octtree(tree);
                        // update the pointer to new address
                        node = &tree->nodes[node_indx];
                    }

                    // take first free location and set the parent of the new
                    // branch
                    int free = tree->firstfree;

                    // free = init_node(&tree->nodes[free], 0);
                    free = init_node(tree, 0);
                    tree->nodes[free].parent = node_indx;
                    // set the new branch as child
                    #pragma omp atomic write
                    node->children[body_pos] = free;
                    #pragma omp atomic write //maybe not needed here
                    tree->nodes[free].ents = 2;

                    // get leaves position in the new branch
                    int old_body_pos=body_pos;
                    body_pos =
                        get_indx_loc(&ent->pos, &volume_center, &border_size);
                    // the center of the leaf is the position of the entity
                    // associated
                    other_indx = get_indx_loc(&tree->nodes[other].center,
                                              &other_center, &other_border);
                    // use the new branch as the current one
                    node_indx = free;
                    omp_set_lock(&tree->nodes[node_indx].writelocks[body_pos]);
                    if(other_indx!= body_pos){
                        omp_set_lock(&tree->nodes[node_indx].writelocks[other_indx]);
                    }
                    omp_unset_lock(&node->writelocks[old_body_pos]);
                    node = &tree->nodes[node_indx];
                    // update first free location
                    // tree->firstfree++;
                }

                // set new parent in the leaves values
                tree->nodes[other].parent = node_indx;
                tree->nodes[id].parent = node_indx;

                // set the leaves as branch children
                #pragma omp atomic write
                node->children[body_pos] = id;
                #pragma omp atomic write
                node->children[other_indx] = other;
                omp_unset_lock(&node->writelocks[body_pos]);
                omp_unset_lock(&node->writelocks[other_indx]);

                allocated = 1;
            } else { // Descend into the tree
                // The current node will have one more body in its subtree
                #pragma omp atomic update
                tree->nodes[node_indx].ents++;
                omp_unset_lock(&node->writelocks[body_pos]);
                // cross the branch
                node_indx = node->children[body_pos];
                node = &tree->nodes[node_indx];
            }
        }
    }

    // Verificare: a sinistra ci sono le foglie ma sono in realtà i pianeti
    tree->nodes[id].center = ent->pos;
    tree->nodes[id].mass = ent->mass;
    tree->nodes[id].ents = 1;
    tree->nodes[id].parent = node_indx;
}


void add_ents(Octtree *tree, Entity *ents, uint ents_sz)
{
    #pragma omp parallel for
    for (int i = 0; i < ents_sz; i++)
    {
        add_ent(tree, &ents[i], i);
    }
}

void center_of_mass(Octtree *tree, Octnode *node) {
    Octnode *child;
    double new_mass;

    if (node->ents == 1)
        return;
    for (int n = 0; n < 8; n++) {
        if (node->children[n] != -1)
            center_of_mass(tree, &tree->nodes[node->children[n]]);
    }
    for (int n = 0; n < 8; n++) {
        if (node->children[n] != -1) {
            child = &tree->nodes[node->children[n]];

            new_mass = node->mass + child->mass;

            node->center.x = (child->center.x * child->mass / new_mass) + (node->center.x * node->mass / new_mass);
            node->center.y = (child->center.y * child->mass / new_mass) + (node->center.y * node->mass / new_mass);
            node->center.z = (child->center.z * child->mass / new_mass) + (node->center.z * node->mass / new_mass);

            node->mass = new_mass;
        }
    }
}

void get_bounding_box(Entity ents[], int ents_sz, double *max_val)
{
    max_v loc_max[thread_count];
    *max_val = 0.0;
# pragma omp parallel num_threads(thread_count) shared(loc_max)
    {
        int id = omp_get_thread_num();
        loc_max[id].max = 0.0;
    # pragma omp for
        for (int i = 0; i < ents_sz; i++) {
            // double l_max = 0.0;
            loc_max[id].max = fabs(ents[i].pos.x) > loc_max[id].max ? fabs(ents[i].pos.x) : loc_max[id].max;
            loc_max[id].max = fabs(ents[i].pos.y) > loc_max[id].max ? fabs(ents[i].pos.y) : loc_max[id].max;
            loc_max[id].max = fabs(ents[i].pos.z) > loc_max[id].max ? fabs(ents[i].pos.z) : loc_max[id].max;
        }
    # pragma omp single
        for (int i = 0; i < thread_count; i++) {
            if (loc_max[i].max > *max_val)
                *max_val = loc_max[i].max;
        }

    }
    //printf("Max: %lf\n", *max_val);
    *max_val *= 2.0;
}

void s_get_bounding_box(Entity ents[], int ents_sz, double *max) {
    *max = 0.0;

    for (int i = 0; i < ents_sz; i++) {
        *max = fabs(ents[i].pos.x) > *max ? fabs(ents[i].pos.x) : *max;
        *max = fabs(ents[i].pos.y) > *max ? fabs(ents[i].pos.y) : *max;
        *max = fabs(ents[i].pos.z) > *max ? fabs(ents[i].pos.z) : *max;
    }
    printf("Max: %lf\n", *max);
    *max *= 2;
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

void create_tree(Entity ents[], int ents_sz, Octtree *tree)
{
    // get bounding box of bodies
    init_tree(tree, ents_sz);
}

//calculate calculation caused by another body
void calculate_acceleration(Octnode *ent, Octnode *node, RVec3 *acc)
{

    RVec3 r_vector;

    r_vector.x = node->center.x - ent->center.x;
    r_vector.y = node->center.y - ent->center.y;
    r_vector.z = node->center.z - ent->center.z;

    double inv_r3 = r_vector.x * r_vector.x + r_vector.y * r_vector.y +
                        r_vector.z * r_vector.z + 0.01;
    inv_r3 = pow(inv_r3, -1.5);

    acc->x += BIG_G * r_vector.x * inv_r3 * node->mass;
    acc->y += BIG_G * r_vector.y * inv_r3 * node->mass;
    acc->z += BIG_G * r_vector.z * inv_r3 * node->mass;
}

//calculate body acceleration
void get_acceleration_rec(Octtree *tree, int node_indx, int id, RVec3 *acc, double border)
{
    double distance;
    Octnode *ent = &tree->nodes[id];
    Octnode *node = &tree->nodes[node_indx];
    distance = get_distance(&node->center, &ent->center);

    if (border / distance < THETA || node->ents == 1)
    {
        calculate_acceleration(ent, node, acc);
    }
    else
    {
        /*
        for (int i = 0; i < 8; i++)
        {
            int indx = node->children[i];
            //if (indx > -1 && indx != id)
            if (indx > -1)
            {
                get_acceleration_rec(tree, indx, id, acc, border / 2);
            }
        }
        */
        double half_b = border / 2;
        int *c = node->children;
        if (c[0] != -1 && c[0] != id)
            get_acceleration_rec(tree, c[0], id, acc, half_b);

        if (c[1] != -1 && c[1] != id)
            get_acceleration_rec(tree, c[1], id, acc, half_b);

        if (c[2] != -1 && c[2] != id)
            get_acceleration_rec(tree, c[2], id, acc, half_b);

        if (c[3] != -1 && c[3] != id)
            get_acceleration_rec(tree, c[3], id, acc, half_b);

        if (c[4] != -1 && c[4] != id)
            get_acceleration_rec(tree, c[4], id, acc, half_b);

        if (c[5] != -1 && c[5] != id)
            get_acceleration_rec(tree, c[5], id, acc, half_b);

        if (c[6] != -1 && c[6] != id)
            get_acceleration_rec(tree, c[6], id, acc, half_b);

        if (c[7] != -1 && c[7] != id)
            get_acceleration_rec(tree, c[7], id, acc, half_b);

    }
}

//call recursive function get_acceleration_rec
void get_acceleration(Octtree *tree, RVec3 *acc, int ents_sz)
{
    // RVec3 acc = {0, 0, 0};
    for (int i = 0; i < ents_sz; i++) {
        acc[i].x = 0;
        acc[i].y = 0;
        acc[i].z = 0;
        get_acceleration_rec(tree, tree->root, i, &acc[i], tree->max);
    }
}

void destroy_mutex(Octtree *tree){
    #pragma omp parallel for
    for (size_t i = 0; i < tree->firstfree; i++)
    {
        for(int j=0; j<8; j++){
            omp_destroy_lock(&tree->nodes[i].writelocks[j]);
        }
    }
    
}

// will calculate the bodies position over time
void propagation(Entity ents[], int ents_sz, int n_steps, float dt, const char *output)
{
    FILE *fpt;
    Octtree tree;
    create_tree(ents, ents_sz, &tree);
    init_node(&tree, 0);
    fpt = fopen(output, "w");
    RVec3 *acc;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    // Initial positions
    for (size_t i = 0; i < ents_sz; i++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i,
                ents[i].pos.x, ents[i].pos.y,
                ents[i].pos.z, ents[i].mass);
    }

    tree.max = 0.0;
    get_bounding_box(ents, ents_sz, &tree.max);
    add_ents(&tree, ents, ents_sz);
    center_of_mass(&tree, &tree.nodes[tree.root]);
    get_acceleration(&tree, acc, ents_sz);

    for (int t = 0; t < n_steps; t++)
    {
        // 1/2 kick
        for (int i = 0; i < ents_sz; i++) {
            ents[i].vel.x += acc[i].x * dt / 2.0;
            ents[i].vel.y += acc[i].y * dt / 2.0;
            ents[i].vel.z += acc[i].z * dt / 2.0;
        }

        // Move bodies
        for (int i = 0; i < ents_sz; i++) {
            ents[i].pos.x += ents[i].vel.x * dt;
            ents[i].pos.y += ents[i].vel.y * dt;
            ents[i].pos.z += ents[i].vel.z * dt;

            fprintf(fpt, "%d,%lf,%lf,%lf,%lf\n", i,
                        ents[i].pos.x, ents[i].pos.y,
                        ents[i].pos.z, ents[i].mass);
        }

        // Build new tree
        //TODO think a strategy to reuse mutex without destroying and recreating them
        destroy_mutex(&tree);
        tree.firstfree=tree.root;
        //function takes depth, this was needed in cuda version, but here is not essential. for now put random value
        //after this call firstfree should be equal to root+1
        init_node(&tree, 0);
        get_bounding_box(ents, ents_sz, &tree.max);
        add_ents(&tree, ents, ents_sz);
        center_of_mass(&tree, &tree.nodes[tree.root]);

        get_acceleration(&tree, acc, ents_sz);

        // 2nd 1/2 kick
        for (int i = 0; i < ents_sz; i++) {
            ents[i].vel.x += acc[i].x * dt / 2.0;
            ents[i].vel.y += acc[i].y * dt / 2.0;
            ents[i].vel.z += acc[i].z * dt / 2.0;
        }
    }

    fclose(fpt);
    free(tree.nodes);
    free(acc);
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
    uint n_ents;
    Entity *ents;
    float start, end, dt;
    int n_steps;
    struct timespec s, e;

    if (argc != 7)
    {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename THREADS_NUM\n", argv[0]);
        return 1;
    }

    thread_count = atoi(argv[6]);

    n_ents = get_entities(argv[1], &ents);

    start = strtof(argv[2], NULL);
    end = strtof(argv[3], NULL);
    dt = strtof(argv[4], NULL);

    n_steps = (end - start) / dt;

    printf("Start: %f, end: %f, delta time: %f, time steps: %d, ents: %d, G: "
           "%lf, threads: %d\n",
           start, end, dt, n_steps, n_ents, BIG_G, thread_count);

    clock_gettime(CLOCK_REALTIME, &s);
    propagation(ents, n_ents, n_steps, dt, argv[5]);
    clock_gettime(CLOCK_REALTIME, &e);

    printf("Completed. Output file: %s\n", argv[5]);

    double time_spent =
        (e.tv_sec - s.tv_sec) + (e.tv_nsec - s.tv_nsec) / BILLION;

    printf("Elapsed wall time: %f s\n", time_spent);

    free(ents);
    return 0;
}
