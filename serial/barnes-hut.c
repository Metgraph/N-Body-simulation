// https://lewiscoleblog.com/barnes-hut
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BILLION 1000000000.0

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

// use int instead of uint for indexs so -1 can be used as a sort of null value
typedef struct
{
    uint ents; // entities in this section (numero figli del nodo - da controllare)
    double mass; //massa dei figli
    RVec3 center;    // mass center media ponderata in base alla massa
    int parent;      // index of parent
    int children[8]; // indexs of children
} Octnode;

// use int for same reason commented above Octnode and to avoid to get number higher than int max value
typedef struct
{
    int sz;        // number of total slot in array
    int firstfree; // first location free
    int root;
    double max; // Root border lenght
    Octnode *nodes;
} Octtree;

//const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;
const double THETA = 0.5; // Theta = 0: senza approssimazione

// TODO put in a common file
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

// get the index of branch where body will be placed.
// Center is the center of volume of branch. It will calculate the center of next branch
// border_size contains the border size of the current branch volume (so the volume is border_size^3)
uint get_indx_loc(RVec3 *pos, RVec3 *center, double *border_size)
{
    int indx;
    int x, y, z;
    double bord_div4 = *border_size / 4;

    z = pos->z >= center->z;
    y = pos->y >= center->y;
    x = pos->x >= center->x;

    indx = z * 4 + y * 2 + x; // Posizione all'interno di uno degli 8 quadranti del cubo

    // Calculate the new center
    center->x += x ? bord_div4 : -(bord_div4);
    center->y += y ? bord_div4 : -(bord_div4);
    center->z += z ? bord_div4 : -(bord_div4);
    *border_size /= 2; // Side's size of inner box

    return indx;
}

void double_Octtree(Octtree *tree)
{
    tree->sz *= 2;
    tree->nodes=realloc(tree->nodes, tree->sz*sizeof(Octnode));
}

// Set empty node
void init_node(Octnode *node)
{
    node->center.x = 0;
    node->center.y = 0;
    node->center.z = 0;
    node->mass = 0;
    node->ents = 0;
    node->parent = -1;
    for (uint i = 0; i < 8; i++)
    {
        node->children[i] = -1;
    }
}

// Add a entity in the tree
// creating all the necessary branches
void add_ent(Octtree *tree, Entity *ent, int id)
{
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

    do
    {
        // center and border_size are updated to the next branch value
        body_pos = get_indx_loc(&ent->pos, &volume_center, &border_size);
        allocated = node->children[body_pos] == -1; // TODO: cambiare con notAllocated
        if (allocated) // Questo caso è se la posizione non è già assegnata non sta coprendo
        {
            tree->nodes[node_indx].children[body_pos] = id;
            tree->nodes[node_indx].ents++;
        }
        else
        {
            // if the location is occupied by a body-leaf
            // [leaf, leaf, leaf, ..., root, branch, branch, ...]
            if (node->children[body_pos] < tree->root)
            {
                // other is the other leaf
                RVec3 other_center = volume_center;
                double other_border = border_size;
                int other = node->children[body_pos];
                int other_indx;

                // Add new body to count in the current node
                tree->nodes[node_indx].ents++;
                do
                {
                    // double space if tree is full
                    if (tree->firstfree >= tree->sz)
                    {

                        printf("Sto allargando la memoria dell'albero");
                        double_Octtree(tree);
                        //update the pointer to new address
                        node = &tree->nodes[node_indx];
                    }

                    // take first free location and set the parent of the new branch
                    int free = tree->firstfree;

                    init_node(&tree->nodes[free]);
                    tree->nodes[free].parent = node_indx;
                    // set the new branch as child
                    node->children[body_pos] = free;
                    tree->nodes[free].ents = 2;

                    // get leaves position in the new branch
                    body_pos = get_indx_loc(&ent->pos, &volume_center, &border_size);
                    // the center of the leaf is the position of the entity associated
                    other_indx = get_indx_loc(&tree->nodes[other].center, &other_center, &other_border);
                   // use the new branch as the current one
                    node_indx = free;
                    node = &tree->nodes[node_indx];
                    // update first free location
                    tree->firstfree++;

                    // if the leaves will be in different position exit the loop
                } while (body_pos == other_indx);

                // set new parent in the leaves values
                tree->nodes[other].parent = node_indx;
                tree->nodes[id].parent = node_indx;

                // set the leaves as branch children
                node->children[body_pos] = id;
                node->children[other_indx] = other;

                allocated = 1;
            }
            else
            {
                // The current node will have one more body in its subtree
                tree->nodes[node_indx].ents++;
                // cross the branch
                node_indx = node->children[body_pos];
                node = &tree->nodes[node_indx];
            }
        }
    } while (!allocated);

    // Verificare: a sinistra ci sono le foglie ma sono in realtà i pianeti
    tree->nodes[id].center = ent->pos;
    tree->nodes[id].mass = ent->mass;
    tree->nodes[id].ents = 1;
    tree->nodes[id].parent = node_indx;
}

// Add the entities in the tree
// The entities are located in the first positions, their position in the tree array is the same position in the ents array.
// The tree is not ready to use, some branch values are not set. Use set_branch_value to complete the tree
void add_ents(Octtree *tree, Entity ents[], int ents_sz)
{
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

void get_bounding_box(Entity ents[], int ents_sz, double *max)
{
    *max = 0.0;

    for (int i = 0; i < ents_sz; i++)
    {
        *max = fabs(ents[i].pos.x) > *max ? fabs(ents[i].pos.x) : *max;
        *max = fabs(ents[i].pos.y) > *max ? fabs(ents[i].pos.y) : *max;
        *max = fabs(ents[i].pos.z) > *max ? fabs(ents[i].pos.z) : *max;
    }
    *max *= 2;
}

// create tree struct and add root node
void init_tree(Entity ents[], int ents_sz, Octtree *tree)
{
    // calculate the minimum quantity of branch required to save ents_sz bodies
    // In pratica è ents_sz * 2/3
    int sz = (ents_sz - 2) / 3 + 1; // = round up (ents_sz-1)/3
    sz *= 2;                        // double the size to leave some space without need to reallocate

    // add the space required for the bodies
    sz += ents_sz;
    tree->firstfree = ents_sz + 1;
    tree->sz = sz;
    tree->root = ents_sz;
    tree->nodes = malloc(sz * sizeof(Octnode));
    // TODO: sostituire coni init node
    Octnode root;
    root.center.x = 0;
    root.center.y = 0;
    root.center.z = 0;
    root.mass = 0;
    root.parent = -1;
    root.ents = 0;
    for(int i=0; i<8; i++){
        root.children[i]=-1;
    }
    tree->nodes[ents_sz] = root;
}

void create_tree(Entity ents[], int ents_sz, Octtree *tree)
{
    // get bounding box of bodies
    init_tree(ents, ents_sz, tree);
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

//calculate body acceleration - TODO: trasformare in iterativo
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
        for (int i = 0; i < 8; i++)
        {
            int indx = node->children[i];
            //if (indx > -1 && indx != id)
            if (indx > -1)
            {
                get_acceleration_rec(tree, indx, id, acc, border / 2);
            }
        }
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

// will calculate the bodies position over time
void propagation(Entity ents[], int ents_sz, int n_steps, float dt, const char *output)
{
    FILE *fpt;
    Octtree tree;
    create_tree(ents, ents_sz, &tree);
    fpt = fopen(output, "w");
    RVec3 *acc;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    for (size_t i = 0; i < ents_sz; i++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i,
                ents[i].pos.x, ents[i].pos.y,
                ents[i].pos.z, ents[i].mass);
    }

    // Initialize system
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
        init_node(&tree.nodes[tree.root]);
        tree.firstfree=tree.root+1;
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

int main(int argc, char *argv[])
{
    uint n_ents;
    Entity *ents;
    float start, end, dt;
    int n_steps;
    struct timespec s, e;

    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
        return 1;
    }

    n_ents = get_entities(argv[1], &ents);

    start = strtof(argv[2], NULL);
    end = strtof(argv[3], NULL);
    dt = strtof(argv[4], NULL);

    n_steps = (end - start) / dt;

    printf("Start: %f, end: %f, delta time: %f, time steps: %d, ents: %d, G: "
           "%lf\n",
           start, end, dt, n_steps, n_ents, BIG_G);

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
