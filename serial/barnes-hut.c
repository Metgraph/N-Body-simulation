// https://lewiscoleblog.com/barnes-hut
//TODO check if algorithm is okay, results are a bit different from exhaustive algorithm, but nothing strange

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

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
    uint ents; // entities in this section
    double mass;
    RVec3 center;    // mass center
    int parent;      // index of parent
    int children[8]; // indexs of children
} Octnode;

// use int for same reason commented above Octnode and to avoid to get number higher than int max value
typedef struct
{
    int sz;        // number of total slot in array
    int firstfree; // first location free
    int root;
    double min; //TODO remove min or update code and use min in the correct way
    double max;
    Octnode *nodes;
} Octtree;

const double BIG_G = 6.67e-11;
const double THETA = 0.5; // 1;

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
    double border=(tree->max*2)/(double)temp;
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

double get_distance(RVec3 *r1, RVec3 *r2)
{
    return sqrt(pow(r1->x - r2->x, 2) + pow(r1->y - r2->y, 2) + pow(r1->z - r2->z, 2));
}

Octnode *allocate_node(int parent)
{
    Octnode *ret;
    // alloc and initialize all to 0 (and then NULL)
    ret = (Octnode *)calloc(1, sizeof(Octnode));
    ret->parent = parent;
    return ret;
}

void update_max_min(double *max, double *min, RVec3 *val)
{
    // update max
    *max = fabs(val->x) > *max ? fabs(val->x) : *max;
    *max = fabs(val->y) > *max ? fabs(val->y) : *max;
    *max = fabs(val->z) > *max ? fabs(val->z) : *max;
    // update min
    *min = val->x < *min ? val->x : *min;
    *min = val->y < *min ? val->y : *min;
    *min = val->z < *min ? val->z : *min;
}

void default_max_min(double *max, double *min)
{
    *max = DBL_MIN;
    *min = DBL_MAX;
}

// get the index of branch where body will be placed.
// Center is the center of volume of branch. It will calculate the center of next branch
// border_size contains the border size of the current branch volume (so the volume is border_size^3)
uint get_indx_loc(RVec3 *pos, RVec3 *center, double *border_size)
{
    int index=0;
    double border4=*border_size/4;
    if(pos->x<center->x){
        center->x-= border4;
    }else{
        center->x+=border4;
        index+=1;
    }
    if(pos->y<center->y){
        center->y-= border4;
    }else{
        center->y+=border4;
        index+=1*2;
    }
    if(pos->z<center->z){
        center->z-= border4;
    }else{
        center->z+=border4;
        index+=1*4;
    }
    *border_size/=2;
    
    return index;
    // int indx;
    // int x, y, z;
    // z = pos->z >= center->z;
    // y = pos->y >= center->y;
    // x = pos->x >= center->x;
    // indx = z * 4 + y * 2 + x;
    // // used to calculate new center
    // double bord_div4 = *border_size / 4;
    // center->x += x ? bord_div4 : -(bord_div4); // double(x)*2*border_size - border_size
    // center->y += y ? bord_div4 : -(bord_div4);
    // center->z += z ? bord_div4 : -(bord_div4);
    // *border_size /= 2;

    // return indx;
}

void double_Octtree(Octtree *tree)
{
    tree->sz *= 2;
    tree->nodes=realloc(tree->nodes, tree->sz*sizeof(Octnode));
}

double border_tree(Octtree *tree){
    return tree->max*2;
}

//set value for a empty node
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

//add a entity in the tree
//it's create all the needed branch
void add_ent(Octtree *tree, Entity *ent, int id)
{
    // allocated is used as a boolean
    int allocated, node_indx, indx;
    Octnode *node;
    double border_size;
    RVec3 volume_center;
    // set init value
    allocated = 0;
    // keep last visited node index
    node_indx = tree->root;
    node = &tree->nodes[node_indx];
    border_size = border_tree(tree);

    // set center of whole volume
    //TODO wrong formula
    volume_center.x = 0;
    volume_center.y = 0;
    volume_center.z = 0;

    do
    {
        // center and border_size are updated to the next branch value
        indx = get_indx_loc(&ent->pos, &volume_center, &border_size);
        allocated = node->children[indx] == -1;
        if (allocated)
        {
            tree->nodes[node_indx].children[indx] = id;
        }
        else
        {
            // if the location is occupied by a leaf (an entity)
            if (node->children[indx] < tree->root)
            {
                // other is the other leaf
                RVec3 other_center = volume_center;
                double other_border = border_size;
                int other = node->children[indx];
                int other_indx;
                do
                {
                    // double space if tree is full
                    if (tree->firstfree >= tree->sz)
                    {
                        
                        // printf("ID: %d, other: %d, border: %lf, other_border: %lf, distance: %lf\n", id, other, border_size, other_border, distance);
                        double_Octtree(tree);
                        //update the pointer to new address
                        node = &tree->nodes[node_indx];
                    }

                    // take first free location and set the parent of the new branch
                    init_node(&tree->nodes[tree->firstfree]);
                    tree->nodes[tree->firstfree].parent = node_indx;
                    // set the new branch as child
                    node->children[indx] = tree->firstfree;

                    // get leaves position in the new branch
                    indx = get_indx_loc(&ent->pos, &volume_center, &border_size);
                    // the center of the leaf is the position of the entity associated
                    other_indx = get_indx_loc(&tree->nodes[other].center, &other_center, &other_border);
                    // double distance = get_distance(&tree->nodes[other].center, &ent->pos);
                    // printf("ID: %d, other: %d, border: %lf, other_border: %lf, distance: %lf, indx: %d, other_indx: %d\n",id, other, border_size, other_border, distance, indx, other_indx);
                    // printf("centerX: %lf, centerY: %lf, centerZ: %lf\notherX: %lf, otherY: %lf, otherZ: %lf\n", volume_center.x, volume_center.y, volume_center.z, other_center.x, other_center.y, other_center.z);
                    // use the new branch as the current one
                    node_indx = tree->firstfree;
                    node = &tree->nodes[node_indx];
                    // update first free location
                    tree->firstfree++;

                    // if the leaves will be in different position exit the loop
                } while (indx == other_indx);

                // set new parent in the leaves values
                tree->nodes[other].parent = node_indx;
                tree->nodes[id].parent = node_indx;

                // set the leaves as branch children
                node->children[indx] = id;
                node->children[other_indx] = other;

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

// calculate number of entities, mass and center of each branch
void set_branch_values(Octtree *tree)
{
    // root is the first node after the leaves with an entity, so from location 0 to root-1 there are all leaves from which
    // will be calculated the branch values

    // avoid to write everytime tree->nodes
    int parent_indx, child_indx;
    double new_mass;
    // Octnode *nodes=tree->nodes;
    Octnode *parent, *child;
    for (int i = 0; i < tree->root; i++)
    {
        child_indx = i;
        while (tree->nodes[child_indx].parent > -1)
        {
            child = &tree->nodes[child_indx];
            parent_indx = child->parent;
            parent = &tree->nodes[parent_indx];

            parent->ents++;
            new_mass = parent->mass + child->mass;
            // calculate center
            //maybe is better check if there are other entities, so to avoid spread of error
            parent->center.x = (child->center.x * child->mass / new_mass) + (parent->center.x * parent->mass / new_mass);
            parent->center.y = (child->center.y * child->mass / new_mass) + (parent->center.y * parent->mass / new_mass);
            parent->center.z = (child->center.z * child->mass / new_mass) + (parent->center.z * parent->mass / new_mass);

            parent->mass = new_mass;

            child_indx = parent_indx;
        }
    }
}

void get_bounding_box(Entity ents[], int ents_sz, double *max, double *min)
{
    default_max_min(max, min);
    for (int i = 0; i < ents_sz; i++)
    {
        update_max_min(max, min, &ents[i].pos);
    }
}

// create tree struct and add root node
void init_tree(Entity ents[], int ents_sz, Octtree *tree)
{
    // get_bounding_box(ents, ents_sz, &tree->max, &tree->min);
    // calculate the minimum quantity of branch required to save ents_sz bodies
    int sz = (ents_sz - 2) / 3 + 1; // = round up (ents_sz-1)/3
    sz *= 2;                        // double the size to leave some space without need to reallocate
    // add the space required for the bodies
    sz += ents_sz;
    tree->firstfree = ents_sz + 1;
    tree->sz = sz;
    tree->root = ents_sz;
    tree->nodes = malloc(sz * sizeof(Octnode));

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

    r_vector.x = ent->center.x - node->center.x;
    r_vector.y = ent->center.y - node->center.y;
    r_vector.z = ent->center.z - node->center.z;

    double r_mag = sqrt(r_vector.x * r_vector.x + r_vector.y * r_vector.y + r_vector.z * r_vector.z);

    double acceleration = -1.0 * BIG_G * (node->mass) / pow(r_mag, 2.0);

    RVec3 r_unit_vector = {r_vector.x / r_mag, r_vector.y / r_mag, r_vector.z / r_mag};

    acc->x += acceleration * r_unit_vector.x;
    acc->y += acceleration * r_unit_vector.y;
    acc->z += acceleration * r_unit_vector.z;
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
        for (int i = 0; i < 8; i++)
        {
            int indx = node->children[i];
            if (indx > -1 && indx != id)
            {
                get_acceleration_rec(tree, indx, id, acc, border / 2);
            }
        }
    }
}

//call recursive function get_acceleration_rec
void get_acceleration(Octtree *tree, int id, RVec3 *acc)
{
    // RVec3 acc = {0, 0, 0};
    acc->x = 0;
    acc->y = 0;
    acc->z = 0;
    get_acceleration_rec(tree, tree->root, id, acc, border_tree(tree));
}

void calculate_propagation(Entity ents[], int ents_sz, Octtree *tree, size_t dt, FILE *fpt)
{

    RVec3 acc;
    // calculate new velocity
    for (int i = 0; i < ents_sz; i++)
    {
        
        get_acceleration(tree, i, &acc);
        ents[i].vel.x+=acc.x*dt;
        ents[i].vel.y+=acc.y*dt;
        ents[i].vel.z+=acc.z*dt;
    }

    // calculate new position
    for (int i = 0; i < ents_sz; i++)
    {
        ents[i].pos.x += ents[i].vel.x * dt;
        ents[i].pos.y += ents[i].vel.y * dt;
        ents[i].pos.z += ents[i].vel.z * dt;
        fprintf(fpt, "%d,%lf,%lf,%lf,%lf,%lf,%lf \n", i, ents[i].pos.x,
                ents[i].pos.y, ents[i].pos.z, ents[i].vel.x, ents[i].vel.y,
                ents[i].vel.z);
    }
}

// will calculate the bodies position over time
void propagation(Entity ents[], int ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output)
{
    FILE *fpt;
    Octtree tree;
    create_tree(ents, ents_sz, &tree);
    fpt = fopen(output, "w");
    for (size_t t = t_start; t < t_end; t += dt)
    {
        // fprintf(fpt, "time: %lu\n", t);
        get_bounding_box(ents, ents_sz, &tree.max, &tree.min);
        add_ents(&tree, ents, ents_sz);
        set_branch_values(&tree);
        // printf("time: %lu\n", t);
        // print_tree(&tree);
        // printf("---------------------------------\n");
        calculate_propagation(ents, ents_sz, &tree, dt, fpt);
        //to reset the tree just reset the root and reset firstfree
        init_node(&tree.nodes[tree.root]);
        tree.firstfree=tree.root+1;

    }

    fclose(fpt);
    free(tree.nodes);
}

int main(int argc, char *argv[])
{
    uint n_ents;
    Entity *ents;
    size_t start, end, dt;
    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
        return 1;
    }

    n_ents = get_entities(argv[1], &ents);
    start = strtoul(argv[2], NULL, 10);
    end = strtoul(argv[3], NULL, 10);
    dt = strtoul(argv[4], NULL, 10);
    propagation(ents, n_ents, start, end, dt, argv[5]);

    // n_ents=get_entities("./tests/sun_earth.csv", &ents);
    // start=0;
    // end=5000;
    // dt=1;
    // propagation(ents, n_ents, start, end, dt, "./test/output.csv");

    free(ents);
    return 0;
}