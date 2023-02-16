// https://lewiscoleblog.com/barnes-hut

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
    double min;
    double max;
    Octnode *nodes;
} Octtree;

const double BIG_G = 6.67e-11;
const double MIN_NODE_DIM = 10;

//TODO put in a common file
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
    *max = val->x > *max ? val->x : *max;
    *max = val->y > *max ? val->y : *max;
    *max = val->z > *max ? val->z : *max;
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
    int indx;
    int x, y, z;
    z = pos->z >= center->z;
    y = pos->y >= center->y;
    x = pos->x >= center->x;
    indx = z * 4 + y * 2 + x;
    // used to calculate new center
    double bord_div4 = *border_size / 4;
    center->x += x ? bord_div4 : -(bord_div4); // double(x)*2*border_size - border_size
    center->y += y ? bord_div4 : -(bord_div4);
    center->z += z ? bord_div4 : -(bord_div4);
    *border_size /= 2;

    return indx;
}

void double_Octtree(Octtree *tree){
    tree->sz*=2;
    realloc(tree->nodes, tree->sz);

}

void init_node(Octnode *node){
    node->center.x=0;
    node->center.y=0;
    node->center.z=0;
    node->mass=0;
    node->ents=0;
    node->parent=-1;
    for (uint i = 0; i < 8; i++)
    {
        node->children[i]=-1;
    }
    
}

void add_ent(Octtree *tree, Entity *ent, int id)
{
    //allocated is used as a boolean
    int allocated, node_indx, indx;
    Octnode *node;
    double border_size;
    RVec3 volume_center;
    // set init value
    allocated = 0;
    // keep last visited node index
    node_indx = tree->root;
    node = &tree->nodes[node_indx];
    border_size = tree->max - tree->min;

    // set center of whole volume
    volume_center.x = border_size / 2;
    volume_center.y = border_size / 2;
    volume_center.z = border_size / 2;

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
                //other is the other leaf
                RVec3 other_center=volume_center;
                double other_border=border_size;
                int other = node->children[indx];
                int other_indx;
                do{
                    //double space if tree is full
                    if(tree->firstfree==tree->sz){
                        double_Octtree(tree);
                    }
                    
                    //take first free location and set the parent of the new branch
                    init_node(&tree->nodes[tree->firstfree]);
                    tree->nodes[tree->firstfree].parent=node_indx;
                    //set the new branch as child
                    node->children[indx]=tree->firstfree;
                    
                    //get leaves position in the new branch
                    indx = get_indx_loc(&ent->pos, &volume_center, &border_size);
                    //the center of the leaf is the position of the entity associated
                    other_indx = get_indx_loc(&tree->nodes[other].center, &other_center, &other_border);

                    //use the new branch as the current one
                    node_indx=tree->firstfree;
                    node=&tree->nodes[node_indx];
                    //update first free location
                    tree->firstfree++;
                    
                //if the leaves will be in different position exit the loop
                }while(indx == other_indx);

                //set new parent in the leaves values
                tree->nodes[other].parent=node_indx;
                tree->nodes[id].parent=node_indx;

                //set the leaves as branch children
                node->children[indx]=id;
                node->children[other_indx]=other;
                

                allocated=1;
            }
            else
            {
                //change current node to the child node
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


//calculate number of entities, mass and center of each branch
void set_branch_values(Octtree *tree)
{
    //root is the first node after the leaves with an entity, so from location 0 to root-1 there are all leaves from which
    //will be calculated the branch values

    //avoid to write everytime tree->nodes
    int parent_indx, child_indx;
    double new_mass;
    // Octnode *nodes=tree->nodes;
    Octnode *parent, *child;
    for (int i = 0; i < tree->root; i++)
    {
        child_indx=i;
        while(tree->nodes[child_indx].parent>-1){
            child=&tree->nodes[child_indx];
            parent_indx=child->parent;
            parent=&tree->nodes[parent_indx];

            parent->ents++;
            new_mass=parent->mass + child->mass;
            //calculate center
            parent->center.x=(child->center.x * child->mass / new_mass) + (parent->center.x * parent->mass / new_mass);
            parent->center.y=(child->center.y * child->mass / new_mass) + (parent->center.y * parent->mass / new_mass);
            parent->center.z=(child->center.z * child->mass / new_mass) + (parent->center.z * parent->mass / new_mass);
            
            parent->mass=new_mass;

            child_indx=parent_indx;
        }
    }
    
}

void get_bounding_box(Entity ents[], int ents_sz, double *max, double *min)
{
    default_max_min(&max, &min);
    for (int i = 0; i < ents_sz; i++)
    {
        update_max_min(&max, &min, &ents[i].pos);
    }
}

// create tree struct and add root node
void init_tree(Entity ents[], int ents_sz, Octtree *tree)
{
    get_bounding_box(ents, ents_sz, &tree->max, &tree->min);
    // calculate the minimum quantity of branch required to save ents_sz bodies
    int sz = (ents_sz - 2) / 3 + 1; // = round up (etns_sz-1)/3
    sz *= 2;                        // double the size to leave some space without need to reallocate
    // add the space required for the bodyes
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
    root.parent = NULL;
    root.ents = 0;
    tree->nodes[ents_sz] = root;
}

void create_tree(Entity ents[], int ents_sz, Octtree *tree)
{
    // get bounding box of bodies
    init_tree(ents, ents_sz, tree);
}

//will calculate the bodies position over time
void propagation(Entity ents[], int ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output)
{
    FILE *fpt;
    // create root and set values
    // root=allocate_node(NULL);
    fpt = fopen(output, "w");
    for (size_t t = t_start; t < t_end; t += dt)
    {
    }

    fclose(fpt);
}

int main(int argc, char *argv[])
{
    uint n_ents;
    Entity *ents;
    size_t start, end, dt;
    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
    }
    n_ents = get_entities(argv[1], &ents);
    start = strtoul(argv[2], NULL, 10);
    end = strtoul(argv[3], NULL, 10);
    dt = strtoul(argv[4], NULL, 10);
    propagation(ents, n_ents, start, end, dt, argv[5]);

    free(ents);
}