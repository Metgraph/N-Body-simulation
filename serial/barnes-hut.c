#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BILLION 1000000000.0

typedef unsigned int uint;

typedef struct {
    double x;
    double y;
    double z;
} RVec3;

typedef struct {
    RVec3 pos;
    RVec3 vel;
    double mass;
} Entity;

typedef struct {
    uint ents;       // Entities in this section
    double mass;     // Total mass of children
    RVec3 center;    // Mass of node
    int parent;      // Index of parent
    int children[8]; // Indeces of children
} Octnode;

typedef struct {
    int sz;         // Total available space for nodes
    int firstfree;  // First free location
    int root;
    double max;     // Lenght of the root border
    Octnode *nodes;
} Octtree;

// const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;
const double THETA = 0.5;

/**
 * Read file and generate an array of Entity
 *
 * @param   filename    Input filename
 * @param   **ents      Array for bodies information storage
 * @return  Total number of bodies
 */
uint get_entities(char filename[], Entity **ents) {
    Entity e_buff;
    int status;
    uint ret_size;
    uint size;
    Entity *ret;
    FILE *f = fopen(filename, "r");

    // Check if file has been open correctly, if not return NULL
    if (!f) {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }

    ret = (Entity *)malloc(1 * sizeof(Entity));

    size = 0;
    ret_size = 1;
    // fscanf return the number of input items successfully matched and assigned
    while ((status =
                fscanf(f, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &e_buff.pos.x,
                       &e_buff.pos.y, &e_buff.pos.z, &e_buff.vel.x,
                       &e_buff.vel.y, &e_buff.vel.z, &e_buff.mass)) == 7) {
        size++;
        if (ret_size < size) {
            ret_size *= 2;
            ret = (Entity *)realloc((void *)ret, ret_size * sizeof(Entity));
        }
        // Save value in first free location
        ret[size - 1] = e_buff;
    }

    // check if while ended because the end of the file has been reached
    if (fgetc(f) != EOF) {
        fprintf(stderr, "Error reading file '%s': file is not well formed\n",
                filename);
        fclose(f);
        return 0;
    }

    *ents = ret;
    fclose(f);
    return size;
}

/*
 * Calculate Euclidean distance beetween two bodies
 *
 * @param r1    Coordinates of first body
 * @param r2    Coordinates of second body
 * @return      Distance
 */
double get_distance(RVec3 *r1, RVec3 *r2) {
    return sqrt(pow(r1->x - r2->x, 2) + pow(r1->y - r2->y, 2) +
                pow(r1->z - r2->z, 2));
}

/*
 * Calculate in wich branch (child index) where the body will be placed into node.
 *
 * @param *pos          Body positions (x, y, z)
 * @param *center       Center of the node
 * @param *border_size  Total lenght of node's border
 * @return              Body location
 */
uint get_indx_loc(RVec3 *pos, RVec3 *center, double *border_size) {
    int indx;
    int x, y, z;
    double bord_div4 = *border_size / 4;

    z = pos->z >= center->z;
    y = pos->y >= center->y;
    x = pos->x >= center->x;

    // Position inside one of the 8 box of the node
    indx = z * 4 + y * 2 + x;

    // Calculate the new center
    center->x += x ? bord_div4 : -(bord_div4);
    center->y += y ? bord_div4 : -(bord_div4);
    center->z += z ? bord_div4 : -(bord_div4);
    *border_size /= 2; // Side's size of inner box

    return indx;
}

/*
 * Realloc nodes space.
 *
 * @param *tree     Tree info struct
 */
void double_Octtree(Octtree *tree) {
    tree->sz *= 2;
    tree->nodes = realloc(tree->nodes, tree->sz * sizeof(Octnode));
}

/*
 * Initialize empty node
 *
 * @param *node     The node to inizialize
 */
void init_node(Octnode *node) {
    node->center.x = 0;
    node->center.y = 0;
    node->center.z = 0;
    node->mass = 0;
    node->ents = 0;
    node->parent = -1;
    for (uint i = 0; i < 8; i++) {
        node->children[i] = -1;
    }
}

/*
 * Add a entity in the tree creating all the necessary branches
 *
 * @param *tree     Tree info struct
 * @param *ent      The body to insert
 * @param id        Body's index into nodes array
 */
void add_ent(Octtree *tree, Entity *ent, int id) {
    // allocated is used as a boolean
    int allocated, node_indx, body_pos;
    Octnode *node;
    double border_size;
    RVec3 volume_center;

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
        if (node->children[body_pos] == -1) {
            tree->nodes[node_indx].children[body_pos] = id;
            tree->nodes[node_indx].ents++;
            allocated = 1;
        } else {
            // if the location is occupied by a body-leaf
            // e.g. [leaf, leaf, leaf, ..., root, branch, branch, ...]
            if (node->children[body_pos] < tree->root) {
                // other is the other leaf
                RVec3 other_center = volume_center;
                double other_border = border_size;
                int other = node->children[body_pos];
                int other_indx = body_pos;

                // Add new body to count in the current node
                tree->nodes[node_indx].ents++;

                // When the leaves will be in different position exit the loop
                while (body_pos == other_indx) {
                    // double up the space if tree is full
                    if (tree->firstfree >= tree->sz) {
                        double_Octtree(tree);
                        // update the pointer to new address
                        node = &tree->nodes[node_indx];
                    }

                    // take first free location and set the parent of the new
                    // branch
                    int free = tree->firstfree;

                    init_node(&tree->nodes[free]);
                    tree->nodes[free].parent = node_indx;
                    // set the new branch as child
                    node->children[body_pos] = free;
                    tree->nodes[free].ents = 2;

                    // get leaves position in the new branch
                    body_pos =
                        get_indx_loc(&ent->pos, &volume_center, &border_size);
                    // the center of the leaf is the position of the entity
                    // associated
                    other_indx = get_indx_loc(&tree->nodes[other].center,
                                              &other_center, &other_border);
                    // use the new branch as the current one
                    node_indx = free;
                    node = &tree->nodes[node_indx];
                    // update first free location
                    tree->firstfree++;
                }

                // set new parent in the leaves values
                tree->nodes[other].parent = node_indx;
                tree->nodes[id].parent = node_indx;

                // set the leaves as branch children
                node->children[body_pos] = id;
                node->children[other_indx] = other;

                allocated = 1;
            } else { // Descend into the tree
                // The current node will have one more body in its subtree
                tree->nodes[node_indx].ents++;
                // cross the branch
                node_indx = node->children[body_pos];
                node = &tree->nodes[node_indx];
            }
        }
    }

    tree->nodes[id].center = ent->pos;
    tree->nodes[id].mass = ent->mass;
    tree->nodes[id].ents = 1;
    tree->nodes[id].parent = node_indx;
}

/*
 * Add the entities in the tree
 *
 * @param *tree     Tree info struct
 * @param ents[]    Array with bodies data
 * @param ents_sz   Total number of bodies
 */
void add_ents(Octtree *tree, Entity ents[], int ents_sz) {
    for (int i = 0; i < ents_sz; i++) {
        add_ent(tree, &ents[i], i);
    }
}

/*
 * Recursively calculate the center of mass for each node of the tree.
 * The calculation is executed with a post-order tree traversal.
 *
 * @param *tree     Tree info struct
 * @param *node     The node for which the calculation is on
 */
void center_of_mass(Octtree *tree, Octnode *node) {
    Octnode *child;
    double new_mass;

    // The node is a leaf
    if (node->ents == 1)
        return;

    for (int n = 0; n < 8; n++) {
        if (node->children[n] != -1)
            center_of_mass(tree, &tree->nodes[node->children[n]]);
    }

    // On return from each child calculate the mass
    for (int n = 0; n < 8; n++) {
        if (node->children[n] != -1) {
            child = &tree->nodes[node->children[n]];

            new_mass = node->mass + child->mass;

            node->center.x = (child->center.x * child->mass / new_mass) +
                             (node->center.x * node->mass / new_mass);

            node->center.y = (child->center.y * child->mass / new_mass) +
                             (node->center.y * node->mass / new_mass);

            node->center.z = (child->center.z * child->mass / new_mass) +
                             (node->center.z * node->mass / new_mass);

            node->mass = new_mass;
        }
    }
}

/*
 * Find the maximum distance beetween the bodies and the center.
 *
 * @params ents[]   Array with bodies data
 * @params ents_sz  Total number of bodies
 * @params *max     Where to save the maximum
 */
void get_bounding_box(Entity ents[], int ents_sz, double *max) {
    *max = 0.0;

    for (int i = 0; i < ents_sz; i++) {
        *max = fabs(ents[i].pos.x) > *max ? fabs(ents[i].pos.x) : *max;
        *max = fabs(ents[i].pos.y) > *max ? fabs(ents[i].pos.y) : *max;
        *max = fabs(ents[i].pos.z) > *max ? fabs(ents[i].pos.z) : *max;
    }

    // The max is double to consider also the space from the center
    // in the opposite direction
    *max *= 2;
}

/*
 * Create tree struct and initialize root node
 *
 * @param ents_sz   Total number of bodies
 * @param *tree     Tree info struct
 */
void create_tree(int ents_sz, Octtree *tree) {
    tree->firstfree = ents_sz + 1;
    tree->sz = ents_sz * 5;
    tree->root = ents_sz;
    tree->nodes = malloc(tree->sz * sizeof(Octnode));
    Octnode root;
    init_node(&root);
    tree->nodes[ents_sz] = root;
}

/*
 * Calculate acceleration beetween two bodies
 *
 * @param *ent  Main body
 * @param *node Node at a sufficient distance for the calculation
 * @param *acc  New acceleration for ent
 */
void calculate_acceleration(Octnode *ent, Octnode *node, RVec3 *acc) {

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

/*
 * Explore the tree to perform acceleration calculations for the selected body_pos
 *
 * @param *tree         Tree info struct
 * @param *node_indx    Index of the body in examination
 * @param *id           Index of main node on wich perform calculation
 * @param *acc          Acceleration for the main node
 * @param *border       Border of current node
 */
void get_acceleration_rec(Octtree *tree, int node_indx, int id, RVec3 *acc,
                          double border) {
    double distance;
    int indx;
    Octnode *ent;
    Octnode *node;

    // Owner of calculation
    ent = &tree->nodes[id];
    // Node in examination
    node = &tree->nodes[node_indx];

    distance = get_distance(&node->center, &ent->center);

    // If distance beetween body and node is less than THETA
    // or the node is a leaf calculate acceleration
    if (border / distance < THETA || node->ents == 1) {
        calculate_acceleration(ent, node, acc);

    } else { // Descend into tree
        for (int i = 0; i < 8; i++) {
            indx = node->children[i];
            if (indx > -1 && indx != id) {
                get_acceleration_rec(tree, indx, id, acc, border / 2);
            }
        }
    }
}

/*
 * Call recursive function get_acceleration_rec for each body
 *
 * @param *tree     Tree info struct
 * @param *acc      Array for the new accelerations
 * @param ents_sz   Total number of bodies
 */
void get_acceleration(Octtree *tree, RVec3 *acc, int ents_sz) {
    // RVec3 acc = {0, 0, 0};
    for (int i = 0; i < ents_sz; i++) {
        acc[i].x = 0;
        acc[i].y = 0;
        acc[i].z = 0;
        // The body in `i` position start the descent from the root
        get_acceleration_rec(tree, tree->root, i, &acc[i], tree->max);
    }
}

/*
 * Main function for execute the simulation
 *
 * @param ents[]     Array with the bodies information
 * @param ents_sz   Total number of bodies
 * @param n_steps   Total of steps for simulation
 * @param dt        Time interval between one step and the next
 * @param *output    Output file name for simulation results (compile with -DRESULTS)
 */
void propagation(Entity ents[], int ents_sz, int n_steps, float dt,
                 const char *output) {
#ifdef RESULTS
    FILE *fpt;
#endif
    Octtree tree;
    create_tree(ents_sz, &tree);
    RVec3 *acc;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

#ifdef RESULTS
    fpt = fopen(output, "w");
    // Initial positions
    for (size_t i = 0; i < ents_sz; i++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i, ents[i].pos.x, ents[i].pos.y,
                ents[i].pos.z, ents[i].mass);
    }
#endif

    // Initialize the tree
    get_bounding_box(ents, ents_sz, &tree.max);
    add_ents(&tree, ents, ents_sz);
    center_of_mass(&tree, &tree.nodes[tree.root]);
    get_acceleration(&tree, acc, ents_sz);

    for (int t = 0; t < n_steps; t++) {
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

#ifdef RESULTS
            fprintf(fpt, "%d,%lf,%lf,%lf,%lf\n", i, ents[i].pos.x,
                    ents[i].pos.y, ents[i].pos.z, ents[i].mass);
#endif
        }

        // Reset tree
        init_node(&tree.nodes[tree.root]);
        tree.firstfree = tree.root + 1;
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

#ifdef RESULTS
    fclose(fpt);
#endif
    free(tree.nodes);
    free(acc);
}

int main(int argc, char *argv[]) {
    uint n_ents;
    Entity *ents;
    float start, end, dt;
    int n_steps;
    struct timespec s, e;

    if (argc != 6) {
        fprintf(stderr,
                "Usage: %s input_filename start_time end_time delta_time "
                "output_filename\n",
                argv[0]);
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
