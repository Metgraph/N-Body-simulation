// https://lewiscoleblog.com/barnes-hut
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// Check L1 chache line size for correct padding size
// On linux: $ getconf LEVEL1_DCACHE_LINESIZE
typedef struct {
    double val;
    double pad[7];
} pad_double;

typedef struct {
    int ents; // how many ents are in the node
    double mass;
    RVec3 center;
    int children[8]; // must be an array of size 8
    omp_lock_t writelocks[8];
    int parent;
    int depth;
} Octnode;

typedef struct {
    int firstfree;
    int root;
    double max;
    Octnode *nodes;
    int sz;
} Octtree;

// const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;
const double THETA = 0.5; // Theta = 0: senza approssimazione
int thread_count;

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
int get_indx_loc(RVec3 *pos, RVec3 *center, double *border_size) {
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
    center->x += x ? bord_div4 : -(bord_div4);
    center->y += y ? bord_div4 : -(bord_div4);
    center->z += z ? bord_div4 : -(bord_div4);
    *border_size /= 2;

    return indx;
}

/*
 * Initialize new node and update `firstfree` node indicator
 *
 * @param *node     The node to inizialize
 */
int init_node(Octtree *tree) {
    int new_node;
#   pragma omp critical
    {
        new_node = tree->firstfree++;
        if (tree->sz <= tree->firstfree) {
            printf("No more space for new nodes! Exiting.\n");
            exit(1);
        }
    }
    tree->nodes[new_node].center.x = 0;
    tree->nodes[new_node].center.y = 0;
    tree->nodes[new_node].center.z = 0;
    tree->nodes[new_node].mass = 0;
    tree->nodes[new_node].ents = 0;
    tree->nodes[new_node].parent = -1;
    for (uint i = 0; i < 8; i++) {
        tree->nodes[new_node].children[i] = -1;
    }
    return new_node;
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
#           pragma omp atomic write
            tree->nodes[node_indx].children[body_pos] = id;

#           pragma omp atomic update
            tree->nodes[node_indx].ents++;

            tree->nodes[id].center = ent->pos;
            tree->nodes[id].mass = ent->mass;
            tree->nodes[id].ents = 1;
            tree->nodes[id].parent = node_indx;

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
#               pragma omp atomic update
                tree->nodes[node_indx].ents++;

                // When the leaves will be in different position exit the loop
                while (body_pos == other_indx) {

                    // take first free location and set the parent of the new
                    // branch
                    int free;

                    free = init_node(tree);
                    tree->nodes[free].parent = node_indx;

                    // set the new branch as child
#                   pragma omp atomic write
                    node->children[body_pos] = free;

#                   pragma omp atomic write
                    tree->nodes[free].ents = 2;

                    // get leaves position in the new branch
                    int old_body_pos = body_pos;
                    body_pos =
                        get_indx_loc(&ent->pos, &volume_center, &border_size);
                    // the center of the leaf is the position of the entity
                    // associated
                    other_indx = get_indx_loc(&tree->nodes[other].center,
                                              &other_center, &other_border);
                    // use the new branch as the current one
                    node_indx = free;
                    omp_set_lock(&tree->nodes[node_indx].writelocks[body_pos]);
                    if (other_indx != body_pos) {
                        omp_set_lock(
                            &tree->nodes[node_indx].writelocks[other_indx]);
                    }
                    omp_unset_lock(&node->writelocks[old_body_pos]);
                    node = &tree->nodes[node_indx];
                }

                // set new parent in the leaves values
                tree->nodes[other].parent = node_indx;

                tree->nodes[id].parent = node_indx;
                tree->nodes[id].center = ent->pos;
                tree->nodes[id].mass = ent->mass;
                tree->nodes[id].ents = 1;

                // set the leaves as branch children
#               pragma omp atomic write
                node->children[body_pos] = id;

#               pragma omp atomic write
                node->children[other_indx] = other;

                omp_unset_lock(&node->writelocks[body_pos]);
                omp_unset_lock(&node->writelocks[other_indx]);

                allocated = 1;
            } else { // Descend into the tree
                // The current node will have one more body in its subtree
#               pragma omp atomic update
                tree->nodes[node_indx].ents++;

                omp_unset_lock(&node->writelocks[body_pos]);

                // cross the branch
                node_indx = node->children[body_pos];
                node = &tree->nodes[node_indx];
            }
        }
    }
}

/*
 * Add the entities in the tree
 *
 * @param *tree     Tree info struct
 * @param *ents     Array with bodies data
 * @param ents_sz   Total number of bodies
 */
void add_ents(Octtree *tree, Entity *ents, uint ents_sz) {
#   pragma omp for nowait
    for (int i = 0; i < ents_sz; i++) {
        add_ent(tree, &ents[i], i);
    }
#   pragma omp barrier
}

/*
 * Update the mass for the node with child data
 *
 * @param *child    The node with wich calculate
 * @param *center   Main node's center
 * @param *mass     Main node's new mass
 */
void update_center_of_mass(Octnode *child, RVec3 *center, double *mass) {
    double new_mass = *mass + child->mass;

    center->x = (child->center.x * child->mass / new_mass) + (center->x * *mass / new_mass);
    center->y = (child->center.y * child->mass / new_mass) + (center->y * *mass / new_mass);
    center->z = (child->center.z * child->mass / new_mass) + (center->z * *mass / new_mass);

    *mass = new_mass;
}

/*
 * The thread calculate the center of mass for the assigned node.
 * The calculation begin when all children has updated.
 *
 * @param *tree     Tree info struct
 */
void center_of_mass(Octtree *tree) {

#   pragma omp for nowait schedule(static)
    for (int n = tree->firstfree - 1; n >= tree->root; n--) {
        Octnode *my_node = &tree->nodes[n];
        RVec3 l_center = my_node->center;
        double l_new_mass = my_node->mass;
        int j;

        int counter = 0;

        j = my_node->children[counter];
        // Loop until all the children have been checked
        while (counter < 8) {
            if (j == -1) { // If child is empty
                j = my_node->children[++counter];
            } else if (tree->nodes[j].mass != 0) {
                update_center_of_mass(&tree->nodes[j], &l_center, &l_new_mass);
                j = my_node->children[++counter];
            }
        }

        // Update my center and mass
        my_node->center = l_center;
        my_node->mass = l_new_mass;
    }
#   pragma omp barrier
}

/*
 * Find the maximum distance beetween the bodies and the center.
 *
 * @params ents[]       Array with bodies data
 * @params ents_sz      Total number of bodies
 * @params *max_val     Where to save the maximum
 * @params *loc_max     Local maximum for each thread
 */
void get_bounding_box(Entity ents[], int ents_sz, double *max_val,
                      pad_double *loc_max) {
    int id = omp_get_thread_num();
    loc_max[id].val = 0.0;

#   pragma omp for
    for (int i = 0; i < ents_sz; i++) {
        loc_max[id].val = fabs(ents[i].pos.x) > loc_max[id].val
                              ? fabs(ents[i].pos.x)
                              : loc_max[id].val;

        loc_max[id].val = fabs(ents[i].pos.y) > loc_max[id].val
                              ? fabs(ents[i].pos.y)
                              : loc_max[id].val;

        loc_max[id].val = fabs(ents[i].pos.z) > loc_max[id].val
                              ? fabs(ents[i].pos.z)
                              : loc_max[id].val;
    }

    // Find the maximum among those calculated by each thread
#   pragma omp single
    {
        int total_threads = omp_get_num_threads();
        for (int i = 0; i < total_threads; i++) {
            if (loc_max[i].val > *max_val)
                *max_val = loc_max[i].val;
        }
        *max_val *= 2.0;
    }
}

/*
 * Create tree struct and initialize root node
 *
 * @param ents_sz   Total number of bodies
 * @param *tree     Tree info struct
 */
void create_tree(int ents_sz, Octtree *tree) {
    tree->firstfree = ents_sz;
    tree->sz = ents_sz * 5;
    tree->nodes = malloc(sizeof(Octnode) * tree->sz);
    tree->max = 0;
    tree->root = ents_sz;
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
    Octnode *ent = &tree->nodes[id];
    Octnode *node = &tree->nodes[node_indx];
    distance = get_distance(&node->center, &ent->center);

    if (border / distance < THETA || node->ents == 1) {
        calculate_acceleration(ent, node, acc);
    } else {
        for (int i = 0; i < 8; i++) {
            int indx = node->children[i];
            if (indx > -1 && indx != id) {
                get_acceleration_rec(tree, indx, id, acc, border / 2);
            }
        }
    }
}

/*
 * Each thread calculate acceleration for each assigned body
 *
 * @param *tree     Tree info struct
 * @param *acc      Array for the new accelerations
 * @param ents_sz   Total number of bodies
 */
void get_acceleration(Octtree *tree, RVec3 *acc, int ents_sz) {
#   pragma omp for
    for (int i = 0; i < ents_sz; i++) {
        acc[i].x = 0;
        acc[i].y = 0;
        acc[i].z = 0;
        // The body in `i` position start the descent from the root
        get_acceleration_rec(tree, tree->root, i, &acc[i], tree->max);
    }
}

/*
 * Initialize children mutex
 *
 * @param *tree     Tree info struct
 */
void init_mutex(Octtree *tree) {
#   pragma omp for
    for (size_t i = 0; i < tree->sz; i++) {
        for (int j = 0; j < 8; j++) {
            omp_init_lock(&tree->nodes[i].writelocks[j]);
        }
    }
}

/*
 * Destroy children mutex
 *
 * @param *tree     Tree info struct
 */
void destroy_mutex(Octtree *tree) {
#   pragma omp for
    for (size_t i = 0; i < tree->sz; i++) {
        for (int j = 0; j < 8; j++) {
            omp_destroy_lock(&tree->nodes[i].writelocks[j]);
        }
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
    Octtree tree;
    RVec3 *acc;
    pad_double loc_max[thread_count];

    create_tree(ents_sz, &tree);
    init_node(&tree);
    tree.max = 0.0;

    for (int i = 0; i < thread_count; i++)
        loc_max[i].val = 0.0;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

#ifdef RESULTS
    FILE *fpt;
    fpt = fopen(output, "w");
    // Initial positions
    for (size_t i = 0; i < ents_sz; i++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i, ents[i].pos.x, ents[i].pos.y,
                ents[i].pos.z, ents[i].mass);
    }
#endif

#   pragma omp parallel num_threads(thread_count)
    {
        init_mutex(&tree);
        get_bounding_box(ents, ents_sz, &tree.max, loc_max);
        add_ents(&tree, ents, ents_sz);
        center_of_mass(&tree);
        get_acceleration(&tree, acc, ents_sz);

        for (int t = 0; t < n_steps; t++) {
            // 1/2 kick
#           pragma omp for
            for (int i = 0; i < ents_sz; i++) {
                ents[i].vel.x += acc[i].x * dt / 2.0;
                ents[i].vel.y += acc[i].y * dt / 2.0;
                ents[i].vel.z += acc[i].z * dt / 2.0;
            }

            // Move bodies
#           pragma omp for
            for (int i = 0; i < ents_sz; i++) {
                ents[i].pos.x += ents[i].vel.x * dt;
                ents[i].pos.y += ents[i].vel.y * dt;
                ents[i].pos.z += ents[i].vel.z * dt;
            }

#ifdef RESULTS
#           pragma omp single
            for (int i = 0; i < ents_sz; i++) {
                fprintf(fpt, "%d,%lf,%lf,%lf,%lf\n", i, ents[i].pos.x,
                        ents[i].pos.y, ents[i].pos.z, ents[i].mass);
            }
#endif

            // Build new tree
#           pragma omp single
            {
                tree.firstfree = tree.root;
                tree.max = 0.0;
                init_node(&tree);
            }

            get_bounding_box(ents, ents_sz, &tree.max, loc_max);
            add_ents(&tree, ents, ents_sz);
            center_of_mass(&tree);
            get_acceleration(&tree, acc, ents_sz);

            // 2nd 1/2 kick
#           pragma omp for
            for (int i = 0; i < ents_sz; i++) {
                ents[i].vel.x += acc[i].x * dt / 2.0;
                ents[i].vel.y += acc[i].y * dt / 2.0;
                ents[i].vel.z += acc[i].z * dt / 2.0;
            }
        }
        destroy_mutex(&tree);
    } // pragma

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

    if (argc != 7) {
        fprintf(stderr,
                "Usage: %s input_filename start_time end_time delta_time "
                "output_filename THREADS_NUM\n",
                argv[0]);
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
