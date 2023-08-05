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
    // useless, but not yet tested with its removal
    omp_lock_t reallocnodes;
} Octtree;

// const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;
const double THETA = 0.5; // Theta = 0: senza approssimazione
int thread_count;

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

    // TODO Check for error in allocation
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
            // TODO Check for error in allocation
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

void print_tree_rec(Octtree *tree, int id, char *space, uint depth) {
    Octnode *node = &tree->nodes[id];
    // how much divide
    uint temp = 1 << depth;
    double border = (tree->max) / (double)temp;
    printf("%sid: %d, (x:%lf, y:%lf, z:%lf), mass:%lf, border: %lf\n", space,
           id, node->center.x, node->center.y, node->center.z, node->mass,
           border);
    if (node->ents > 1) {

        int i;
        for (i = depth * 4; i < depth * 4 + 4; i++) {
            space[i] = ' ';
        }
        space[i] = '\0';
        for (int i = 0; i < 8; i++) {
            if (node->children[i] > -1) {
                print_tree_rec(tree, node->children[i], space, depth + 1);
            }
        }
        space[depth * 4] = '\0';
    }
}

// used for debug
void print_tree(Octtree *tree) {
    uint sz_space = 4 * 40;
    char *space = malloc(sz_space * sizeof(char));
    space[0] = '\0';
    print_tree_rec(tree, tree->root, space, 0);
    free(space);
}

void print_tree_indented(Octtree *tree, Octnode *node, int tabs, int pos,
                         int *body_count) {
    printf("Body position %d -> ", pos);
    if (node->ents == 1) {
        int i;
        for (i = 0; i < 500; i++) {
            if (tree->nodes[i].parent == node->parent)
                break;
        }
        printf("Total ents: %d, parent: %d, total mass: %lf, bpos: %d.\n",
               node->ents, node->parent, node->mass, i);
        *body_count += 1;
        return;
    }
    printf("Total ents: %d, parent: %d, total mass: %lf\n", node->ents,
           node->parent, node->mass);
    for (int i = 0; i < 8; i++) {
        if (node->children[i] == -1) {
            for (int j = 0; j < tabs; j++)
                printf("\t|");
            printf("Child %d is empty\n", i);
        } else {
            for (int j = 0; j < tabs; j++)
                printf("\t|");
            printf("Child %d: ", i);
            print_tree_indented(tree, &tree->nodes[node->children[i]], tabs + 1,
                                node->children[i], body_count);
        }
    }
}

double get_distance(RVec3 *r1, RVec3 *r2) {
    return sqrt(pow(r1->x - r2->x, 2) + pow(r1->y - r2->y, 2) +
                pow(r1->z - r2->z, 2));
}

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

// now init_node manages firstfree variable and memory reallocation
int init_node(Octtree *tree) {
    // omp_set_lock(&tree->reallocnodes);
    int new_node;
// #pragma omp atomic capture
#pragma omp critical
    {
        new_node = tree->firstfree++;
        if (tree->sz <= tree->firstfree) {
            printf("No more space for new nodes! Exiting.\n");
            exit(1);
            // tree->sz *= 2;
            // tree->nodes = realloc(tree->nodes, tree->sz*sizeof(Octnode));
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
#pragma omp atomic update
                tree->nodes[node_indx].ents++;

                // When the leaves will be in different position exit the loop
                while (body_pos == other_indx) {

                    // take first free location and set the parent of the new
                    // branch
                    int free;

                    free = init_node(tree);
                    tree->nodes[free].parent = node_indx;

// set the new branch as child
#pragma omp atomic write
                    node->children[body_pos] = free;
#pragma omp atomic write // maybe not needed here
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
}

void add_ents(Octtree *tree, Entity *ents, uint ents_sz) {
#pragma omp for nowait
    for (int i = 0; i < ents_sz; i++) {
        add_ent(tree, &ents[i], i);
    }
#pragma omp barrier
}

void update_center_of_mass(Octnode *child, RVec3 *center, double *mass) {
    double new_mass = *mass + child->mass;

    center->x = (child->center.x * child->mass / new_mass) +
                (center->x * *mass / new_mass);
    center->y = (child->center.y * child->mass / new_mass) +
                (center->y * *mass / new_mass);
    center->z = (child->center.z * child->mass / new_mass) +
                (center->z * *mass / new_mass);

    *mass = new_mass;
}

void center_of_mass(Octtree *tree) {

#pragma omp for nowait
    for (int n = tree->firstfree - 1; n >= tree->root; n--) {
        Octnode *my_node = &tree->nodes[n];
        RVec3 l_center = my_node->center;
        double l_new_mass = my_node->mass;
        int j;

        int counter = 0;

        j = my_node->children[counter];
        while (counter < 8) {
            if (j == -1) {
                j = my_node->children[++counter];
            } else if (tree->nodes[j].mass != 0) {
                update_center_of_mass(&tree->nodes[j], &l_center, &l_new_mass);
                j = my_node->children[++counter];
            }
        }

        my_node->center = l_center;
        my_node->mass = l_new_mass;
    }
#pragma omp barrier
}

void get_bounding_box(Entity ents[], int ents_sz, double *max_val,
                      pad_double *loc_max) {
    int id = omp_get_thread_num();
    loc_max[id].val = 0.0;

#pragma omp for
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

#pragma omp single
    {
        int total_threads = omp_get_num_threads();
        for (int i = 0; i < total_threads; i++) {
            if (loc_max[i].val > *max_val)
                *max_val = loc_max[i].val;
        }
        *max_val *= 2.0;
    }
}

void create_tree(int ents_sz, Octtree *tree) {
    tree->firstfree = ents_sz;
    omp_init_lock(&tree->reallocnodes);
    tree->sz = ents_sz * 5;
    tree->nodes = malloc(sizeof(Octnode) * tree->sz);
    tree->max = 0;
    tree->root = ents_sz;
}

// calculate calculation caused by another body
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

// calculate body acceleration
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

// call recursive function get_acceleration_rec
void get_acceleration(Octtree *tree, RVec3 *acc, int ents_sz) {
#pragma omp for
    for (int i = 0; i < ents_sz; i++) {
        acc[i].x = 0;
        acc[i].y = 0;
        acc[i].z = 0;
        get_acceleration_rec(tree, tree->root, i, &acc[i], tree->max);
    }
}

void init_mutex(Octtree *tree) {
#pragma omp for
    for (size_t i = 0; i < tree->sz; i++) {
        for (int j = 0; j < 8; j++) {
            omp_init_lock(&tree->nodes[i].writelocks[j]);
        }
    }
}

void destroy_mutex(Octtree *tree) {
#pragma omp for
    for (size_t i = 0; i < tree->sz; i++) {
        for (int j = 0; j < 8; j++) {
            omp_destroy_lock(&tree->nodes[i].writelocks[j]);
        }
    }
}

void get_energy(Entity *ents, int sz, double *KE, double *PE,
                pad_double *local_KE, pad_double *local_PE) {

    // calculate KE
    int id = omp_get_thread_num();
    local_KE[id].val = 0.0;

#pragma omp for
    for (int i = 0; i < sz; i++) {
        RVec3 vel = ents[i].vel;
        double mass = ents[i].mass;

        local_KE[id].val += vel.x * vel.x * mass;
        local_KE[id].val += vel.x * vel.x * mass;
        local_KE[id].val += vel.x * vel.x * mass;
    }

#pragma omp single
    {
        *KE = 0.0;
        for (int i = 0; i < thread_count; i++)
            *KE += local_KE[i].val;
        *KE *= 0.5;
    }

    // calculate PE
    local_PE[id].val = 0.0;

#pragma omp for
    for (int i = 0; i < sz; i++) {
        RVec3 *i_pos = &ents[i].pos;
        for (int j = i; j < sz; j++) {
            double dx, dy, dz, D;
            RVec3 *j_pos = &ents[j].pos;

            dx = i_pos->x - j_pos->x;
            dy = i_pos->y - j_pos->y;
            dz = i_pos->z - j_pos->z;
            D = sqrt(dx * dx + dy * dy + dz * dz);
            D = D > 0 ? 1.0 / D : D;

            local_PE[id].val += -(ents[i].mass * ents[j].mass) * D;
        }
    }

#pragma omp single
    {
        *PE = 0.0;
        for (int i = 0; i < thread_count; i++)
            *PE += local_PE[i].val;
        *PE *= BIG_G;
    }
}

void save_energy(const char *output, int n_steps, double *KE, double *PE) {
    FILE *fpt;
    char energy_file[100];
    int l;

    *energy_file = '\0';
    l = strlen(output);
    strcat(energy_file, output);
    energy_file[l - 4] = '\0';
    strcat(energy_file, "_energy.csv");
    fpt = fopen(energy_file, "w");
    for (int i = 0; i < n_steps; i++)
        fprintf(fpt, "%lf,%lf,%lf\n", KE[i], PE[i], KE[i] + PE[i]);
    fclose(fpt);
}

// will calculate the bodies position over time
void propagation(Entity ents[], int ents_sz, int n_steps, float dt,
                 const char *output) {
    FILE *fpt;
    Octtree tree;
    RVec3 *acc;
    pad_double loc_max[thread_count];
    pad_double local_KE[thread_count];
    pad_double local_PE[thread_count];
    double *KE, *PE;

    create_tree(ents_sz, &tree);
    init_node(&tree);
    tree.max = 0.0;

    fpt = fopen(output, "w");

    for (int i = 0; i < thread_count; i++)
        loc_max[i].val = 0.0;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    KE = malloc(n_steps * sizeof(double));
    if (KE == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    PE = malloc(n_steps * sizeof(double));
    if (PE == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    // Initial positions
    for (size_t i = 0; i < ents_sz; i++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i, ents[i].pos.x, ents[i].pos.y,
                ents[i].pos.z, ents[i].mass);
    }

#pragma omp parallel num_threads(thread_count)
    {
        init_mutex(&tree);
        get_bounding_box(ents, ents_sz, &tree.max, loc_max);
        add_ents(&tree, ents, ents_sz);
        center_of_mass(&tree);
        get_acceleration(&tree, acc, ents_sz);

        for (int t = 0; t < n_steps; t++) {
            // 1/2 kick
#pragma omp for
            for (int i = 0; i < ents_sz; i++) {
                ents[i].vel.x += acc[i].x * dt / 2.0;
                ents[i].vel.y += acc[i].y * dt / 2.0;
                ents[i].vel.z += acc[i].z * dt / 2.0;
            }

            // Move bodies
#pragma omp for
            for (int i = 0; i < ents_sz; i++) {
                ents[i].pos.x += ents[i].vel.x * dt;
                ents[i].pos.y += ents[i].vel.y * dt;
                ents[i].pos.z += ents[i].vel.z * dt;
            }

#pragma omp single
            for (int i = 0; i < ents_sz; i++) {
                fprintf(fpt, "%d,%lf,%lf,%lf,%lf\n", i, ents[i].pos.x,
                        ents[i].pos.y, ents[i].pos.z, ents[i].mass);
            }

            // Build new tree
#pragma omp single
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
#pragma omp for
            for (int i = 0; i < ents_sz; i++) {
                ents[i].vel.x += acc[i].x * dt / 2.0;
                ents[i].vel.y += acc[i].y * dt / 2.0;
                ents[i].vel.z += acc[i].z * dt / 2.0;
            }
            get_energy(ents, ents_sz, KE + t, PE + t, local_KE, local_PE);
        }
        destroy_mutex(&tree);
    } // pragma

    fclose(fpt);
    save_energy(output, n_steps, KE, PE);
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
