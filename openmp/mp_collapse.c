// https://en.wikipedia.org/wiki/N-body_simulation
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

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

const double BIG_G = 6.67e-11;
int thread_count;

void get_entities(char filename[], Entity **ents, uint *n_ents);
void count_entities_file(char *filename, uint *n);
void propagation(Entity ents[], uint ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output);

int main(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename THREADS_NUM\n", argv[0]);
        exit(1);
    }

    thread_count = atoi(argv[6]);

    uint n_ents;
    Entity *ents;
    size_t start, end, dt;

    get_entities(argv[1], &ents, &n_ents);
    start = strtoul(argv[2], NULL, 10);
    end = strtoul(argv[3], NULL, 10);
    dt = strtoul(argv[4], NULL, 10);
    propagation(ents, n_ents, start, end, dt, argv[5]);

    free(ents);
}

/**
 * Estimate the number of bodies by counting the lines of the file
 */
void count_entities_file(char *filename, uint *n){
    FILE *file;
    char c;

    file = fopen(filename, "r");
    if (file == NULL){
        fprintf(stderr, "Error opening file '%s'\n", filename);
        exit(1);
    }

    *n = 0;
    for (c = getc(file); c != EOF; c = getc(file))
        if (c == '\n')
            (*n)++;
    fclose(file);
    if (n == 0){
        fprintf(stderr, "No bodies found into file. Closing\n");
        exit(1);
    } else {
        // For prevent files that do not have a newline character at the end
        (*n)++;
    }
}

/**
 * Read file and generate an array of Entity
 */
void get_entities(char filename[], Entity **ents, uint *n_ents) {
    Entity e_buff;
    Entity *ret_ptr;
    FILE *file;

    count_entities_file(filename, n_ents);

    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        exit(1);
    }

    *ents = malloc(*n_ents * sizeof(Entity));
    if (*ents == NULL){
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    ret_ptr = *ents;
    while ((fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &e_buff.pos.x,
                    &e_buff.pos.y, &e_buff.pos.z, &e_buff.vel.x,
                    &e_buff.vel.y, &e_buff.vel.z, &e_buff.mass)) == 7) {
        *ret_ptr++ = e_buff;
    }
    // check if while ended because the end of the file has been reached
    if (fgetc(file) != EOF) {
        fprintf(stderr, "Error while reading file '%s': file is not well formed\n", filename);
        fclose(file);
        exit(1);
    }
    // Update n_ents with the correct number of scanned lines
    *n_ents = ret_ptr - *ents;
    fclose(file);
}

void propagation(Entity ents[], uint ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output) {
    FILE *fpt;
    RVec3 *a_g;
    omp_lock_t *locks;
//    for (int i = 0; i < ents_sz; i++){
//        printf("%d,%lf,%lf,%lf,%lf,%lf,%lf \n", i,
//                    ents[i].pos.x, ents[i].pos.y, ents[i].pos.z,
//                    ents[i].vel.x,ents[i].vel.y, ents[i].vel.z);
//    }

	a_g = calloc(ents_sz, sizeof(RVec3)); // Fare malloc e poi far azzerare ai thread?
    if (a_g == NULL){
        fprintf(stderr, "Error during memory allocation into propagation function\n");
        exit(2);
    }

	locks = calloc(ents_sz, sizeof(omp_lock_t)); // Fare malloc e poi far azzerare ai thread?
    if (a_g == NULL){
        fprintf(stderr, "Error during memory allocation into propagation function\n");
        exit(2);
    }


#   pragma omp parallel num_threads(thread_count) // Initialize locks
    for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++) {
        omp_init_lock(&locks[m1_idx]);
    }

    fpt = fopen(output, "w");
    /* printf("Thread richiesti: %d\n", thread_count); */
    /* printf("Thread totali: %d\n", omp_get_num_threads()); */
    for(size_t t = t_start; t<t_end; t+=dt) {
#       pragma omp parallel num_threads(thread_count) // Fork unico per tutti i cicli

#       pragma omp parallel for // Azzeramento prima di ogni ciclo
        for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++) {
            a_g[m1_idx].x = 0;
            a_g[m1_idx].y = 0;
            a_g[m1_idx].z = 0;
        }
#       pragma omp parallel for collapse(2)
        for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++) {
            for (size_t m2_idx = 0; m2_idx < ents_sz; m2_idx++) {
                if (m2_idx != m1_idx) {
                    RVec3 r_vector;

                    r_vector.x = ents[m1_idx].pos.x - ents[m2_idx].pos.x;
                    r_vector.y = ents[m1_idx].pos.y - ents[m2_idx].pos.y;
                    r_vector.z = ents[m1_idx].pos.z - ents[m2_idx].pos.z;
                    /* printf("m2: %f, %f, %f\n", ents[m2_idx].pos.x, ents[m2_idx].pos.y, ents[m2_idx].pos.z); */
                    /* printf("%lu vs %lu - %f, %f, %f\n", m1_idx, m2_idx, r_vector.x, r_vector.y, r_vector.z); */
                    //distanza tra i due corpi
                    double r_mag = sqrt(r_vector.x * r_vector.x + r_vector.y * r_vector.y + r_vector.z * r_vector.z);

                    double acceleration = -1.0 * BIG_G * (ents[m2_idx].mass) / pow(r_mag, 2.0);

                    RVec3 r_unit_vector = {r_vector.x / r_mag, r_vector.y / r_mag, r_vector.z / r_mag};

                    omp_set_lock(&locks[m1_idx]);
                    a_g[m1_idx].x += acceleration * r_unit_vector.x;
                    a_g[m1_idx].y += acceleration * r_unit_vector.y;
                    a_g[m1_idx].z += acceleration * r_unit_vector.z;
                    omp_unset_lock(&locks[m1_idx]);
                }
            }
        }
        /* printf("Sono il thread: %d\n", omp_get_thread_num()); */

#       pragma omp parallel for
        for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
            ents[entity_idx].vel.x += a_g[entity_idx].x * dt;
            ents[entity_idx].vel.y += a_g[entity_idx].y * dt;
            ents[entity_idx].vel.z += a_g[entity_idx].z * dt;
            ents[entity_idx].pos.x += ents[entity_idx].vel.x * dt;
            ents[entity_idx].pos.y += ents[entity_idx].vel.y * dt;
            ents[entity_idx].pos.z += ents[entity_idx].vel.z * dt;
        }
#       pragma omp single nowait // Il primo thread che arriva lo esegue, nowait senza barriera
        for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
            fprintf(fpt, "%lu,%lf,%lf,%lf,%lf,%lf,%lf\n", entity_idx,
                    ents[entity_idx].pos.x, ents[entity_idx].pos.y, ents[entity_idx].pos.z,
                    ents[entity_idx].vel.x,ents[entity_idx].vel.y, ents[entity_idx].vel.z);
        }
    }

#   pragma omp parallel num_threads(thread_count)
    for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++) {
        omp_destroy_lock(&locks[m1_idx]);
    }

	free(a_g);
    free(locks);
    fclose(fpt);
}

