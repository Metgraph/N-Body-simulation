#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

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

// const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;
int thread_count;

void get_entities(char filename[], Entity **ents, uint *n_ents);
void count_entities_file(char *filename, uint *n);
void propagation(Entity ents[], uint ents_sz, int n_steps, float dt,
                 const char *output);
void acceleration(uint ents_sz, Entity *ents, RVec3 *acc);

int main(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename THREADS_NUM\n", argv[0]);
        exit(1);
    }

    thread_count = atoi(argv[6]);

    uint n_ents;
    Entity *ents;
    float start, end, dt;
    int n_steps;
    struct timespec s, e;

    get_entities(argv[1], &ents, &n_ents);

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
}

/**
 * Estimate the number of bodies by counting the lines of the file
 */
void count_entities_file(char *filename, uint *n) {
    FILE *file;
    char c;

    file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        exit(1);
    }

    *n = 0;
    for (c = getc(file); c != EOF; c = getc(file))
        if (c == '\n')
            (*n)++;
    fclose(file);
    if (n == 0) {
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
    if (*ents == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    ret_ptr = *ents;
    while ((fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &e_buff.pos.x,
                   &e_buff.pos.y, &e_buff.pos.z, &e_buff.vel.x, &e_buff.vel.y,
                   &e_buff.vel.z, &e_buff.mass)) == 7) {
        *ret_ptr++ = e_buff;
    }
    // check if while ended because the end of the file has been reached
    if (fgetc(file) != EOF) {
        fprintf(stderr,
                "Error while reading file '%s': file is not well formed\n",
                filename);
        fclose(file);
        exit(1);
    }
    // Update n_ents with the correct number of scanned lines
    *n_ents = ret_ptr - *ents;
    fclose(file);
}

void acceleration(uint ents_sz, Entity *ents, RVec3 *acc) {

    // Azzero prima per evitare race conditions con il collapse
#   pragma omp for
    for (size_t m1 = 0; m1 < ents_sz; m1++) {
        /* printf("Thread %d calculate m: %lu\n", omp_get_thread_num(), m1); */
        acc[m1].x = 0;
        acc[m1].y = 0;
        acc[m1].z = 0;
    }

#   pragma omp for collapse(2)
    for (size_t m1 = 0; m1 < ents_sz; m1++) {
        for (size_t m2 = 0; m2 < ents_sz; m2++) {
        /* printf("Thread %d calculate m1 m2: %lu, %lu\n", omp_get_thread_num(), m1, m2); */
            RVec3 r_vector;

            r_vector.x = ents[m2].pos.x - ents[m1].pos.x;
            r_vector.y = ents[m2].pos.y - ents[m1].pos.y;
            r_vector.z = ents[m2].pos.z - ents[m1].pos.z;

            double inv_r3 = r_vector.x * r_vector.x + r_vector.y * r_vector.y +
                            r_vector.z * r_vector.z + 0.01;
            inv_r3 = pow(inv_r3, -1.5);

            double ax = BIG_G * r_vector.x * inv_r3 * ents[m2].mass;
            double ay = BIG_G * r_vector.y * inv_r3 * ents[m2].mass;
            double az = BIG_G * r_vector.z * inv_r3 * ents[m2].mass;

            // Con atomic si hanno prestazioni simili al lock
#           pragma omp atomic
            acc[m1].x += ax;

#           pragma omp atomic
            acc[m1].y += ay;

#           pragma omp atomic
            acc[m1].z += az;

            /*
            // Con critical peggiora terribilmente
#           pragma omp critical
            acc[m1].x += ax;

#           pragma omp critical
            acc[m1].y += ay;

#           pragma omp critical
            acc[m1].z += az;
            */

        }
    }
}

void propagation(Entity *ents, uint ents_sz, int n_steps, float dt,
                 const char *output) {
    FILE *fpt;
    RVec3 *acc;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    fpt = fopen(output, "w");

    // Initial acceleration
#   pragma omp parallel num_threads(thread_count) // Fork unico
    {
        acceleration(ents_sz, ents, acc);
    }
    /* printf("\n"); */
    //printf("Inizializzazione terminata\n");

    // Print to file initial state
    for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", entity_idx,
                ents[entity_idx].pos.x, ents[entity_idx].pos.y,
                ents[entity_idx].pos.z, ents[entity_idx].mass);
    }

    // Con il for unico e utilizzando omp for invece di omp parallel for
    // Vengono forkati i thread solo all'inizio del ciclo e riutilizzati
    // invece di avere un fork ad ogni for
#   pragma omp parallel num_threads(thread_count)
    for (int t = 0; t < n_steps; t++) {

        // First 1/2 kick
#       pragma omp for
        for (size_t m1 = 0; m1 < ents_sz; m1++) {
            ents[m1].vel.x += acc[m1].x * dt / 2.0;
            ents[m1].vel.y += acc[m1].y * dt / 2.0;
            ents[m1].vel.z += acc[m1].z * dt / 2.0;
        }

        // Move bodies
#       pragma omp for
        for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
            ents[entity_idx].pos.x += ents[entity_idx].vel.x * dt;
            ents[entity_idx].pos.y += ents[entity_idx].vel.y * dt;
            ents[entity_idx].pos.z += ents[entity_idx].vel.z * dt;

        }

        // Save positions
        // Il primo thread che arriva lo esegue, nowait senza barriera
        // Non mi preoccupo di scrittura di dati svagliati perchè ogni blocco
        // ha una barriera implicita alla sua fine e la prossima istruzione
        // è accelerations che oltre a modificare le accelerazioni (che qui non servono)
        // non ha altri side effects
#       pragma omp single nowait
        for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
            //printf("Write file thread: %d\n", omp_get_thread_num());

            fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", entity_idx,
                ents[entity_idx].pos.x, ents[entity_idx].pos.y,
                ents[entity_idx].pos.z, ents[entity_idx].mass);
        }

        // Update accelerations
        acceleration(ents_sz, ents, acc);
        /* printf("\n"); */

        // Second 1/2 kick
#       pragma omp for
        for (size_t m1 = 0; m1 < ents_sz; m1++) {
            //printf("Thread %d\n", omp_get_thread_num());
            ents[m1].vel.x += acc[m1].x * dt / 2.0;
            ents[m1].vel.y += acc[m1].y * dt / 2.0;
            ents[m1].vel.z += acc[m1].z * dt / 2.0;
        }
    }
    fclose(fpt);
    free(acc);
}


