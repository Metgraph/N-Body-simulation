// https://en.wikipedia.org/wiki/N-body_simulation
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

// const double BIG_G = 6.67e-11;
const double BIG_G = 1.0;

void get_entities(char filename[], Entity **ents, uint *n_ents);
void count_entities_file(char *filename, uint *n);
void propagation(Entity ents[], uint ents_sz, int n_steps, float dt,
                 const char *output);
void acceleration(uint ents_sz, Entity *ents, RVec3 *acc);

int main(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr,
                "Usage: %s input_filename start_time end_time delta_time "
                "output_filename\n",
                argv[0]);
        exit(1);
    }

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
}

/**
 * Estimate the number of bodies by counting the lines of the file
 *
 * @params *filename    Input filename
 * @params *n           Pointer on wich to save the count
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
 *
 * @param   filename    Input filename
 * @param   **ents      Array for bodies information storage
 * @param   *n_ents     Storage for bodies count
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

/*
 * Calculate bodies accelerations
 *
 * @param ents_sz   Number of bodies in the simulation
 * @param *ents     Array with bodies informations
 * @param *acc      Array for store new acceleration
 */
void acceleration(uint ents_sz, Entity *ents, RVec3 *acc) {

    for (size_t m1 = 0; m1 < ents_sz; m1++) {

        acc[m1].x = 0;
        acc[m1].y = 0;
        acc[m1].z = 0;

        for (size_t m2 = 0; m2 < ents_sz; m2++) {
            RVec3 r_vector;

            // Distance beetween body m1 and body m2
            r_vector.x = ents[m2].pos.x - ents[m1].pos.x;
            r_vector.y = ents[m2].pos.y - ents[m1].pos.y;
            r_vector.z = ents[m2].pos.z - ents[m1].pos.z;

            double inv_r3 = r_vector.x * r_vector.x + r_vector.y * r_vector.y +
                            r_vector.z * r_vector.z + 0.01;
            inv_r3 = pow(inv_r3, -1.5);

            // Update acceleration
            acc[m1].x += BIG_G * r_vector.x * inv_r3 * ents[m2].mass;
            acc[m1].y += BIG_G * r_vector.y * inv_r3 * ents[m2].mass;
            acc[m1].z += BIG_G * r_vector.z * inv_r3 * ents[m2].mass;
        }
    }
}

/*
 * Main function for execute the simulation
 *
 * @param *ents     Array with the bodies information
 * @param ents_sz   Total number of bodies
 * @param n_steps   Total of steps for simulation
 * @param dt        Time interval beetween one step and the next
 * @param *output    Output file name for simulation results (compile with -DRESULTS)
 */
void propagation(Entity *ents, uint ents_sz, int n_steps, float dt,
                 const char *output) {
#ifdef RESULTS
    FILE *fpt;
#endif
    RVec3 *acc;

    acc = malloc(ents_sz * sizeof(RVec3));
    if (acc == NULL) {
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }

    // Calculate initial acceleration
    acceleration(ents_sz, ents, acc);

#ifdef RESULTS
    fpt = fopen(output, "w");
    // Print to file initial state
    for (size_t i = 0; i < ents_sz; i++) {
        fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i,
                ents[i].pos.x, ents[i].pos.y,
                ents[i].pos.z, ents[i].mass);
    }
#endif

    for (int t = 0; t < n_steps; t++) {

        // First 1/2 kick
        for (size_t i = 0; i < ents_sz; i++) {
            ents[i].vel.x += acc[i].x * dt / 2.0;
            ents[i].vel.y += acc[i].y * dt / 2.0;
            ents[i].vel.z += acc[i].z * dt / 2.0;
        }

        // Move bodies
        for (size_t i = 0; i < ents_sz; i++) {
            ents[i].pos.x += ents[i].vel.x * dt;
            ents[i].pos.y += ents[i].vel.y * dt;
            ents[i].pos.z += ents[i].vel.z * dt;

#ifdef RESULTS
            // Save positions
            fprintf(fpt, "%lu,%lf,%lf,%lf,%lf\n", i,
                    ents[i].pos.x, ents[i].pos.y,
                    ents[i].pos.z, ents[i].mass);
#endif
        }

        // Update accelerations
        acceleration(ents_sz, ents, acc);

        // Second 1/2 kick
        for (size_t m1 = 0; m1 < ents_sz; m1++) {
            ents[m1].vel.x += acc[m1].x * dt / 2.0;
            ents[m1].vel.y += acc[m1].y * dt / 2.0;
            ents[m1].vel.z += acc[m1].z * dt / 2.0;
        }
    }
#ifdef RESULTS
    fclose(fpt);
#endif
    free(acc);
}

