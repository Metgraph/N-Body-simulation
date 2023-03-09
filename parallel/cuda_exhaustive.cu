// https://en.wikipedia.org/wiki/N-body_simulation
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

typedef unsigned int uint;

typedef struct {
    double x;
    double y;
    double z;
} RVec3;

const double BIG_G = 6.67e-11;

__host__ void get_entities(char filename[], uint *n_ents, double **masses, double **positions, double **velocities);
__host__ void count_entities_file(char *filename, uint *n);
__global__ void propagation(double *masses, double *positions, double *velocities, uint ents_sz, double dt);
__host__ void safe_malloc_double(double **pointer, int quantity);
__global__ void update_positions(double *positions, double *velocities, uint ents_sz, double dt);
__host__ void run(double *masses, double *positions, double *velocities, uint ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output);
__host__ void write_positions(double *positions, double *velocities, uint ents_sz, FILE *fpt);

int main(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
        exit(1);
    }

    uint n_ents;
    size_t start, end, dt;

    double *masses;
    double *positions;
    double *velocities;

    get_entities(argv[1], &n_ents, &masses, &positions, &velocities);
    start = strtoul(argv[2], NULL, 10);
    end = strtoul(argv[3], NULL, 10);
    dt = strtoul(argv[4], NULL, 10);
    run(masses, positions, velocities, n_ents, start, end, dt, argv[5]);

    free(masses);
    free(positions);
    free(velocities);
}

__host__
void cuda_check_error(cudaError_t error) {
    if (error != cudaSuccess) {
        fprintf(stderr, "Error: %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
}

__host__
void run(double *h_masses, double *h_positions, double *h_velocities, uint ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output) {
    FILE *fpt;
    cudaError_t error;
    double *d_masses;
    double *d_positions;
    double *d_velocities;

    error = cudaMalloc((void **)&d_masses, ents_sz * sizeof(double));
    cuda_check_error(error);
    error = cudaMalloc((void **)&d_positions, ents_sz * 3 * sizeof(double));
    cuda_check_error(error);
    error = cudaMalloc((void **)&d_velocities, ents_sz * 3 * sizeof(double));
    cuda_check_error(error);

    error = cudaMemcpy(d_masses, h_masses, ents_sz * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error);
    error = cudaMemcpy(d_positions, h_positions, ents_sz * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error);
    error = cudaMemcpy(d_velocities, h_velocities, ents_sz * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error);

//    printf("\nDEBUG: valori letti - parte c\n");
//    for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
//        fprintf(stderr,  "%lu,%lf,%lf,%lf,%lf,%lf,%lf\n",
//            entity_idx,
//            h_positions[entity_idx*3    ],
//            h_positions[entity_idx*3 + 1],
//            h_positions[entity_idx*3 + 2],
//            h_velocities[entity_idx*3    ],
//            h_velocities[entity_idx*3 + 1],
//            h_velocities[entity_idx*3 + 2]
//            );
//    }

    fpt = fopen(output, "w");
    for(size_t t = t_start; t < t_end; t += dt) { //Ciclo sul tempo
        propagation<<<1, 1>>>(d_masses, d_positions, d_velocities, ents_sz, (double) dt);
        cudaDeviceSynchronize();
        update_positions<<<1, 1>>>(d_positions, d_velocities, ents_sz, (double) dt);
        cudaDeviceSynchronize();

        // memcopy da device a host
        error = cudaMemcpy(h_positions, d_positions, ents_sz * 3 * sizeof(double), cudaMemcpyDeviceToHost);
        cuda_check_error(error);
        error = cudaMemcpy(h_velocities, d_velocities, ents_sz * 3 * sizeof(double), cudaMemcpyDeviceToHost);
        cuda_check_error(error);

        write_positions(h_positions, h_velocities, ents_sz, fpt);
    }

    cudaFree(d_masses);
    cudaFree(d_positions);
    cudaFree(d_velocities);

    fclose(fpt);
}

/**
 * Estimate the number of bodies by counting the lines of the file
 */
__host__
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
 * Allocate array for doubles and check for errors
 */
__host__
void safe_malloc_double(double **pointer, int quantity) {
    *pointer = (double *) malloc(quantity * sizeof(double));
    if (*pointer == NULL){
        fprintf(stderr, "Error during memory allocation\n");
        exit(2);
    }
}

/**
 * Read file and generate an arrays for masses, positions and velocities for the bodies
 */
__host__
void get_entities(char filename[], uint *n_ents, double **masses, double **positions, double **velocities) {
    FILE *file;
    int i;

    count_entities_file(filename, n_ents);

    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        exit(1);
    }

    safe_malloc_double(masses, *n_ents);
    safe_malloc_double(positions, *n_ents * 3);
    safe_malloc_double(velocities, *n_ents * 3);
    fprintf(stderr, "Aperto il file");
    i = 0;

    while ((fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
                        &((*positions)[i*3  ] ),
                        &((*positions)[i*3+1] ),
                        &((*positions)[i*3+2] ),
                        &((*velocities)[i*3  ]),
                        &((*velocities)[i*3+1]),
                        &((*velocities)[i*3+2]),
                        &((*masses)[i])       )) == 7) {
        i++;
    }
    // check if while ended because the end of the file has been reached
    if (fgetc(file) != EOF) {
        fprintf(stderr, "Error while reading file '%s': file is not well formed\n", filename);
        fclose(file);
        exit(1);
    }
    // Update n_ents with the correct number of scanned lines
    *n_ents = i;
    fclose(file);
    fprintf(stderr, "Scritti %d su memoria\n", *n_ents);
}

__global__
void propagation(double *masses, double *positions, double *velocities, uint ents_sz, double dt) {
    //const unsigned int tid = threadIdx.x;
    /* printf("Thread id: %u\n", tid); */
    /* printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf \n", masses[0], positions[0], positions[1], positions[2], velocities[0], velocities[1], velocities[3]); */
    printf("DEBUG: dati inizio kernel\n");
    for (int i=0; i < ents_sz; i++){
            printf("%u,%lf,%lf,%lf,%lf,%lf,%lf \n",
                i,
                positions[i*3],
                positions[i*3+1],
                positions[i*3+2],
                velocities[i*3],
                velocities[i*3+1],
                velocities[i*3+2]
                );
    }
    printf("DEBUG: delta: %f, #ents: %u\n", dt, ents_sz);

    for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++) { //Fisso un corpo
        RVec3 a_g = {0, 0, 0}; //vettore di appoggio
        RVec3 r_vector; // Spostato da sotto - lo alloco una sola volta
        for (size_t m2_idx = 0; m2_idx < ents_sz; m2_idx++) {
            if (m2_idx != m1_idx) {

                r_vector.x = positions[m1_idx*3  ] - positions[m2_idx*3  ];
                r_vector.y = positions[m1_idx*3+1] - positions[m2_idx*3+1];
                r_vector.z = positions[m1_idx*3+2] - positions[m2_idx*3+2];
                printf("DEBUG: Rvector: x: %f, y: %f, z: %f\n", r_vector.x, r_vector.y, r_vector.x);

                // distanza tra i due corpi
                double r_mag = sqrt(r_vector.x * r_vector.x + r_vector.y * r_vector.y + r_vector.z * r_vector.z);
                double pow_R_mag = pow(r_mag, 2.0);
                double acceleration = -1.0 * BIG_G * masses[m2_idx] / pow_R_mag;
                printf("DEBUG: BIG G: %f, masses m2: %f\n", BIG_G, masses[m2_idx]);
                printf("DEBUG: Pow_R_mag: %f, R_mag: %f, Acc: %f\n", pow_R_mag, r_mag, acceleration);

                r_vector.x = r_vector.x / r_mag; // Uso la stessa variabile invece di una nuova
                r_vector.y = r_vector.y / r_mag;
                r_vector.z = r_vector.z / r_mag;
                printf("DEBUG: updated vector Rvector: x: %f, y: %f, z: %f\n", r_vector.x, r_vector.y, r_vector.x);

                a_g.x += acceleration * r_vector.x;
                a_g.y += acceleration * r_vector.y;
                a_g.z += acceleration * r_vector.z;
                printf("DEBUG: a_g fine ciclo interno: x: %f, y: %f, z: %f\n", a_g.x, a_g.y, a_g.z);
            }
        }
        velocities[m1_idx*3    ] += a_g.x * dt;
        velocities[m1_idx*3 + 1] += a_g.y * dt;
        velocities[m1_idx*3 + 2] += a_g.z * dt;
        printf("DEBUG: Update velocities, end kernel: %f, %f, %f\n", velocities[m1_idx*3], velocities[m1_idx*3 + 1], velocities[m1_idx*3 + 2]);
        }
}

__global__
void update_positions(double *positions, double *velocities, uint ents_sz, double dt) {
    printf("DEBUG: Into Kernel for updating new position\n");
    for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
        positions[entity_idx*3    ] += velocities[entity_idx*3    ] * dt;
        positions[entity_idx*3 + 1] += velocities[entity_idx*3 + 1] * dt;
        positions[entity_idx*3 + 2] += velocities[entity_idx*3 + 2] * dt;
        printf("%lu,%lf,%lf,%lf,%lf,%lf,%lf\n",
            entity_idx,
            positions[entity_idx*3    ],
            positions[entity_idx*3 + 1],
            positions[entity_idx*3 + 2],
            velocities[entity_idx*3    ],
            velocities[entity_idx*3 + 1],
            velocities[entity_idx*3 + 2]
            );
    }
}

__host__
void write_positions(double *positions, double *velocities, uint ents_sz, FILE *fpt){
    for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
        fprintf(fpt,  "%lu,%lf,%lf,%lf,%lf,%lf,%lf\n",
            entity_idx,
            positions[entity_idx*3    ],
            positions[entity_idx*3 + 1],
            positions[entity_idx*3 + 2],
            velocities[entity_idx*3    ],
            velocities[entity_idx*3 + 1],
            velocities[entity_idx*3 + 2]
            );
    }
}
