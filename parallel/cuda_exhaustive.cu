// https://en.wikipedia.org/wiki/N-body_simulation
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#define CUDA_CHECK_RETURN(value){ \
    cudaError_t_m_cudaStat = value; \
    if(_m_cudaStat!=cudaSuccess){ \
    fprintf(stderr,"Error %s at line %d infile %s\n",\
        cudaGetErrorString(_m_cudaStat), __LINE__,__FILE__); \
    exit(1); \
    }}

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
void print_error(cudaError_t error, int i) {
    if (error != cudaSuccess) {
        fprintf(stderr, "Error: %s, id: %d\n", cudaGetErrorString(error), i);
    }
}

__host__
void run(double *h_masses, double *h_positions, double *h_velocities, uint ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output) {
    FILE *fpt;
    cudaError_t error;
    double *d_masses;
    double *d_positions;
    double *d_velocities;

    // TODO: check errori
    error = cudaMalloc((void **)&d_masses, ents_sz * sizeof(double));
    print_error(error, 1);
    error = cudaMalloc((void **)&d_positions, ents_sz * 3 * sizeof(double));
    print_error(error, 2);
    error = cudaMalloc((void **)&d_velocities, ents_sz * 3 * sizeof(double));
    print_error(error, 3);

    error = cudaMemcpy(d_masses, h_masses, ents_sz, cudaMemcpyHostToDevice);
    print_error(error, 4);
    error = cudaMemcpy(d_positions, h_positions, ents_sz * 3, cudaMemcpyHostToDevice);
    print_error(error, 5);
    error = cudaMemcpy(d_velocities, h_velocities, ents_sz * 3, cudaMemcpyHostToDevice);
    print_error(error, 6);

    fpt = fopen(output, "w");
    for(size_t t = t_start; t<t_end; t+=dt) { //Ciclo sul tempo
        propagation<<<1, 1>>>(d_masses, d_positions, d_velocities, ents_sz, (double) dt);
        cudaDeviceSynchronize();
        update_positions<<<1, 1>>>(d_positions, d_velocities, ents_sz, (double) dt);
        cudaDeviceSynchronize();

        // memcopy da device a host
        error = cudaMemcpy(h_masses, d_masses, ents_sz, cudaMemcpyDeviceToHost);
        print_error(error, 7);
        error = cudaMemcpy(h_positions, d_positions, ents_sz * 3, cudaMemcpyDeviceToHost);
        print_error(error, 8);
        error = cudaMemcpy(h_velocities, d_velocities, ents_sz * 3, cudaMemcpyDeviceToHost);
        print_error(error, 9);

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
    const unsigned int tid = threadIdx.x;
    printf("Thread id: %u\n", tid);
    printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf \n", masses[0], positions[0], positions[1], positions[2], velocities[0], velocities[1], velocities[3]);
//    for (int i=0; i < ents_sz; i++){
//            printf("%u,%lf,%lf,%lf,%lf,%lf,%lf \n", i,
//                positions[i*3],
//                positions[i*3+1],
//                positions[i*3+2],
//                velocities[i*3],
//                velocities[i*3+1],
//                velocities[i*3+2]
//                );
//    }

    for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++) { //Fisso un corpo
        RVec3 a_g = {0, 0, 0}; //vettore di appoggio
        RVec3 r_vector; // Spostato da sotto - lo alloco una sola volta
        for (size_t m2_idx = 0; m2_idx < ents_sz; m2_idx++) {
            if (m2_idx != m1_idx) {

                r_vector.x = positions[m1_idx*3  ] - positions[m2_idx*3  ];
                r_vector.y = positions[m1_idx*3+1] - positions[m2_idx*3+1];
                r_vector.z = positions[m1_idx*3+2] - positions[m2_idx*3+2];

                // distanza tra i due corpi
                double r_mag = sqrt(r_vector.x * r_vector.x + r_vector.y * r_vector.y + r_vector.z * r_vector.z);

                double acceleration = -1.0 * BIG_G * masses[m2_idx] / pow(r_mag, 2.0);

                r_vector.x = r_vector.x / r_mag; // Uso la stessa variabile invece di una nuova
                r_vector.y = r_vector.y / r_mag;
                r_vector.z = r_vector.z / r_mag;

                a_g.x += acceleration * r_vector.x;
                a_g.y += acceleration * r_vector.y;
                a_g.z += acceleration * r_vector.z;
            }
        }
        velocities[m1_idx*3    ] += a_g.x * dt;
        velocities[m1_idx*3 + 1] += a_g.y * dt;
        velocities[m1_idx*3 + 2] += a_g.z * dt;
        }
}

__global__
void update_positions(double *positions, double *velocities, uint ents_sz, double dt) {
    for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
        positions[entity_idx*3    ] += velocities[entity_idx*3    ] * dt;
        positions[entity_idx*3 + 1] += velocities[entity_idx*3 + 1] * dt;
        positions[entity_idx*3 + 2] += velocities[entity_idx*3 + 2] * dt;
    }
}

__host__
void write_positions(double *positions, double *velocities, uint ents_sz, FILE *fpt){
    for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++) {
        fprintf(fpt,  "%u,%lf,%lf,%lf,%lf,%lf,%lf \n",
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
