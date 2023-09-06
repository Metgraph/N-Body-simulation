// https://en.wikipedia.org/wiki/N-body_simulation
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime_api.h>

typedef unsigned int uint;

#define BILLION 1000000000.0

//__constant__ double BIG_G = 6.67e-11; // Nella read-only cache e condivisa dai thread del blocco
__constant__ double BIG_G = 1.0; // Nella read-only cache e condivisa dai thread del blocco

__host__ void get_entities(char filename[], uint *n_ents, double **positions, double **velocities);
__host__ void count_entities_file(char *filename, uint *n);
__global__ void acceleration(double *positions, double *acc, uint ents_sz, int step);
__host__ void safe_malloc_double(double **pointer, int quantity);
__global__ void update_positions(double4 *positions, double3 *velocities, uint ents_sz, double dt, int step);
__global__ void update_velocities(double *acc, double *vel, uint ents_sz, double dt);
__host__ void propagation(double *positions, double *velocities, uint ents_sz, int n_steps, double dt, const char *output);

uint grid_s;
uint block_s;


int main(int argc, char *argv[]) {
    if (argc < 8) {
        fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename grid_size block_size\n", argv[0]);
        exit(1);
    }

    grid_s = atoi(argv[6]);
    block_s = atoi(argv[7]);

    uint n_ents;
    float start, end, dt;
    int n_steps;
    struct timespec s, e;

    double *positions;
    double *velocities;

    get_entities(argv[1], &n_ents, &positions, &velocities);

    start = strtof(argv[2], NULL);
    end = strtof(argv[3], NULL);
    dt = strtof(argv[4], NULL);

    n_steps = (end - start) / dt;

    printf("Start: %f, end: %f, delta time: %f, time steps: %d, ents: %d\n", start, end, dt, n_steps, n_ents);

    clock_gettime(CLOCK_REALTIME, &s);
    propagation(positions, velocities, n_ents, n_steps, dt, argv[5]);
    clock_gettime(CLOCK_REALTIME, &e);

    printf("Completed. Output file: %s\n", argv[5]);

    double time_spent =
        (e.tv_sec - s.tv_sec) + (e.tv_nsec - s.tv_nsec) / BILLION;

    printf("Elapsed wall time: %f s\n", time_spent);


    free(positions);
    free(velocities);
}


/*
 *  Check error after cuda functions.
 *
 *  @param error        Cuda error message
 *  @param *message     Message to print if there's an error
 */
__host__
void cuda_check_error(cudaError_t error, const char *message) {
    if (error != cudaSuccess) {
        fprintf(stderr, "Error: %s\nLine: %s", cudaGetErrorString(error), message);
        exit(EXIT_FAILURE);
    }
}

/*
 * Main function for execute the simulation
 *
 * @param *h_positions  Bodies positions e masses
 * @param *h_velocities Initial velocities
 * @param ents_sz       Total number of bodies
 * @param n_steps       Total of steps for simulation
 * @param dt            Time interval between one step and the next
 * @param *output       Output file name for simulation results (compile with -DRESULTS)
 */
__host__
void propagation(double *h_positions, double *h_velocities, uint ents_sz, int n_steps, double dt, const char *output) {
#ifdef RESULTS
    FILE *fpt;
#endif
    cudaError_t error;
    double *d_positions;
    double *d_velocities;
    double *d_accelerations;

    dim3 grid(grid_s, 1, 1);
    dim3 block(block_s, 1, 1);
    int shsz = block.x; // Shared memory size

    if (block_s * grid_s < ents_sz) {
        fprintf(stderr, "Insufficent threads!\nTotal cuda threads requested: %d - Total bodies: %d\n", block_s * grid_s, ents_sz);
        exit(1);
    }

    // Memory allocations and copies
    error = cudaMalloc((void **)&d_positions, ents_sz * 4 * sizeof(double) * (n_steps + 1));
    cuda_check_error(error, "Device positions malloc\n");

    error = cudaMalloc((void **)&d_velocities, ents_sz * 3 * sizeof(double));
    cuda_check_error(error, "Device velocities malloc\n");

    error = cudaMalloc((void **)&d_accelerations, ents_sz * 3 * sizeof(double));
    cuda_check_error(error, "Device accelerations malloc\n");

    error = cudaMemcpy(d_positions, h_positions, ents_sz * 4 * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error, "Host do Device copy positions\n");

    error = cudaMemcpy(d_velocities, h_velocities, ents_sz * 3 * sizeof(double), cudaMemcpyHostToDevice);
    cuda_check_error(error, "Host to Device copy velocities\n");


    acceleration<<<grid, block, shsz * sizeof(double) * 4>>>(d_positions, d_accelerations, ents_sz, 0);
    cudaDeviceSynchronize();

    for (int t = 1; t <= n_steps; t++) {

        // First 1/2 kick
        update_velocities<<<grid, block>>>(d_accelerations, d_velocities, ents_sz, dt);
        cudaDeviceSynchronize();

        // Drift
        update_positions<<<grid, block>>>((double4 *)d_positions, (double3 *)d_velocities, ents_sz, dt, t);
        cudaDeviceSynchronize();

        acceleration<<<grid, block, shsz * sizeof(double) * 4>>>(d_positions, d_accelerations, ents_sz, t);
        cudaDeviceSynchronize();

        // Second 1/2 kick
        update_velocities<<<grid, block>>>(d_accelerations, d_velocities, ents_sz, dt);
        cudaDeviceSynchronize();
    }

#ifdef RESULTS
    double *h_pos;
    safe_malloc_double(&h_pos, ents_sz * 4 * (n_steps+1));
    error = cudaMemcpy(h_pos, d_positions, (size_t)(ents_sz * 4 * sizeof(double) * (n_steps+1)), cudaMemcpyDeviceToHost);
    cuda_check_error(error, "Device to Host copy positions\n");
    double4 *h_poss = (double4 *) h_pos;

    fpt = fopen(output, "w");
    for (int i = 0; i < ents_sz * (n_steps+1); i++)
        fprintf(fpt, "%d,%lf,%lf,%lf,%lf\n", i%ents_sz, h_poss[i].x, h_poss[i].y, h_poss[i].z, h_poss[i].w);
    fclose(fpt);
    free(h_pos);
#endif

    cudaFree(d_positions);
    cudaFree(d_velocities);
    cudaFree(d_accelerations);

}

/**
 * Estimate the number of bodies by counting the lines of the file
 *
 * @params *filename    Input filename
 * @params *n           Pointer on wich to save the count
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
 *
 * @param **pointer     Destination pointer
 * @param quantity      Size of malloc
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
 * Read file and generate an array of Entity
 *
 * @param   filename        Input filename
 * @param   *n_ents         Storage for bodies count
 * @param   **positions     Storage for bodies positions
 * @param   **velocities    Storage for bodies velocities
 */
__host__
void get_entities(char filename[], uint *n_ents, double **positions, double **velocities) {
    FILE *file;
    int i;

    count_entities_file(filename, n_ents);

    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error while opening file '%s'\n", filename);
        exit(1);
    }

    safe_malloc_double(positions, *n_ents * 4);
    safe_malloc_double(velocities, *n_ents * 3);
    i = 0;

    while ((fscanf(file, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
                        &((*positions)[i*4  ] ),
                        &((*positions)[i*4+1] ),
                        &((*positions)[i*4+2] ),
                        &((*velocities)[i*3  ]),
                        &((*velocities)[i*3+1]),
                        &((*velocities)[i*3+2]),
                        &((*positions)[i*4+3]))) == 7) {
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
}

/*
 * Calculate bodies accelerations
 *
 * @param positions     Bodies positions
 * @param *acc          Array for store new acceleration
 * @param ents_sz       Number of bodies in the simulation
 * @param step          Number of current iteration into simulation
 */
__global__
void acceleration(double *positions, double *acc, uint ents_sz, int step) {
    extern __shared__ double4 shPositions[];
    double4 *gPositions;
    double3 *gAcc;
    double3 r_vector, l_acc;
    double inv_r3;
    int myId;
    int i;
    int iteration;
    int offset;
    double3 myBodyPos;

    // Move pointer to current positions of iteration
    positions = positions + (ents_sz * 4 * step);
    // Cast to double3-4 for more clear code later
    gPositions = (double4 *) positions;
    gAcc = (double3 *) acc;

    // Local acceleration
    l_acc = {0.0, 0.0, 0.0};
    myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        // Keep thread's body positions to avoid access to main memroy
        myBodyPos.x = gPositions[myId].x;
        myBodyPos.y = gPositions[myId].y;
        myBodyPos.z = gPositions[myId].z;
    }

    for (i = 0, iteration = 0; i < ents_sz; i += blockDim.x, iteration++) {
        // Offset for wich chunck of bodies retrieve in current iteration
        offset = blockDim.x * iteration;
        if (offset + threadIdx.x < ents_sz){
            // Each thread load a body into shared memory
            shPositions[threadIdx.x] = gPositions[offset + threadIdx.x];
        }
        __syncthreads(); // To ensure all threads load a body

        // Loop inside shared memory size && just on loaded bodies
        if (myId < ents_sz) {
            for (int j = 0; j < blockDim.x && offset+j < ents_sz; j++){
                // Calculate partial acceleration
                r_vector.x = shPositions[j].x - myBodyPos.x;
                r_vector.y = shPositions[j].y - myBodyPos.y;
                r_vector.z = shPositions[j].z - myBodyPos.z;

                inv_r3 = r_vector.x * r_vector.x + r_vector.y * r_vector.y +
                            r_vector.z * r_vector.z + 0.01;
                inv_r3 = pow(inv_r3, -1.5);

                // w attribute contain body's mass
                l_acc.x += BIG_G * r_vector.x * inv_r3 * shPositions[j].w;
                l_acc.y += BIG_G * r_vector.y * inv_r3 * shPositions[j].w;
                l_acc.z += BIG_G * r_vector.z * inv_r3 * shPositions[j].w;
            }
        }
        __syncthreads();
    }
    // Save new acceleration to global memory
    if (myId < ents_sz) {
        gAcc[myId].x = l_acc.x;
        gAcc[myId].y = l_acc.y;
        gAcc[myId].z = l_acc.z;
    }
}

/*
 * Update bodies velocities
 *
 * @param *acc      Accelerations array
 * @param *vel      Velocities array
 * @param ents_sz   Total number of bodies
 * @param dt        Time interval beetween one step and the next
 */
__global__
void update_velocities(double *acc, double *vel, uint ents_sz, double dt) {
    int myId;

    myId = blockIdx.x * blockDim.x + threadIdx.x;

    if (myId < ents_sz) {
        vel[myId*3    ] += acc[myId*3    ] * dt / 2.0;
        vel[myId*3 + 1] += acc[myId*3 + 1] * dt / 2.0;
        vel[myId*3 + 2] += acc[myId*3 + 2] * dt / 2.0;
    }
}

/*
 * Update bodies positions
 *
 * @param *gPos             Bodies positions buffer
 * @param *gVelocities      Velocities array
 * @param ents_sz           Total number of bodies
 * @param dt                Time interval beetween one step and the next
 * @param step              Number of current iteration into simulation
 */
__global__
void update_positions(double4 *gPos, double3 *gVelocities, uint ents_sz, double dt, int step) {
    int myId;
    int new_pos;
    int old_pos;

    myId = blockIdx.x * blockDim.x + threadIdx.x;
    // Calculate position into buffer to write new position
    new_pos = step * ents_sz + myId;
    // Previous position offset
    old_pos = (step - 1) * ents_sz + myId;

    if (myId < ents_sz) {
        gPos[new_pos].x = gPos[old_pos].x + gVelocities[myId].x * dt;
        gPos[new_pos].y = gPos[old_pos].y + gVelocities[myId].y * dt;
        gPos[new_pos].z = gPos[old_pos].z + gVelocities[myId].z * dt;
        gPos[new_pos].w = gPos[old_pos].w; // Mass of the body
    }
}

