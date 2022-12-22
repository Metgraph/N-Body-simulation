// https://rosettacode.org/wiki/N-body_problem#C
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

//TODO choose value for GravConstant
const double GravConstant = 6.67e-11;
// const double GravConstant = 0.01;

// Sum vectors
RVec3 vec3_sum(RVec3 *v1, RVec3 *v2) {
	RVec3 sum = {v1->x + v2->x, v1->y + v2->y, v1->z + v2->z};
	return sum;
}

RVec3 vec3_sub(RVec3 *v1, RVec3 *v2) {
	RVec3 sub = {v1->x - v2->x, v1->y - v2->y, v1->z - v2->z};
	return sub;
}

RVec3 vec3_scale(double scalar, RVec3 *v) {
	RVec3 prod = {v->x * scalar, v->y * scalar, v->z * scalar};
	return prod;
}

double vec3_mod(RVec3 *v) {
	return sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
}

// Read file and generate an array of Entity
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

void resolveCollisions(Entity *ents, uint ents_sz) {
	uint i, j;

	for (i = 0; i < ents_sz - 1; i++)
		for (j = i + 1; j < ents_sz; j++) {
			if (ents[i].pos.x == ents[j].pos.x &&
				ents[i].pos.y == ents[j].pos.y &&
				ents[i].pos.z == ents[j].pos.z) {
				RVec3 temp = ents[i].vel;
				ents[i].vel = ents[j].vel;
				ents[j].vel = temp;
			}
		}
}

void computeAccelerations(Entity *ents, RVec3 *accels, uint ents_sz) {
	uint i, j;
	RVec3 accel, sub_ji, accel_from_j;
	for (i = 0; i < ents_sz; i++) {
		accels[i].x = 0;
		accels[i].y = 0;
		accels[i].z = 0;
		for (j = 0; j < ents_sz; j++) {
			if (i != j) {
				sub_ji = vec3_sub(&ents[j].pos, &ents[i].pos);
				accel_from_j = vec3_scale(
					GravConstant * ents[j].mass / pow(vec3_mod(&sub_ji), 3),
					&sub_ji);
				accels[i] = vec3_sum(&accels[i], &accel_from_j);
			}
		}
	}
}

void computeVelocities(Entity *ents, RVec3 *accels, uint ents_sz, size_t dt) {
	uint i;
	RVec3 diff_vel;
	for (i = 0; i < ents_sz; i++) {
		diff_vel = vec3_scale((double)dt, &accels[i]);
		ents[i].vel = vec3_sum(&ents[i].vel, &diff_vel);
	}
}

void computePositions(Entity *ents, RVec3 *accels, uint ents_sz, size_t dt) {
	uint i;
	RVec3 space_from_acc, space_from_vel, movement;
	for (i = 0; i < ents_sz; i++) {
		space_from_acc = vec3_scale(0.5 * (double)(dt * dt), &accels[i]);
		space_from_vel = vec3_scale((double)dt, &ents[i].vel);
		movement = vec3_sum(&space_from_acc, &space_from_vel);
		ents[i].pos = vec3_sum(&ents[i].pos, &movement);
	}
}

void propagation(Entity ents[], uint ents_sz, size_t t_start, size_t t_end,
				 size_t dt, const char *output) {
	const double BIG_G = 6.67e-11;
	RVec3 *accels = (RVec3 *)malloc(ents_sz * sizeof(RVec3));
	FILE *fpt;
	fpt = fopen(output, "w");
	for (size_t t = t_start; t < t_end; t += dt) {
		computeAccelerations(ents, accels, ents_sz);
		computePositions(ents, accels, ents_sz, dt);
		computeVelocities(ents, accels, ents_sz, dt);
		resolveCollisions(ents, ents_sz);
        for(uint i=0; i<ents_sz; i++)
            fprintf(fpt, "%u,%lf,%lf,%lf,%lf,%lf,%lf \n", i, ents[i].pos.x,
            ents[i].pos.y, ents[i].pos.z, ents[i].vel.x,ents[i].vel.y,
            ents[i].vel.z);
	}

	fclose(fpt);
    free(accels);
}

int main(int argc, char *argv[]) {
	uint n_ents;
	Entity *ents;
    size_t start, end, dt;
	if (argc != 6) {
		fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
	}
	n_ents = get_entities(argv[1], &ents);
    start=strtoul(argv[2], NULL, 10);
    end=strtoul(argv[3], NULL, 10);
    dt=strtoul(argv[4], NULL, 10);
	propagation(ents, n_ents, start, end, dt, argv[5]);
    free(ents);
}
