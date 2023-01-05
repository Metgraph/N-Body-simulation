// https://en.wikipedia.org/wiki/N-body_simulation
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

typedef unsigned int uint;

typedef struct
{
	double x;
	double y;
	double z;
} RVec3;

typedef struct
{
	RVec3 pos;
	RVec3 vel;
	double mass;
} Entity;

// TODO choose value for GravConstant
const double BIG_G = 6.67e-11;
// const double GravConstant = 0.01;

// Read file and generate an array of Entity
uint get_entities(char filename[], Entity **ents)
{
	Entity e_buff;
	int status;
	uint ret_size;
	uint size;
	Entity *ret;
	FILE *f = fopen(filename, "r");

	// Check if file has been open correctly, if not return NULL
	if (!f)
	{
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
					   &e_buff.vel.y, &e_buff.vel.z, &e_buff.mass)) == 7)
	{
		size++;
		if (ret_size < size)
		{
			// TODO Check for error in allocation
			ret_size *= 2;
			ret = (Entity *)realloc((void *)ret, ret_size * sizeof(Entity));
		}
		// Save value in first free location
		ret[size - 1] = e_buff;
	}

	// check if while ended because the end of the file has been reached
	if (fgetc(f) != EOF)
	{
		fprintf(stderr, "Error reading file '%s': file is not well formed\n",
				filename);
		fclose(f);
		return 0;
	}

	*ents = ret;
	fclose(f);
	return size;
}

void propagation(Entity ents[], uint ents_sz, size_t t_start, size_t t_end, size_t dt, const char *output)
{
	FILE *fpt;

	// alloc and initialize all value to 0 (0 in floating point is represented by all 0)
	RVec3 *a_g = (RVec3*)calloc(ents_sz, sizeof(RVec3));
	fpt = fopen(output, "w");
	for (size_t t = t_start; t < t_end; t += dt)
	{
#pragma omp parallel
		{
#pragma omp for collapse(2)
			for (size_t m1_idx = 0; m1_idx < ents_sz; m1_idx++)
			{

				for (size_t m2_idx = 0; m2_idx < ents_sz; m2_idx++)
				{
					if (m2_idx != m1_idx)
					{
						RVec3 r_vector;

						r_vector.x = ents[m1_idx].pos.x - ents[m2_idx].pos.x;
						r_vector.y = ents[m1_idx].pos.y - ents[m2_idx].pos.y;
						r_vector.z = ents[m1_idx].pos.z - ents[m2_idx].pos.z;

						double r_mag = sqrt(r_vector.x * r_vector.x + r_vector.y * r_vector.y + r_vector.z * r_vector.z);

						double acceleration = -1.0 * BIG_G * (ents[m2_idx].mass) / pow(r_mag, 2.0);

						RVec3 r_unit_vector = {r_vector.x / r_mag, r_vector.y / r_mag, r_vector.z / r_mag};

						a_g[m1_idx].x += acceleration * r_unit_vector.x;
						a_g[m1_idx].y += acceleration * r_unit_vector.y;
						a_g[m1_idx].z += acceleration * r_unit_vector.z;
					}
				}
			}

#pragma omp for
			for (size_t m_idx = 0; m_idx < ents_sz; m_idx++)
			{
				ents[m_idx].vel.x += a_g[m_idx].x * dt;
				ents[m_idx].vel.y += a_g[m_idx].y * dt;
				ents[m_idx].vel.z += a_g[m_idx].z * dt;
				ents[m_idx].pos.x += ents[m_idx].vel.x * dt;
				ents[m_idx].pos.y += ents[m_idx].vel.y * dt;
				ents[m_idx].pos.z += ents[m_idx].vel.z * dt;
			}
		}

		for (size_t entity_idx = 0; entity_idx < ents_sz; entity_idx++)
		{
			fprintf(fpt, "%u,%lf,%lf,%lf,%lf,%lf,%lf \n", entity_idx, ents[entity_idx].pos.x,
					ents[entity_idx].pos.y, ents[entity_idx].pos.z, ents[entity_idx].vel.x, ents[entity_idx].vel.y,
					ents[entity_idx].vel.z);
		}
	}
	free(a_g);
	fclose(fpt);
}

int main(int argc, char *argv[])
{

	uint n_ents;
	Entity *ents;
	size_t start, end, dt;
	if (argc != 6)
	{
		fprintf(stderr, "Usage: %s input_filename start_time end_time delta_time output_filename\n", argv[0]);
	}
	n_ents = get_entities(argv[1], &ents);
	start = strtoul(argv[2], NULL, 10);
	end = strtoul(argv[3], NULL, 10);
	dt = strtoul(argv[4], NULL, 10);
	propagation(ents, n_ents, start, end, dt, argv[5]);

	free(ents);
}
