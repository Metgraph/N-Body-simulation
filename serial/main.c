#include "stdio.h"
#include "stdlib.h"
#include "math.h"

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
	RVec3 vec;
	double mass;
} Entity;

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
					   &e_buff.pos.y, &e_buff.pos.z, &e_buff.vec.x,
					   &e_buff.vec.y, &e_buff.vec.z, &e_buff.mass)) == 7)
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

void propagation(Entity ents[], uint ents_sz, size_t t_start, size_t t_end, size_t dt)
{
	const double BIG_G = 6.67e-11;
	size_t t;
	uint i, j;
	RVec3 vec3, unit_vec3, a_g;
	double r_mag, acceleration;
	for (t = t_start; t < t_end; t += dt)
	{
		for (i = 0; i < ents_sz; i++)
		{
			for (j = 0; j < ents_sz; j++)
			{
				if (i != j)
				{
					vec3.x = ents[i].pos.x - ents[j].pos.x;
					vec3.y = ents[i].pos.y - ents[j].pos.y;
					vec3.z = ents[i].pos.z - ents[j].pos.z;

					// distance
					r_mag = sqrt(vec3.x * vec3.x + vec3.y * vec3.y + vec3.z * vec3.z);

					acceleration = -1.0 * BIG_G * (ents[j].mass) / (r_mag * r_mag);

					unit_vec3.x = vec3.x / r_mag;
					unit_vec3.y = vec3.y / r_mag;
					unit_vec3.z = vec3.z / r_mag;

					a_g.x = unit_vec3.x * acceleration;
					a_g.y = unit_vec3.y * acceleration;
					a_g.z = unit_vec3.z * acceleration;
				}
			}
			ents[i].vec.x += a_g.x * dt;
			ents[i].vec.y += a_g.y * dt;
			ents[i].vec.z += a_g.z * dt;
		}

		for (uint i = 0; i < ents_sz; i++)
		{
			ents[i].pos.x += ents[i].vec.x * dt;
			ents[i].pos.y += ents[i].vec.y * dt;
			ents[i].pos.z += ents[i].vec.z * dt;
		}
	}
}

int main(int argc, char *argv[])
{
	uint n_ents;
	Entity *ents;
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s filename\n", argv[0]);
	}
	n_ents = get_entities(argv[1], &ents);
}
