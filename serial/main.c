#include "stdio.h"
#include "stdlib.h"

typedef unsigned int uint;

typedef struct {
	double x;
	double y;
	double z;
} RVec3;

typedef struct {
	RVec3 pos;
	RVec3 vec;
	double mass;
} Entity;


//Read file and generate an array of Entity
uint get_entities(char filename[], Entity** ents) {
	FILE* f = fopen(filename, 'r');

	// Check if file has been open correctly, if not return NULL
	if (!f) {
		fprintf(stderr, "Error opening file '%s'\n", filename);
		return NULL;
	}

	// TODO Check for error in allocation
	Entity* ret = (Entity*)malloc(1 * sizeof(Entity));

	Entity e_buff;
	int status;
	uint size = 0;
	uint ret_size = 1;
	// fscanf return the number of input items successfully matched and assigned
	while ((status =
				fscanf(f, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", &e_buff.pos.x,
					   &e_buff.pos.y, &e_buff.pos.z, &e_buff.vec.x,
					   &e_buff.vec.y, &e_buff.vec.z, &e_buff.mass)) == 7) {
		size++;
		if (ret_size <= size) {
			// TODO Check for error in allocation
			ret_size *= 2;
			ret = (Entity*)realloc((void*)ret, ret_size);
		}
		// Save value in first free location
		ret[size - 1] = e_buff;
	}

	// check if while ended because the end of the file has been reached
	if (fgetc(f) != EOF) {
		fprintf(stderr, "Error reading file '%s': file is not well formed\n",
				filename);
		fclose(f);
		return NULL;
	}
    
    *ents=ret;
	fclose(f);
    return size;
}

int main(int argc, char* argv[]){
    if(argc!=2){
        fprintf(stderr, "Usage: %s filename\n", argv[0]);
    }

}
