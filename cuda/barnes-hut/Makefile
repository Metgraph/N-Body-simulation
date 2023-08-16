objects = barnes-hut.o bh_tree.o

all: $(objects)
	nvcc -g -G -Werror all-warnings $(objects) -o bh

%.o: %.cu
	nvcc -g -G -Werror all-warnings -I. -dc $< -o $@

clean:
	rm -f *.o app
