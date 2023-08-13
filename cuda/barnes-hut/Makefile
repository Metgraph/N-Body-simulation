objects = barnes-hut.o bh_tree.o

all: $(objects)
	nvcc -Werror all-warnings $(objects) -o bh

%.o: %.cu
	nvcc -x cu -Werror all-warnings -I. -dc $< -o $@

clean:
	rm -f *.o app
