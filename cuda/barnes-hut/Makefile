objects = barnes-hut.o bh_tree.o

all: $(objects)
	nvcc --generate-line-info -Werror all-warnings $(objects) -o bh

%.o: %.cu
	nvcc --generate-line-info -x cu -Werror all-warnings -I. -dc $< -o $@

clean:
	rm -f *.o app
