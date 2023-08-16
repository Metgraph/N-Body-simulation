objects = barnes-hut.o bh_tree.o

all: $(objects)
	nvcc -O0 -Xcicc -O0 -Xptxas -O0 --generate-line-info -Werror all-warnings $(objects) -o bh

%.o: %.cu
	nvcc -O0 -Xcicc -O0 -Xptxas -O0 --generate-line-info -x cu -Werror all-warnings -I. -dc $< -o $@

clean:
	rm -f *.o app
