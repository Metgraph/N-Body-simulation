.PHONY : all serial
BUILD_DIR ?= build
COMPILER ?= gcc

.PHONY: all clean

exhaustive: 
	${COMPILER} -o ${BUILD_DIR}/exhaustive serial/exhaustive.c -lm
	${COMPILER} -o ${BUILD_DIR}/openmp1loop openmp/openmp1loop.c -lm -fopenmp
	${COMPILER} -o ${BUILD_DIR}/openmp2loop openmp/openmp2loop.c -lm -fopenmp

debug:
	${COMPILER} -g -o ${BUILD_DIR}/exhaustive serial/exhaustive.c -lm
	${COMPILER} -g -o ${BUILD_DIR}/openmp1loop openmp/openmp1loop.c -lm -fopenmp
	${COMPILER} -g -o ${BUILD_DIR}/openmp2loop openmp/openmp2loop.c -lm -fopenmp

barneshut:
	${COMPILER} -g -o ${BUILD_DIR}/barnes-hut serial/barnes-hut.c -lm
	${COMPILER} -g -o ${BUILD_DIR}/barnes-hut-omp openmp/barnes-hut.c -lm -fopenmp


clean:
	rm -rf build/*
