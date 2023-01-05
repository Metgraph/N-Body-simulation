.PHONY : all serial
BUILD_DIR ?= build
COMPILER ?= gcc

.PHONY: all clean

serial: 
	${COMPILER} -o ${BUILD_DIR}/main serial/main.c -lm
	${COMPILER} -o ${BUILD_DIR}/openmp1loop openmp/openmp1loop.c -lm -fopenmp
	${COMPILER} -o ${BUILD_DIR}/openmp2loop openmp/openmp2loop.c -lm -fopenmp

debug:
	${COMPILER} -g -o ${BUILD_DIR}/main serial/main.c -lm
	${COMPILER} -o ${BUILD_DIR}/openmp1loop openmp/openmp1loop.c -lm -fopenmp
	${COMPILER} -g -o ${BUILD_DIR}/openmp2loop openmp/openmp2loop.c -lm -fopenmp
clean:
	rm -rf build/*
