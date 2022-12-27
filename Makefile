.PHONY : all serial
BUILD_DIR ?= build
COMPILER ?= gcc

.PHONY: all clean

serial: 
	${COMPILER} -o ${BUILD_DIR}/main serial/main.c -lm
	${COMPILER} -o ${BUILD_DIR}/openmp openmp/openmp.c -lm -fopenmp

debug:
	${COMPILER} -g -o ${BUILD_DIR}/main serial/main.c -lm
	${COMPILER} -g -o ${BUILD_DIR}/openmp openmp/openmp.c -lm -fopenmp
clean:
	rm -rf build/*
