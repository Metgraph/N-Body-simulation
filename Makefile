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

clean:
	rm -rf build/*

cleanclust:
	rm -rf exhaustive.*
	rm -rf results/*.csv

buildclust:
	nvcc -o exhaustive -Xcompiler -Wall cuda_exhaustive.cu

clusterall:
	rm -rf exhaustive.*
	rm -rf results/*.csv
	nvcc -o exhaustive -Xcompiler -Wall cuda_exhaustive.cu
