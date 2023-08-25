.PHONY : all serial
BUILD_DIR ?= build
COMPILER ?= gcc
CUDA_COMPILER ?= nvcc

.PHONY: all clean

ifeq ($(MAKECMDGOALS), play)
ifndef FILE
$(error Usage: 'make play FILE=<path to csv>')
endif
endif

exhaustive:
	${COMPILER} -Wall -o ${BUILD_DIR}/exhaustive_s serial/exhaustive.c -lm

exhaustivedebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/exhaustive_d serial/exhaustive.c -lm

manualcollapse:
	${COMPILER} -Wall  -o ${BUILD_DIR}/manual_collapse openmp/manual_collapse.c -lm -fopenmp

manualcollapsedebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/manual_collapse_d openmp/manual_collapse.c -lm -fopenmp

collapse:
	${COMPILER} -Wall -o ${BUILD_DIR}/collapse openmp/collapse.c -lm -fopenmp

collapsedebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/collapse_d openmp/collapse.c -lm -fopenmp

parallelfor:
	${COMPILER} -Wall -o ${BUILD_DIR}/parallel_for openmp/parallel_for.c -lm -fopenmp

parallelfordebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/parallel_for_d openmp/parallel_for.c -lm -fopenmp

barneshut:
	${COMPILER} -Wall -o ${BUILD_DIR}/barnes-hut serial/barnes-hut.c -lm
	${COMPILER} -Wall -o ${BUILD_DIR}/barnes-hut-omp openmp/barnes-hut.c -lm -fopenmp

barneshutdebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/barnes-hut serial/barnes-hut.c -lm
	${COMPILER} -Wall -g -o ${BUILD_DIR}/barnes-hut-omp openmp/barnes-hut.c -lm -fopenmp

cuda_barneshut:
	${CUDA_COMPILER} -o ${BUILD_DIR}/barnes-hut-cuda cuda/definitive_barnes-hut/barnes-hut.cu -lm -O2

cuda_barneshutdebug:
	${CUDA_COMPILER} -g -G -o ${BUILD_DIR}/barnes-hut-cuda cuda/definitive_barnes-hut/barnes-hut.cu -lm

play:
	python3 utils/pygame-show.py $(FILE)

clean:
	rm -rf build/*
