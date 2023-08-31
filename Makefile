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
	${COMPILER} -Wall -o ${BUILD_DIR}/exh_serial serial/exhaustive.c -lm -O2

exhaustive_debug:
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/exh_serial_debug serial/exhaustive.c -lm

exhaustive_openmp:
	${COMPILER} -Wall -o ${BUILD_DIR}/exh_mp openmp/exhaustive.c -lm -fopenmp -O2

exhaustive_openmp_debug:
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/exh_mp_debug openmp/exhaustive.c -lm -fopenmp

barneshut:
	${COMPILER} -Wall -o ${BUILD_DIR}/barnes-hut serial/barnes-hut.c -lm -O2

barneshut_debug:
	${COMPILER} -DRESULTS -Wall -o ${BUILD_DIR}/barnes-hut-debug serial/barnes-hut.c -lm

barneshut_openmp:
	${COMPILER} -Wall -o ${BUILD_DIR}/barnes-hut-omp openmp/barnes-hut.c -lm -fopenmp -O2

barneshut_openmp_debug:
	${COMPILER} -DRESULTS -Wall -o ${BUILD_DIR}/barnes-hut-omp-debug openmp/barnes-hut.c -lm -fopenmp

all:
	${COMPILER} -Wall -o ${BUILD_DIR}/exh_serial serial/exhaustive.c -lm
	${COMPILER} -Wall -o ${BUILD_DIR}/exh_mp openmp/exhaustive.c -lm -fopenmp
	${COMPILER} -Wall -o ${BUILD_DIR}/barnes-hut serial/barnes-hut.c -lm
	${COMPILER} -Wall -o ${BUILD_DIR}/barnes-hut-omp openmp/barnes-hut.c -lm -fopenmp

all_debug:
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/exh_serial_debug serial/exhaustive.c -lm
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/exh_mp_debug openmp/exhaustive.c -lm -fopenmp
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/barnes-hut-debug serial/barnes-hut.c -lm
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/barnes-hut-omp-debug openmp/barnes-hut.c -lm -fopenmp

barneshut_cuda:
	${CUDA_COMPILER} -o ${BUILD_DIR}/barnes-hut-cuda cuda/barnes-hut.cu -lm -O2

barneshutdebug_cuda:
	${CUDA_COMPILER} -g -G -o ${BUILD_DIR}/barnes-hut-cuda cuda/barnes-hut.cu -lm

exhaustive_cuda:
	${CUDA_COMPILER} -o ${BUILD_DIR}/exhaustive-cuda cuda/exhaustive.cu -lm -O2

play:
	python3 utils/pygame-show.py $(FILE)

clean:
	rm -rf build/*
