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

exhaustive_serial:
	${COMPILER} -Wall -o ${BUILD_DIR}/exh_serial serial/exhaustive.c -lm

exhaustive_serial_debug:
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/exh_serial_debug serial/exhaustive.c -lm

exhaustive_openmp:
	${COMPILER} -Wall -o ${BUILD_DIR}/exh_mp openmp/exhaustive.c -lm -fopenmp

exhaustive_openmp_debug:
	${COMPILER} -DRESULTS -Wall -g -o ${BUILD_DIR}/exh_mp_debug openmp/exhaustive.c -lm -fopenmp

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
