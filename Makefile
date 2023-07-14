.PHONY : all serial
BUILD_DIR ?= build
COMPILER ?= gcc

.PHONY: all clean

ifeq ($(MAKECMDGOALS), play)
ifndef FILE
$(error Usage: 'make play FILE=<path to csv>')
endif
endif

exhaustive:
	${COMPILER} -Wall -o ${BUILD_DIR}/exhaustive serial/exhaustive.c -lm

exhaustivedebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/exhaustive_d serial/exhaustive.c -lm

manualcollapse:
	${COMPILER} -Wall  -o ${BUILD_DIR}/manual_collapse openmp/manual_collapse.c -lm -fopenmp

manualcollapsedebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/manual_collapse_d openmp/manual_collapse.c -lm -fopenmp

collapse:
	${COMPILER} -Wall -o ${BUILD_DIR}/collapse openmp/mp_collapse.c -lm -fopenmp

collapsedebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/collapse_d openmp/mp_collapse.c -lm -fopenmp

parallelfor:
	${COMPILER} -Wall -o ${BUILD_DIR}/parallel_for openmp/mp_parallel_for.c -lm -fopenmp

parallelfordebug:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/parallel_for_d openmp/mp_parallel_for.c -lm -fopenmp

barneshut:
	${COMPILER} -Wall -g -o ${BUILD_DIR}/barnes-hut serial/barnes-hut.c -lm

play:
	python3 utils/pygame-show.py $(FILE)

clean:
	rm -rf build/*
