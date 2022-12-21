.PHONY : all serial
BUILD_DIR ?= build
COMPILER ?= gcc

.PHONY: all clean

serial: 
	${COMPILER} -o ${BUILD_DIR}/main serial/main.c -lm

clean:
	rm -rf build/*