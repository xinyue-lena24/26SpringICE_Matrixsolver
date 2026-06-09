CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -pedantic -O0 -g
CPPFLAGS = -Isrc/core -Isrc/util -Isrc/algorithm
LDFLAGS = -lm

CORE_SRCS = src/core/matrix_core.c src/core/matrix_ops.c
UTIL_SRCS = src/util/matrix_rand.c src/util/timer.c
ALGORITHM_SRCS = src/algorithm/matrix_lu.c src/algorithm/matrix_solve.c
LIB_SRCS = $(CORE_SRCS) $(UTIL_SRCS) $(ALGORITHM_SRCS)

TARGETS = main_basic main_timing_ops main_timing_mul

all: $(TARGETS)

main_basic: main_basic.c $(LIB_SRCS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ main_basic.c $(LIB_SRCS) $(LDFLAGS)

main_timing_ops: main_timing_ops.c $(LIB_SRCS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ main_timing_ops.c $(LIB_SRCS) $(LDFLAGS)

main_timing_mul: main_timing_mul.c $(LIB_SRCS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ main_timing_mul.c $(LIB_SRCS) $(LDFLAGS)

results:
	mkdir -p results

run: run-basic

run-basic: main_basic
	./main_basic

run-ops: main_timing_ops results
	./main_timing_ops

run-mul: main_timing_mul results
	./main_timing_mul

run-all: all results
	./main_basic
	./main_timing_ops
	./main_timing_mul

clean:
	rm -f $(TARGETS) *.o *.exe *.out
	rm -f results/*.csv

.PHONY: all results run run-basic run-ops run-mul run-all clean
