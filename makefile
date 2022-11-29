#CC=g++
CC=gcc
CFLAG= -Wall -I. -O3

TARGETS=fft fft-mpi fft-omp

all: $(TARGETS)

# Run all
runAll: all
	clear
	./main

# Original unoptimized version
fft: main.o microtime.o
	$(CC) -o $@ $^ -lm

# fft.o: main.cpp microtime.h
fft.o: main.c microtime.h
	$(CC)  $(CFLAG) -c $<

# MPI
fft-mpi: main-mpi.o microtime.o
	$(CC) -o $@ $^ -lm

fft-mpi.o: main-mpi.c microtime.h
	$(CC)  $(CFLAG) -c $<

# OpenMP
fft-omp: main-omp.o microtime.o
	$(CC) -o $@ $^ -lm

fft-omp.o: main-omp.c microtime.h
	$(CC)  $(CFLAG) -c $<
	

microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
