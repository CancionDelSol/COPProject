#CC=g++
#CC=gcc
MM=mpicc

CFLAG= -Wall -I. -O3

TARGETS=fft-mpi fft-omp

all: $(TARGETS)

# Run all
runAll: all
	clear
	./main

# MPI
fft-mpi: main_mpi.o microtime.o
	$(MM) -o $@ $^ -lm

fft-mpi.o:  main_mpi.c microtime.h
	$(MM)  $(CFLAG) -c $<

# OpenMP
fft-omp: main-omp.o microtime.o
	$(CC) -o $@ $^ -lm

fft-omp.o: main-omp.c microtime.h
	$(CC)  $(CFLAG) -c $<


microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
