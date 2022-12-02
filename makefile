CC=g++
MM=mpicxx
CFLAG= -Wall -I. -O3
TARGETS= fft_mpi fft_omp

all: $(TARGETS)

# Run all
runAll: all
	clear
	./main

# MPI
fft_mpi: fft_mpi.o microtime.o
	$(MM) -o $@ $^ -lm

fft_mpi.o:  main_mpi.c microtime.h
	$(MM) -o $@ $(CFLAG) -c $<

# OpenMP
fft_omp: main_omp.o microtime.o
	$(CC) -o $@ $^ -lm

fft_omp.o: main_omp.c microtime.h
	$(CC)  $(CFLAG) -c $<

# Timing
microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
