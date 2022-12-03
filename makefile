CC=g++
MM=mpicxx
CFLAG= -Wall -I. -O3
TARGETS= fft_mpi fft_omp
PROCESS_COUNT= 5

all: $(TARGETS)

runMPI: fft_mpi
	mpirun -np $(PROCESS_COUNT) ./fft_mpi

# MPI
fft_mpi: fft_mpi.o microtime.o
	$(MM) -o $@ $^ -lm

fft_mpi.o:  main_mpi.c microtime.h
	$(MM) -o $@ $(CFLAG) -c $<

# OpenMP
fft_omp: main_omp.o microtime.o
	$(CC) -o $@ $^ -lm -fopenmp

fft_omp.o: main_omp.c microtime.h
	$(CC) $(CFLAG) -fopenmp -c $<

# Timing
microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
