CC=g++
MM=mpicxx
CFLAG= -Wall -I. -O3
TARGETS= fft_mpi fft_omp
PROCESS_COUNT= 5

all: $(TARGETS)

runMPI: fft_mpi
	mpirun -np $(PROCESS_COUNT) ./fft_mpi

runOMP: fft_omp
	./fft_omp

# MPI
fft_mpi: fft_mpi.o microtime.o
	$(MM) -o $@ $^ -lm

fft_mpi.o:  main_mpi.c microtime.h
	$(MM) -o $@ $(CFLAG) -c $<

# OpenMP
fft_omp: fft_omp.o microtime.o
	$(CC) -lm -fopenmp -o $@ $^ 

fft_omp.o: main_omp.c microtime.h
	$(CC) -fopenmp $(CFLAG) -o $@ -c $<

# Timing
microtime.o: microtime.c microtime.h
	$(CC) $(CFLAG) -c $<

clean:
	rm -f *.o *~ core $(TARGETS)
