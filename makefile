CC=g++
MM=mpicxx
CFLAG= -Wall -I. -O3
TARGETS= fft_mpi fft_omp fft_mpi_col

all: $(TARGETS)

runMPI: fft_mpi
	#mpirun -np 2 ./fft_mpi
	#mpirun -np 4 ./fft_mpi
	#mpirun -np 5 ./fft_mpi

runMPICol: fft_mpi_col
	mpirun -np 2 ./fft_mpi_col
	mpirun -np 4 ./fft_mpi_col
	mpirun -np 5 ./fft_mpi_col

runOMP: fft_omp
	./fft_omp 2
	./fft_omp 4
	./fft_omp 5

# MPI-PEER
fft_mpi: fft_mpi.o microtime.o
	$(MM) -o $@ $^ -lm

fft_mpi.o:  main_mpi.c microtime.h
	$(MM) -o $@ $(CFLAG) -c $<

# MPI-COLLECTIVE
fft_mpi_col: fft_mpi_col.o microtime.o
	$(MM) -o $@ $^ -lm

fft_mpi_col.o:  main_mpi_col.c microtime.h
	$(MM) $(CFLAG) -o $@ -c $<

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
