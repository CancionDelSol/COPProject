clear
export TMPDIR=\tmp
mpicc -g -Wall -o3 -o main-mpi main-mpi.c microtime.c -lm
mpiexec -n 2 ./main-mpi

# Project Team Members:

- Roger Johnson
- Doug Woodall
- Shaib Alagily
- Peyman Samimi
- Sadaf Charkhabi

# Project Description

We aim to perform a Fourier Transform on simple harmonic functions, knowing what the output should be. The FT is used to transform a composite function into its constituent frequencies.

In the continuous domain, this is completed through the use of an integral over infinite range. In computational science, this must be discretized. The algorithm involves a nested loop, with the same element count for each loop.

The inner loop performs a multiplication between an element inside the source array (base function output) and an exponential (Euler representation of the harmonic function). Each element of the output array is independent of the element prior or subsequent, thus rendering it parallelizable.

# Serial Implementation

Run from the console:

```
cd Serial
make fft
./fft
```

# OpenMP

Run from the console:

```
make runOMP
```

# MPI Gather

Run from the console:

```
make runMPICol
```

# MPI Send/Receive

```
make runMPI
```

Results will be displayed in console.
