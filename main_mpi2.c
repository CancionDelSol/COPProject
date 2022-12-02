#include "utils.h"
#include "fourier-mpi.h"
#include <math.h>
#include <stdlib.h>
// #include <iostream>
#include <stdbool.h> // added
#include <mpi.h>

#define ARRAY_LENGTH 1000

// Helper function to compare doubles
bool approximatelyEqual(double a, double b, double epsilon)
{
    return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

void compareMax(double *signalArray, int sampleSize, int expectedIndex, double expectedValue, double epsilon)
{
    // Initialize max variables
    int greatestValIndex = 0;
    double greatestValue = -1000;

    // Get max variables from data
    for (int i = 0; i < sampleSize; i++)
    {
        if (signalArray[i] > greatestValue)
        {
            greatestValIndex = i;
            greatestValue = signalArray[i];
        }
    }

    // Report expected results from test
    if (greatestValIndex == expectedIndex)
    {
        // std::cout << "PASS - ";
        printf("PASS - ");
    }
    else
    {
        // std::cout << "FAIL - ";
        printf("FAIL - ");
    }
    // std::cout << "Index of Greatest Value (expected, actual): " << expectedIndex << ", " << greatestValIndex << '\n';

    printf("Index of Greatest Value (expected, actual): %d, %d\n", expectedIndex, greatestValIndex);

    if (approximatelyEqual(expectedValue, greatestValue, epsilon))
    {
        // std::cout << "PASS - ";
        printf("PASS - ");
    }
    else
    {
        // std::cout << "FAIL - ";
        printf("FAIL - ");
    }
    // std::cout << "Value of Greatest Value (expected, actual): " << expectedValue << ", " << greatestValue << '\n'<< '\n';
    printf("Value of Greatest Value (expected, actual): %lf, %lf\n ", expectedValue, greatestValue);
}

// Test #1 - Generate a simple input signal of small sample size (4)
void testSmallSampleSize(int my_rank, int comm_sz, MPI_Comm comm)
{
    int sampleSize = 4;
    int l_sz = sampleSize / comm_sz;
    int l_idx_s = my_rank * l_sz;
    int l_idx_e = l_idx_s + l_sz;

    printf("\nProcess: %d ===> start index: %d end index: %d local size: %d", my_rank, l_idx_s, l_idx_e, l_sz);

    double *inputSignal = NULL;
    double *output = NULL;

    inputSignal = malloc(sizeof(double) * sampleSize);
    output = malloc(sizeof(double) * sampleSize);
    for (int i = 0; i < sampleSize; i++)
    {
        output[i] = 0;
    }

    inputSignal[0] = 1;
    inputSignal[1] = 4;
    inputSignal[2] = 9;
    inputSignal[3] = 16;

    // MPI_Bcast(inputSignal, sampleSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(output, sampleSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Expected test results
    int EXPECTED_INDEX = 0;
    double EXPECTED_VALUE = 30;
    double EPSILON = 0.1;

    MPI_Gather(inputSignal, ARRAY_LENGTH, MPI_DOUBLE,
               output, ARRAY_LENGTH, MPI_DOUBLE, 0, comm);
    printf("\nTest #1 - testSmallSampleSize");

    GetFourierTransform(inputSignal, sampleSize, output, l_idx_s, l_idx_e, my_rank);
    printf("\nProccess %d ===> Printing the Out Array Elements: ", my_rank);

    for (int i = 0; i < sampleSize; i++)
        printf("%lf ", output[i]);

    // if (my_rank == 0)
    compareMax(output, sampleSize, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    MPI_Barrier(comm);
    // free(output);
    // free(inputSignal);
}

// Test #2 - Generate a sine wave of 4 periods with sample size 1000
void testSineWave(int my_rank, int comm_sz, MPI_Comm comm)
{
    int l_sz = ARRAY_LENGTH / comm_sz;
    int l_idx_s = my_rank * l_sz;
    int l_idx_e = l_idx_s + l_sz;

    printf("\nTEST-2  Process: %d ===> start index: %d end index: %d local size: %d", my_rank, l_idx_s, l_idx_e, l_sz);

    // Generate test sine wave input signal
    // double *sinWave = (double *)malloc(sizeof(double) * ARRAY_LENGTH);
    double *inputSignal = malloc(sizeof(double) * ARRAY_LENGTH);
    double *output = malloc(sizeof(double) * ARRAY_LENGTH);

    // Expected test results
    int EXPECTED_INDEX = 996;
    double EXPECTED_VALUE = 0.000205102;
    double EPSILON = 0.001;

    // Total of four periods divided into ARRAY_LENGTH pieces
    double stepSize = (4.0 * TWOPI) / ARRAY_LENGTH;

    for (int i = 0; i < ARRAY_LENGTH; i++)
    {
        // sinWave[i] = sin(stepSize * i);
        inputSignal[i] = sin(stepSize * i);
    }

    // std::cout << "Test #2 - testSineWave" << '\n';
    printf("\nTest #2 - testSineWave");

    if (my_rank == 0)
        MPI_Gather(inputSignal, ARRAY_LENGTH, MPI_DOUBLE,
                   output, ARRAY_LENGTH, MPI_DOUBLE, 0, comm);

    GetFourierTransform(inputSignal, ARRAY_LENGTH, output, l_idx_s, l_idx_e, my_rank);

    compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    // MPI_Barrier(comm);
    // free(inputSignal);
    // free(output);
}

// Test #3 - Generate an arbitrary continuous signal with sample size 1000
void testArbitraryContinuousSignal(int my_rank, int comm_sz, MPI_Comm comm)
{
    int l_sz = ARRAY_LENGTH / comm_sz;
    int l_idx_s = my_rank * l_sz;
    int l_idx_e = l_idx_s + l_sz;

    printf("\nTEST-3  Process: %d ===> start index: %d end index: %d local size: %d", my_rank, l_idx_s, l_idx_e, l_sz);

    // input_signal(t) = 5 + 2cos(2*PI*t - 90 deg) + 3cos(4*PIT*t)
    // input_signal(t) components = DC + 1Hz + 2Hz
    // sample rate = 4 Hz
    double *inputSignal = malloc(sizeof(double) * ARRAY_LENGTH);
    double *output = malloc(sizeof(double) * ARRAY_LENGTH);

    // Expected test results
    int EXPECTED_INDEX = 0;
    double EXPECTED_VALUE = 5000;
    double EPSILON = 0.1;

    // Generate input signal
    for (int i = 0; i < ARRAY_LENGTH; i++)
    {
        inputSignal[i] = 5 + 2 * cos((PI / 2) * i - (PI / 2)) + 3 * cos(PI * i);
    }

    // if (my_rank == 0)
    MPI_Gather(inputSignal, ARRAY_LENGTH, MPI_DOUBLE,
               output, ARRAY_LENGTH, MPI_DOUBLE, 0, comm);

    // std::cout << "Test #3 - testArbitraryContinuousSignal" << '\n';
    printf("\nTest #3 - testArbitraryContinuousSignal");
    GetFourierTransform(inputSignal, ARRAY_LENGTH, output, l_idx_s, l_idx_e, my_rank);

    if (my_rank == 0)
    {
        compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);
    }

    MPI_Barrier(comm);
    // free(output);
    // free(inputSignal);
}

int main(int argc, char **argv)
{
    int my_rank, comm_sz;
    /* Let the system do what it needs to start up MPI */
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    printf("%d\n", comm_sz);
    MPI_Comm comm = MPI_COMM_WORLD;
    printf("\n Process %d is running", my_rank);

    testSmallSampleSize(my_rank, comm_sz, comm);
    testSineWave(my_rank, comm_sz, comm);
    testArbitraryContinuousSignal(my_rank, comm_sz, comm);

    // MPI_Barrier(comm);

    MPI_Finalize();
}
