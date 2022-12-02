#include "utils.h"
#include "fourier_mpi.h"
#include <math.h>
#include <stdlib.h>
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
        printf("PASS - ");
    }
    else
    {
        printf("FAIL - ");
    }

    printf("Index of Greatest Value (expected, actual): %d, %d\n", expectedIndex, greatestValIndex);

    if (approximatelyEqual(expectedValue, greatestValue, epsilon))
    {
        printf("PASS - ");
    }
    else
    {
        printf("FAIL - ");
    }
    printf("Value of Greatest Value (expected, actual): %lf, %lf\n ", expectedValue, greatestValue);
}

// Test #1 - Generate a simple input signal of small sample size (4)
void testSmallSampleSize(int my_rank, int processCount, MPI_Comm comm)
{
    int sampleSize = 4;
    int l_sz = sampleSize / (processCount - 1);
    int l_idx_s = (my_rank-1) * l_sz;
    int l_idx_e = l_idx_s + l_sz;

    double *inputSignal = NULL;
    double *output = NULL;

    // Allocate memory for the input and output arrays
    inputSignal = (double*)malloc(sizeof(double) * sampleSize);
    output = (double*)malloc(sizeof(double) * sampleSize);
    for (int i = 0; i < sampleSize; i++)
        output[i] = 0;

    inputSignal[0] = 1;
    inputSignal[1] = 4;
    inputSignal[2] = 9;
    inputSignal[3] = 16;

    // Expected test results
    int EXPECTED_INDEX = 0;
    double EXPECTED_VALUE = 30;
    double EPSILON = 0.001;

    if (my_rank == 0)
        printf("\nTest #1 - testSmallSampleSize\n");

    if (my_rank == 0) {
        // Gather the results from the other processes
        for (int curProcess = 1; curProcess < processCount; curProcess++) {
            double response[l_sz] = { 0 };
            MPI_Status status;
            MPI_Recv(&response, l_sz, MPI_DOUBLE, curProcess, 1, comm, &status);

            int startIndex = (curProcess - 1) * l_sz;
            for (int i = 0; i < l_sz; i++) {
                output[startIndex + i] = response[i];
            }
        }
    } else {
        // Get the section of the output array for this process
        GetFourierTransform(inputSignal, sampleSize, output, l_idx_s, l_idx_e, my_rank);

        // Send results back to process 0
        double msg[l_sz] = { 0 };
        for (int i = 0; i < l_sz; i++) {
            msg[i] = output[l_idx_s + i];
        }
        MPI_Send(&msg, l_sz, MPI_DOUBLE, 0, 1, comm);
    }

    if (my_rank == 0)
        compareMax(output, sampleSize, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    // Free the allocated memory
    free(output);
    free(inputSignal);
}

// Test #2 - Generate a sine wave of 4 periods with sample size 1000
void testSineWave(int my_rank, int processCount, MPI_Comm comm)
{
    int l_sz = ARRAY_LENGTH / (processCount - 1);
    int l_idx_s = (my_rank-1) * l_sz;
    int l_idx_e = l_idx_s + l_sz;

    // Generate test sine wave input signal
    double *inputSignal = (double*)malloc(sizeof(double) * ARRAY_LENGTH);
    double *output = (double*)malloc(sizeof(double) * ARRAY_LENGTH);

    // Expected test results
    int EXPECTED_INDEX = 996;
    double EXPECTED_VALUE = 0.000205102;
    double EPSILON = 0.001;

    // Total of four periods divided into ARRAY_LENGTH pieces
    double stepSize = (4.0 * TWOPI) / ARRAY_LENGTH;

    for (int i = 0; i < ARRAY_LENGTH; i++)
        inputSignal[i] = sin(stepSize * i);

    if (my_rank == 0)
        printf("\nTest #2 - testSineWave\n");

    if (my_rank == 0) {
        // Gather the results from the other processes
        for (int curProcess = 1; curProcess < processCount; curProcess++) {
            double response[l_sz] = { 0 };
            MPI_Status status;
            MPI_Recv(&response, l_sz, MPI_DOUBLE, curProcess, 1, comm, &status);

            int startIndex = (curProcess - 1) * l_sz;
            for (int i = 0; i < l_sz; i++) {
                output[startIndex + i] = response[i];
            }
        }
    }else {
        // Get the section of the output array for this process
        GetFourierTransform(inputSignal, ARRAY_LENGTH, output, l_idx_s, l_idx_e, my_rank);

        // Send results back to process 0
        double msg[l_sz] = { 0 };
        for (int i = 0; i < l_sz; i++) {
            msg[i] = output[l_idx_s + i];
        }
        MPI_Send(&msg, l_sz, MPI_DOUBLE, 0, 1, comm);
    }

    if (my_rank == 0)
        compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    // Free the allocated memory
    free(output);
    free(inputSignal);

}

// Test #3 - Generate an arbitrary continuous signal with sample size 1000
void testArbitraryContinuousSignal(int my_rank, int processCount, MPI_Comm comm)
{
    int l_sz = ARRAY_LENGTH / processCount;
    int l_idx_s = my_rank * l_sz;
    int l_idx_e = l_idx_s + l_sz;

    printf("\nTEST-3  Process: %d ===> start index: %d end index: %d local size: %d", my_rank, l_idx_s, l_idx_e, l_sz);

    // input_signal(t) = 5 + 2cos(2*PI*t - 90 deg) + 3cos(4*PIT*t)
    // input_signal(t) components = DC + 1Hz + 2Hz
    // sample rate = 4 Hz
    double *inputSignal = (double*)malloc(sizeof(double) * ARRAY_LENGTH);
    double *output = (double*)malloc(sizeof(double) * ARRAY_LENGTH);

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
    // This processes rank and process count
    int my_rank, processCount;

    // Initialize MPI framework
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    // Use world scope communication
    MPI_Comm comm = MPI_COMM_WORLD;

    // Small sample size test
    MPI_Barrier(comm);
    if (my_rank == 0) {
        double time1 = microtime();
        testSmallSampleSize(my_rank, processCount, comm);
        double time2 = microtime();

        double t = time2 - time1;
        printf("\nTime = %g us\n", t);
        printf("Timer Resolution = %g us\n", getMicrotimeResolution());
    } else {
        testSmallSampleSize(my_rank, processCount, comm);
    }
    
    // Sine wave test
    MPI_Barrier(comm);
    if (my_rank == 0) {
        double time1 = microtime();
        testSineWave(my_rank, processCount, comm);
        double time2 = microtime();

        double t = time2 - time1;
        printf("\nTime = %g us\n", t);
        printf("Timer Resolution = %g us\n", getMicrotimeResolution());
    } else {
        testSineWave(my_rank, processCount, comm);
    }
    
    /*
    // Arbitrary continuous signal test
    MPI_Barrier(comm);
    if (my_rank == 0) {
        double time1 = microtime();
        testArbitraryContinuousSignal(my_rank, processCount, comm);
        double time2 = microtime();

        double t = time2 - time1;
        printf("\nTime = %g us\n", t);
        printf("Timer Resolution = %g us\n", getMicrotimeResolution());
    } else {
        testArbitraryContinuousSignal(my_rank, processCount, comm);
    }
    */

    // Pause before finalizing
    MPI_Barrier(comm);

    // Finalize
    MPI_Finalize();
}