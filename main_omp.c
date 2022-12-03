#include "utils.h"
#include "fourier_omp.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

#define ARRAY_LENGTH 1000
#define CHUNK_SIZE 4

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
void testSmallSampleSize()
{
    double *inputSignal = NULL;
    double *output = NULL;
    int arraySize = 4;

    // Allocate memory for the input and output arrays
    inputSignal = (double*)malloc(sizeof(double) * arraySize);
    output = (double*)malloc(sizeof(double) * arraySize);
    for (int i = 0; i < arraySize; i++)
        output[i] = 0;

    inputSignal[0] = 1;
    inputSignal[1] = 4;
    inputSignal[2] = 9;
    inputSignal[3] = 16;

    // Expected test results
    int EXPECTED_INDEX = 0;
    double EXPECTED_VALUE = 30;
    double EPSILON = 0.001;

    printf("\nTest #1 - testSmallSampleSize\n");

    GetFourierTransform(inputSignal, arraySize, output);

    compareMax(output, arraySize, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    // Free the allocated memory
    free(output);
    free(inputSignal);
}

// Test #2 - Generate a sine wave of 4 periods with sample size 1000
void testSineWave()
{
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

    printf("\nTest #2 - testSineWave\n");

    GetFourierTransform(inputSignal, ARRAY_LENGTH, output);

    compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    // Free the allocated memory
    free(output);
    free(inputSignal);

}

// Test #3 - Generate an arbitrary continuous signal with sample size 1000
void testArbitraryContinuousSignal()
{
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
        inputSignal[i] = 5 + 2 * cos((PI / 2) * i - (PI / 2)) + 3 * cos(PI * i);

    printf("\nTest #3 - testArbitraryContinuousSignal\n");

    GetFourierTransform(inputSignal, ARRAY_LENGTH, output);

    compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);

    free(output);
    free(inputSignal);

}

int main(int argc, char **argv)
{
    // For time outputs
    double time1;
    double time2;
    double t;

    // This processes rank and process count
    int my_rank, processCount;

    // Set our number of threads to the
    //  number of available processors
    omp_set_num_threads(omp_get_num_procs());

    // Small sample size test
    time1 = microtime();
    testSmallSampleSize();
    time2 = microtime();

    t = time2 - time1;
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());
    
    // Sine wave test
    time1 = microtime();
    testSineWave();
    time2 = microtime();

    t = time2 - time1;
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());
    
    // Arbitrary continuous signal test
    time1 = microtime();
    testArbitraryContinuousSignal();
    time2 = microtime();

    t = time2 - time1;
    printf("\nTime = %g us\n", t);
    printf("Timer Resolution = %g us\n", getMicrotimeResolution());

}