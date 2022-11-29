#include "utils.h"
#include "fourier.h"
#include <math.h>
#include <stdlib.h>
// #include <iostream>
#include <stdbool.h> // added

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
void testSmallSampleSize()
{
    int sampleSize = 4;

    double *inputSignal = (double *)malloc(sizeof(double) * sampleSize);
    double *output = (double *)malloc(sizeof(double) * sampleSize);

    // Expected test results
    int EXPECTED_INDEX = 0;
    double EXPECTED_VALUE = 30;
    double EPSILON = 0.1;

    // Generate input signal
    inputSignal[0] = 1;
    inputSignal[1] = 4;
    inputSignal[2] = 9;
    inputSignal[3] = 16;

    printf("\nTest #1 - testSmallSampleSize");
    GetFourierTransform(inputSignal, sampleSize, output);
    compareMax(output, sampleSize, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);
}

// Test #2 - Generate a sine wave of 4 periods with sample size 1000
void testSineWave()
{
    // Generate test sine wave input signal
    double *sinWave = (double *)malloc(sizeof(double) * ARRAY_LENGTH);
    double *output = (double *)malloc(sizeof(double) * ARRAY_LENGTH);

    // Expected test results
    int EXPECTED_INDEX = 996;
    double EXPECTED_VALUE = 0.000205102;
    double EPSILON = 0.001;

    // Total of four periods divided into ARRAY_LENGTH pieces
    double stepSize = (4.0 * TWOPI) / ARRAY_LENGTH;

    for (int i = 0; i < ARRAY_LENGTH; i++)
    {
        sinWave[i] = sin(stepSize * i);
    }

    // std::cout << "Test #2 - testSineWave" << '\n';
    printf("\nTest #2 - testSineWave");
    GetFourierTransform(sinWave, ARRAY_LENGTH, output);
    compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);
}

// Test #3 - Generate an arbitrary continuous signal with sample size 1000
void testArbitraryContinuousSignal()
{
    // input_signal(t) = 5 + 2cos(2*PI*t - 90 deg) + 3cos(4*PIT*t)
    // input_signal(t) components = DC + 1Hz + 2Hz
    // sample rate = 4 Hz
    double *inputSignal = (double *)malloc(sizeof(double) * ARRAY_LENGTH);
    double *output = (double *)malloc(sizeof(double) * ARRAY_LENGTH);

    // Expected test results
    int EXPECTED_INDEX = 0;
    double EXPECTED_VALUE = 5000;
    double EPSILON = 0.1;

    // Generate input signal
    for (int i = 0; i < ARRAY_LENGTH; i++)
    {
        inputSignal[i] = 5 + 2 * cos((PI / 2) * i - (PI / 2)) + 3 * cos(PI * i);
    }

    // std::cout << "Test #3 - testArbitraryContinuousSignal" << '\n';
    printf("\nTest #3 - testArbitraryContinuousSignal");
    GetFourierTransform(inputSignal, ARRAY_LENGTH, output);
    compareMax(output, ARRAY_LENGTH, EXPECTED_INDEX, EXPECTED_VALUE, EPSILON);
}

int main(int argc, char **argv)
{
    testSmallSampleSize();
    testSineWave();
    testArbitraryContinuousSignal();
}
