#include "utils.h"
#include "fourier.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>

#define ARRAY_LENGTH 1000

/*
 * Main entry
 */
int main(int argc, char** argv) {
    // Here is the example usage
    double* sinWave = (double*)malloc(sizeof(double) * ARRAY_LENGTH);
    double* output = (double*)malloc(sizeof(double) * ARRAY_LENGTH);

    // Total of four periods divided into ARRAY_LENGTH pieces
    double stepSize = (4.0 * TWOPI) / ARRAY_LENGTH;

    for (int i = 0; i < ARRAY_LENGTH; i++) {
        sinWave[i] = sin(stepSize * i);
    }

    GetFourierTransform(sinWave, ARRAY_LENGTH, output);

    int greatestValIndex = 0;
    double greatestValue = -1000;
    for (int i = 0; i < ARRAY_LENGTH; i++) {
        if (output[i] > greatestValue) {
            greatestValIndex = i;
            greatestValue = output[i];
        } 
    }
    std::cout << greatestValIndex << '\n';
    std::cout << greatestValue << '\n';
}