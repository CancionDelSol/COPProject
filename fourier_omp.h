#ifndef FOURIER_H
#define FOURIER_H

#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "microtime.h"
#include <omp.h>

/*
 * General Notes on fourier analysis
 * I'm assuming for this project, we
 * can just do a discrete version of
 * the fourier transform and not
 * worry about parallelizing the FFT
 * algorithm, that would take forever
 * Here is the pseudocode though to get the fourier transform at index n:
 *  float GetFourierAtF_n(float* baseFunction)
 *  {
 *    float sum = 0.0;
 *    for (int k = 0; k < N; k++) { // function f(k) over range N
 *      float f_kAtkEqualsN = GetBaseFunctionValueAtXEqualsk(k);
 *      sum += f_kAtkEqualsN * WaveFunctionAtk(k);
 *    }
 *  }
 *  float* GetFourierTransform(float* baseFunction, float* outputArray, int length) {
 *    for (int n = 0; n < length; n++) {
 *      otuputArray[n] = GetFourierAtF_n(baseFunction);
 *    }
 *  }
 */

// We are only accounting for non-imaginary components
// This means we can break down the exponential into it's
// cosine and sine components and then toss out the sin portion
double FourierInnerSum(double baseFunctionAtK, /* The original function evaluated at index k */
                       int N,                  /* The total length of both arrays */
                       int k,                  /* The index k from the original function */
                       int n                   /* The index n of the output function */
)
{

    double termOne = cos((TWOPI / N) * k * n);
    return baseFunctionAtK * termOne;
}

void GetFourierTransform(
    double *baseFunction,        /* The array of values for the original function. */
    int baseFunctionArrayLength, /* The length of the array */
    double *outputArray          /* The output array for the calculated fourier transform */
)
{
    int N = baseFunctionArrayLength;
    int n,k;
    #pragma omp parallel for schedule(dynamic,32) private(n, k) shared(outputArray)
    for (n = 0; n < N; n++) {
        double sum = 0.0;
        for (k = 0; k <N ; k++) {
            double baseFunctionAtK = baseFunction[k];
            sum += FourierInnerSum(baseFunctionAtK, N, k, n);
        }
        outputArray[n] = sum;
    }
}
#endif