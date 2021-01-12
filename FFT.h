//#include "msp430_math.h"

/*
 * FFT.h
 * Discrete Fourier transform using Cooley-Turkey FFT algorithm
 *  Created on: Nov 16, 2020
 *      Author: Josef Matuska
 */

#ifndef FFT_H_
#define FFT_H_

struct complex
{
    float re;
    float im;
};

struct complex eulers(float x);

struct complex c_multiply(struct complex a, struct complex b);

void FFT_algorithm(struct complex * y, float * x, int n, int s);


#endif /* FFT_H_ */
