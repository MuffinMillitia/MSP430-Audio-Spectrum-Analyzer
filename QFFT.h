#define GLOBAL_Q 12
#include "QmathLib.h"
/*
 * QFFT.h
 * Fixed-Point Discrete Fourier transform using Cooley-Turkey FFT algorithm
 *  Created on: Nov 16, 2020
 *      Author: Josef Matuska
 */

#ifndef QFFT_H_
#define QFFT_H_

// represents complex number
struct qcomplex
{
    _q re;
    _q im;
};

// calculates Euler's formula
struct qcomplex qeulers(_q x);

// multiplies two complex numbers together
struct qcomplex qc_multiply(struct qcomplex a, struct qcomplex b);

// calculates the fast Fourier transform recursively
void qFFT_algorithm(struct qcomplex * euler_lookup, int size, struct qcomplex * y, _q * x, int n, int s);

// populates the Hamming window lookup table
void populate_hamming_lookup(_q * table, int size);

// populates Euler's coefficient lookup table
void populate_euler_lookup(struct qcomplex * table, int size);

#endif /* QFFT_H_ */
