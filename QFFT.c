#include <msp430.h>
#include "QFFT.h"

#define PI 3.14159265

/*
 * QFFT.c
 * Fixed Point Discrete Fourier transform using Cooley-Turkey FFT algorithm
 *  Created on: Nov 16, 2020
 *      Author: Josef Matuska
 */

struct qcomplex qeulers(_q x)
{
    struct qcomplex result;

    result.re = _Qcos(x);
    result.im = _Qsin(x);

    return result;
}


struct qcomplex qc_multiply(struct qcomplex a, struct qcomplex b)
{
    struct qcomplex result;

    // multiply two complex numbers
    result.re = _Qmpy(a.re, b.re) - _Qmpy(a.im, b.im);
    result.im = _Qmpy(a.re, b.im) + _Qmpy(a.im, b.re);

    return result;
}

void qFFT_algorithm(struct qcomplex * euler_lookup, int size, struct qcomplex * y, _q * x, int n, int s)
{
    if (n == 1)
    {
        y[0].re = x[0];                                 // FFT trivial base case for size 1
    }
    else
    {
        qFFT_algorithm(euler_lookup, size, y, x, n/2, 2*s);                  // Traverses left branch of tree
        qFFT_algorithm(euler_lookup, size, y+n/2, x+s, n/2, 2*s);            // Traverses right branch of tree

        int k;                                          // index
        struct qcomplex t;                              // previous value of y

        struct qcomplex c1;                             // struct qcomplex constant to reduce computations when adding struct qcomplex numbers
        struct qcomplex c2;

        for (k = 0; k < n/2; k++)
        {
            t = y[k];

            //c1 = qeulers(_Q(-2 * PI * k / n));        // method without lookup tables
            c1 = *(euler_lookup+size*k+n);              // method with lookup tables
            c2 = qc_multiply(c1, y[k+n/2]);


            // Update transformed values recursively
            y[k].re = _Qdiv2(t.re + c2.re);
            y[k].im = _Qdiv2(t.im + c2.im);

            y[k+n/2].re = _Qdiv2(t.re - c2.re);
            y[k+n/2].im = _Qdiv2(t.im - c2.im);
        }
    }
}

void populate_hamming_lookup(_q * table, int size)
{
    int i;
    for(i = 0; i < size; i++)
    {
        // calculate Hamming window
        table[i] = _Q(.54) - _Qmpy(_Q(.46),_Qcos(_Q((2*3.14159*i)/size)));
    }
}

void populate_euler_lookup(struct qcomplex * table, int size)
{
    int k;
    for(k = 0; k < size/2; k++)
    {
        int n;
        for(n = 1; n <= size; n++)
        {
            // Calculate Euler's coefficients for every possible value used in FFT
            (table+size*k+n)->re = _Qcos(_Q(-2 * PI * k / n));
            (table+size*k+n)->im = _Qsin(_Q(-2 * PI * k / n));
        }
    }
}

