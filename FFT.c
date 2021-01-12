#include <msp430.h>
#include "FFT.h"
#include "math.h"

#define PI 3.14159265
/*
 * FFT.c
 * Discrete Fourier transform using Cooley-Turkey FFT algorithm
 *  Created on: Nov 16, 2020
 *      Author: Josef
 */

struct complex eulers(float x)
{
    struct complex result;
    struct complex *pr = &result;
    float *px = &x;

    result.re = cos(x);
    result.im = sin(x);


    return result;
}


struct complex c_multiply(struct complex a, struct complex b)
{
    struct complex result;
    struct complex *pr = &result;
    struct complex *pa = &a;
    struct complex *pb = &b;

    result.re = a.re * b.re - a.im * b.im;
    result.im = a.re * b.im + a.im * b.re;

    return result;
}


void FFT_algorithm(struct complex * y, float * x, int n, int s)
{
    if (n == 1)
    {
        y[0].re = x[0];                         // FFT trivial base case for size 1
    }
    else
    {
        FFT_algorithm(y, x, n/2, 2*s);                    // Traverses left branch of tree
        FFT_algorithm(y+n/2, x+s, n/2, 2*s);            // Traverses right branch of tree

        int k;                                          // index
        struct complex t;                              // previous value of y
        struct complex *pt = &t;

        struct complex c1;                             // struct complex constant to reduce computations when adding struct complex numbers
        struct complex c2;

        struct complex *p1 = &c1;
        struct complex *p2 = &c2;


        for (k = 0; k < n/2; k++)
        {
            t = y[k];

            c1 = eulers(-2 * PI * k / n);
            c2 = c_multiply(c1, y[k+n/2]);


            y[k].re = (t.re + c2.re)/2;
            y[k].im = (t.im + c2.im)/2;

            y[k+n/2].re = (t.re - c2.re)/2;
            y[k+n/2].im = (t.im - c2.im)/2;
        }
    }
}


