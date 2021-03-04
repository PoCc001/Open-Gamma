/**
* Copyright Johannes Kloimböck 2021.
* Distributed under the Boost Software License, Version 1.0.
* (See accompanying file LICENSE or copy at
* https://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef FASTGAMMA_H
#define FASTGAMMA_H

#include<math.h>
#include<inttypes.h>

// FOR THE LOGGAMMA FUNCTIONS TO WORK, DOUBLE AND FLOAT ARE
// EXPECTED TO BE 64 BITS AND 32 BITS WIDE RESPECTIVELY!

/* 
 * useful constants
 * they don't need many digits,
 * as the approximation functions
 * won't be able to compute gamma and
 * loggamma to such a high precision
 * anyways
 */
#define TWO_PI_D 6.283185
#define TWO_PI_F (float)(TWO_PI_D)
#define LOG_2_D 0.693147
#define LOG_2_F (float)(LOG_2_D)
#define RECIP_E_D 0.367879
#define RECIP_E_F (float)(RECIP_E_D)

#define DO_INLINE 1 // change to 0 to prevent inlining the following functions

#if DO_INLINE != 0
#define FUNC_INLINE inline
#else
#define FUNC_INLINE
#endif

#define USE_APPROX_ARRAY 1 // change to 0 to make fast_loggamma even faster but more incorrect

// array to correct the approximation of the logarithm in fast_loggamma
double LD_ARRAY [] = {0.0, 0.04439411, 0.08746284, 0.129283,
	0.169925, 0.2094534, 0.2479275, 0.2854022, 0.3219281, 
	0.357552, 0.3923174, 0.4262647, 0.4594316, 0.4918531,
	0.523562, 0.5545889, 0.5849625, 0.6147098, 0.6438562,
	0.672425, 0.7004397, 0.7279205, 0.7548875, 0.7813597,
	0.8073549, 0.83289, 0.857981, 0.882643, 0.9068906,
	0.9307373, 0.9541963, 0.9772799};

/*
 * Two unions to enable bit manipulation
 * on floating-point arguments, so that a
 * quick, rough estimate of the natural log
 * can be calculated.
 */
union double_to_i64 {
	double d;
	uint64_t i;
} double_i64;

union float_to_i32 {
	float f;
	uint32_t i;
} float_i32;

/*
 * Approximates the gamma function for the double-precision
 * input variable x. The greater the argument is the more
 * precision the result is going to have. The function DOES
 * NOT perform any argument checks! This is done on purpose.
 */
FUNC_INLINE double fast_gamma(double x) {
	x -= 1.0;
	// standard Stirling formula
	return sqrt(TWO_PI_D * x) * pow(x * RECIP_E_D, x);
}

/*
 * Approximates the gamma function for the single-precision
 * input variable x. The greater the argument is the more
 * precision the result is going to have. The function DOES
 * NOT perform any argument checks! This is done on purpose.
 */
FUNC_INLINE float fast_gammaf(float x) {
	x -= 1.0f;
	// standard Stirling formula
	return sqrtf(TWO_PI_F * x) * powf(x * RECIP_E_F, x);
}

/*
 * Approximates the natural logarithm of the gamma function
 * for the double-precision input variable x. The greater
 * the argument is the more precision the result is going
 * to have. The function DOES NOT perform any argument
 * checks! This is done on purpose.
 */
FUNC_INLINE double fast_loggamma(double x) {
	union double_to_i64 ux;
	ux.d = x;
#if USE_APPROX_ARRAY == 0
	ux.i &= (uint64_t)(0x3fffffffffffffffULL);
#else
	ux.i -= (uint64_t)(0x3ff0000000000000ULL);
	unsigned int array_index = (unsigned int)(ux.i >> 47) & 0b11111;
#endif
	ux.i >>= 52;
	double log = (double)(ux.i);
#if USE_APPROX_ARRAY != 0
	log += LD_ARRAY[array_index];
#endif
	log *= LOG_2_D; // a rough estimate of the natural log of x
	return x * (log - 1.0);
}

/*
 * Approximates the natural logarithm of the gamma function
 * for the single-precision input variable x. The greater
 * the argument is the more precision the result is going
 * to have. The function DOES NOT perform any argument
 * checks! This is done on purpose.
 */
FUNC_INLINE float fast_loggammaf(float x) {
	union float_to_i32 ux;
	ux.f = x;
#if USE_APPROX_ARRAY == 0
	ux.i &= (uint32_t)(0x3fffffffUL);
#else
	ux.i -= (uint32_t)(0x3f800000UL);
	unsigned int array_index = (unsigned int)(ux.i >> 18) & 0b11111;
#endif
	ux.i >>= 23;
	float log = (float)(ux.i);
#if USE_APPROX_ARRAY != 0
	log += (float)(LD_ARRAY[array_index]);
#endif
	log *= LOG_2_F; // a rough estimate of the natural log of x
	return x * log;
}

#endif
