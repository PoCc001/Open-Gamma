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

/*
 * Two unions to enable bit manipulation
 * on floating-point arguments, so that a
 * quick, rough estimate of the natural log
 * can be calculated.
 */
union double_to_ui64 {
	double d;
	uint64_t ui;
} double_ui64;

union float_to_ui32 {
	float f;
	uint32_t ui;
} float_u32;

/*
 * Approximates the gamma function for the double-precision
 * input variable x. The greater the argument is the more
 * precision the result is going to have. The function DOES
 * NOT perform any argument checks! This is done on purpose.
 */
double fast_gamma(double x) {
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
float fast_gammaf(float x) {
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
double fast_loggamma(double x) {
	union double_to_ui64 ux;
	ux.d = x;
	ux.ui -= 0x3ff0000000000000ULL;
	ux.ui >>= 52;
	double log = (double)(ux.ui);
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
float fast_loggammaf(float x) {
	union float_to_ui32 ux;
	ux.f = x;
	ux.ui = 0x3f800000UL;
	ux.ui >>= 23;
	float log = (float)(ux.ui);
	log *= LOG_2_F; // a rough estimate of the natural log of x
	return x * (log - 1.0f);
}

#endif
