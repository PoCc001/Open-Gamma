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

#define USE_APPROX_ARRAYS 1 // change to 0 to make fast_loggamma even faster but more incorrect

// arrays to correct the approximation of the logarithm in fast_loggamma
#if USE_APPROX_ARRAYS != 0
double LD_ARRAY_D [128] = {
0.000000, 0.011227, 0.022368, 0.033423, 0.044394,
0.055282, 0.066089, 0.076816, 0.087463, 0.098032,
0.108524, 0.118941, 0.129283, 0.139551, 0.149747,
0.159871, 0.169925, 0.179909, 0.189825, 0.199672,
0.209453, 0.219169, 0.228819, 0.238405, 0.247928,
0.257388, 0.266787, 0.276124, 0.285402, 0.294621,
0.303781, 0.312883, 0.321928, 0.330917, 0.339850,
0.348728, 0.357552, 0.366322, 0.375039, 0.383704,
0.392317, 0.400879, 0.409391, 0.417853, 0.426265,
0.434628, 0.442943, 0.451211, 0.459432, 0.467606,
0.475733, 0.483816, 0.491853, 0.499846, 0.507795,
0.515700, 0.523562, 0.531381, 0.539159, 0.546894,
0.554589, 0.562242, 0.569856, 0.577429, 0.584962,
0.592457, 0.599913, 0.607330, 0.614710, 0.622052,
0.629357, 0.636625, 0.643856, 0.651052, 0.658211,
0.665336, 0.672425, 0.679480, 0.686501, 0.693487,
0.700440, 0.707359, 0.714245, 0.721099, 0.727920,
0.734710, 0.741467, 0.748193, 0.754888, 0.761551,
0.768184, 0.774787, 0.781360, 0.787903, 0.794416,
0.800900, 0.807355, 0.813781, 0.820179, 0.826549,
0.832890, 0.839204, 0.845490, 0.851749, 0.857981,
0.864186, 0.870365, 0.876517, 0.882643, 0.888743,
0.894818, 0.900867, 0.906891, 0.912889, 0.918863,
0.924812, 0.930737, 0.936638, 0.942514, 0.948367,
0.954196, 0.960002, 0.965784, 0.971544, 0.977280,
0.982994, 0.988685, 0.994353
};


float LD_ARRAY_F [128] = {
0.000000f, 0.011227f, 0.022368f, 0.033423f, 0.044394f,
0.055282f, 0.066089f, 0.076816f, 0.087463f, 0.098032f,
0.108524f, 0.118941f, 0.129283f, 0.139551f, 0.149747f,
0.159871f, 0.169925f, 0.179909f, 0.189825f, 0.199672f,
0.209453f, 0.219169f, 0.228819f, 0.238405f, 0.247928f,
0.257388f, 0.266787f, 0.276124f, 0.285402f, 0.294621f,
0.303781f, 0.312883f, 0.321928f, 0.330917f, 0.339850f,
0.348728f, 0.357552f, 0.366322f, 0.375039f, 0.383704f,
0.392317f, 0.400879f, 0.409391f, 0.417853f, 0.426265f,
0.434628f, 0.442943f, 0.451211f, 0.459432f, 0.467606f,
0.475733f, 0.483816f, 0.491853f, 0.499846f, 0.507795f,
0.515700f, 0.523562f, 0.531381f, 0.539159f, 0.546894f,
0.554589f, 0.562242f, 0.569856f, 0.577429f, 0.584962f,
0.592457f, 0.599913f, 0.607330f, 0.614710f, 0.622052f,
0.629357f, 0.636625f, 0.643856f, 0.651052f, 0.658211f,
0.665336f, 0.672425f, 0.679480f, 0.686501f, 0.693487f,
0.700440f, 0.707359f, 0.714245f, 0.721099f, 0.727920f,
0.734710f, 0.741467f, 0.748193f, 0.754888f, 0.761551f,
0.768184f, 0.774787f, 0.781360f, 0.787903f, 0.794416f,
0.800900f, 0.807355f, 0.813781f, 0.820179f, 0.826549f,
0.832890f, 0.839204f, 0.845490f, 0.851749f, 0.857981f,
0.864186f, 0.870365f, 0.876517f, 0.882643f, 0.888743f,
0.894818f, 0.900867f, 0.906891f, 0.912889f, 0.918863f,
0.924812f, 0.930737f, 0.936638f, 0.942514f, 0.948367f,
0.954196f, 0.960002f, 0.965784f, 0.971544f, 0.977280f,
0.982994f, 0.988685f, 0.994353f
};

#endif

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
#if USE_APPROX_ARRAYS == 0
	ux.i &= (uint64_t)(0x3fffffffffffffffULL);
#else
	ux.i -= (uint64_t)(0x3ff0000000000000ULL);
	unsigned int array_index = (unsigned int)(ux.i >> 45) & 0b1111111;
#endif
	ux.i >>= 52;
	double log = (double)(ux.i);
#if USE_APPROX_ARRAYS != 0
	log += LD_ARRAYS_D[array_index];
#endif
	log *= LOG_2_D; // a rough estimate of the natural log of x
#if USE_APPROX_ARRAYS == 0
	return x * log;
#else
	return x * (log - 1.0);
#endif
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
#if USE_APPROX_ARRAYS == 0
	ux.i &= (uint32_t)(0x3fffffffUL);
#else
	ux.i -= (uint32_t)(0x3f800000UL);
	unsigned int array_index = (unsigned int)(ux.i >> 16) & 0b1111111;
#endif
	ux.i >>= 23;
	float log = (float)(ux.i);
#if USE_APPROX_ARRAYS != 0
	log += LD_ARRAY_F[array_index]);
#endif
	log *= LOG_2_F; // a rough estimate of the natural log of x
#if USE_APPROX_ARRAYS == 0
	return x * log;
#else
	return x * (log - 1.0f);
#endif
}

#endif
