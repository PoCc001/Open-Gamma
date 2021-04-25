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
#define LOG_2_D 0.6931472
#define LOG_2_F (float)(LOG_2_D)
#define RECIP_E_D 0.3678794
#define RECIP_E_F (float)(RECIP_E_D)

#define DO_INLINE 1

#if DO_INLINE != 0
#define FUNC_INLINE inline
#else
#define FUNC_INLINE
#endif

#define USE_APPROX_ARRAYS 1

float LD_ARRAY[256] = {
0.0000000f, 0.0056245f, 0.0112273f, 0.0168083f, 0.0223678f, 
0.0279060f, 0.0334230f, 0.0389190f, 0.0443941f, 0.0498485f, 
0.0552824f, 0.0606959f, 0.0660892f, 0.0714624f, 0.0768156f, 
0.0821490f, 0.0874628f, 0.0927571f, 0.0980321f, 0.1032878f, 
0.1085245f, 0.1137422f, 0.1189411f, 0.1241213f, 0.1292830f, 
0.1344263f, 0.1395514f, 0.1446582f, 0.1497471f, 0.1548181f, 
0.1598713f, 0.1649069f, 0.1699250f, 0.1749257f, 0.1799091f, 
0.1848753f, 0.1898246f, 0.1947569f, 0.1996723f, 0.2045711f, 
0.2094534f, 0.2143191f, 0.2191685f, 0.2240017f, 0.2288187f, 
0.2336197f, 0.2384047f, 0.2431740f, 0.2479275f, 0.2526654f, 
0.2573878f, 0.2620948f, 0.2667865f, 0.2714630f, 0.2761244f, 
0.2807708f, 0.2854022f, 0.2900188f, 0.2946207f, 0.2992080f, 
0.3037807f, 0.3083390f, 0.3128830f, 0.3174126f, 0.3219281f, 
0.3264295f, 0.3309169f, 0.3353904f, 0.3398500f, 0.3442959f, 
0.3487282f, 0.3531468f, 0.3575520f, 0.3619438f, 0.3663222f, 
0.3706874f, 0.3750394f, 0.3793784f, 0.3837043f, 0.3880173f, 
0.3923174f, 0.3966048f, 0.4008794f, 0.4051415f, 0.4093909f, 
0.4136279f, 0.4178525f, 0.4220648f, 0.4262648f, 0.4304526f, 
0.4346282f, 0.4387919f, 0.4429435f, 0.4470832f, 0.4512111f, 
0.4553272f, 0.4594316f, 0.4635244f, 0.4676056f, 0.4716752f, 
0.4757334f, 0.4797803f, 0.4838158f, 0.4878400f, 0.4918531f, 
0.4958550f, 0.4998459f, 0.5038257f, 0.5077946f, 0.5117527f, 
0.5156998f, 0.5196363f, 0.5235620f, 0.5274770f, 0.5313815f, 
0.5352754f, 0.5391588f, 0.5430318f, 0.5468945f, 0.5507468f, 
0.5545889f, 0.5584207f, 0.5622424f, 0.5660540f, 0.5698556f, 
0.5736472f, 0.5774288f, 0.5812006f, 0.5849625f, 0.5887146f, 
0.5924570f, 0.5961898f, 0.5999128f, 0.6036263f, 0.6073303f, 
0.6110248f, 0.6147098f, 0.6183855f, 0.6220518f, 0.6257088f, 
0.6293566f, 0.6329952f, 0.6366246f, 0.6402449f, 0.6438562f, 
0.6474584f, 0.6510517f, 0.6546360f, 0.6582115f, 0.6617781f, 
0.6653359f, 0.6688850f, 0.6724253f, 0.6759570f, 0.6794801f, 
0.6829946f, 0.6865005f, 0.6899980f, 0.6934870f, 0.6969675f, 
0.7004397f, 0.7039036f, 0.7073591f, 0.7108064f, 0.7142455f, 
0.7176764f, 0.7210992f, 0.7245139f, 0.7279205f, 0.7313190f, 
0.7347096f, 0.7380923f, 0.7414670f, 0.7448338f, 0.7481928f, 
0.7515441f, 0.7548875f, 0.7582232f, 0.7615512f, 0.7648716f, 
0.7681843f, 0.7714895f, 0.7747871f, 0.7780771f, 0.7813597f, 
0.7846348f, 0.7879026f, 0.7911629f, 0.7944159f, 0.7976615f, 
0.8008999f, 0.8041310f, 0.8073549f, 0.8105716f, 0.8137812f, 
0.8169836f, 0.8201790f, 0.8233672f, 0.8265485f, 0.8297227f, 
0.8328900f, 0.8360504f, 0.8392038f, 0.8423503f, 0.8454901f, 
0.8486229f, 0.8517490f, 0.8548684f, 0.8579810f, 0.8610869f, 
0.8641861f, 0.8672787f, 0.8703647f, 0.8734441f, 0.8765169f, 
0.8795832f, 0.8826430f, 0.8856964f, 0.8887432f, 0.8917837f, 
0.8948178f, 0.8978455f, 0.9008668f, 0.9038818f, 0.9068906f, 
0.9098931f, 0.9128893f, 0.9158794f, 0.9188632f, 0.9218409f, 
0.9248125f, 0.9277780f, 0.9307373f, 0.9336907f, 0.9366379f, 
0.9395792f, 0.9425145f, 0.9454438f, 0.9483672f, 0.9512847f, 
0.9541963f, 0.9571020f, 0.9600019f, 0.9628960f, 0.9657843f, 
0.9686668f, 0.9715436f, 0.9744146f, 0.9772799f, 0.9801396f, 
0.9829936f, 0.9858419f, 0.9886847f, 0.9915218f, 0.9943534f, 
0.9971795f
};

float P2_ARRAY[256] = {
1.0000000f, 1.0027113f, 1.0054299f, 1.0081559f, 1.0108893f,
1.0136301f, 1.0163783f, 1.0191340f, 1.0218971f, 1.0246678f,
1.0274459f, 1.0302316f, 1.0330249f, 1.0358257f, 1.0386341f,
1.0414501f, 1.0442738f, 1.0471051f, 1.0499441f, 1.0527908f,
1.0556452f, 1.0585073f, 1.0613772f, 1.0642549f, 1.0671404f,
1.0700337f, 1.0729349f, 1.0758439f, 1.0787608f, 1.0816856f,
1.0846184f, 1.0875591f, 1.0905077f, 1.0934644f, 1.0964291f,
1.0994018f, 1.1023826f, 1.1053714f, 1.1083684f, 1.1113735f,
1.1143867f, 1.1174082f, 1.1204378f, 1.1234756f, 1.1265216f,
1.1295759f, 1.1326385f, 1.1357094f, 1.1387886f, 1.1418762f,
1.1449721f, 1.1480765f, 1.1511892f, 1.1543104f, 1.1574401f,
1.1605782f, 1.1637249f, 1.1668800f, 1.1700438f, 1.1732161f,
1.1763970f, 1.1795865f, 1.1827847f, 1.1859916f, 1.1892071f,
1.1924314f, 1.1956644f, 1.1989062f, 1.2021567f, 1.2054161f,
1.2086843f, 1.2119614f, 1.2152474f, 1.2185422f, 1.2218460f,
1.2251588f, 1.2284805f, 1.2318113f, 1.2351511f, 1.2384999f,
1.2418578f, 1.2452248f, 1.2486010f, 1.2519863f, 1.2553808f,
1.2587844f, 1.2621974f, 1.2656195f, 1.2690510f, 1.2724917f,
1.2759418f, 1.2794012f, 1.2828700f, 1.2863482f, 1.2898359f,
1.2933330f, 1.2968396f, 1.3003556f, 1.3038813f, 1.3074164f,
1.3109612f, 1.3145156f, 1.3180796f, 1.3216533f, 1.3252366f,
1.3288297f, 1.3324325f, 1.3360451f, 1.3396675f, 1.3432997f,
1.3469418f, 1.3505937f, 1.3542555f, 1.3579273f, 1.3616090f,
1.3653007f, 1.3690024f, 1.3727142f, 1.3764360f, 1.3801679f,
1.3839099f, 1.3876620f, 1.3914244f, 1.3951969f, 1.3989797f,
1.4027727f, 1.4065760f, 1.4103896f, 1.4142136f, 1.4180479f,
1.4218926f, 1.4257477f, 1.4296133f, 1.4334894f, 1.4373760f,
1.4412731f, 1.4451808f, 1.4490991f, 1.4530280f, 1.4569676f,
1.4609178f, 1.4648787f, 1.4688504f, 1.4728329f, 1.4768261f,
1.4808302f, 1.4848452f, 1.4888710f, 1.4929077f, 1.4969554f,
1.5010141f, 1.5050837f, 1.5091644f, 1.5132562f, 1.5173590f,
1.5214730f, 1.5255982f, 1.5297345f, 1.5338820f, 1.5380408f,
1.5422108f, 1.5463922f, 1.5505849f, 1.5547889f, 1.5590044f,
1.5632313f, 1.5674696f, 1.5717195f, 1.5759808f, 1.5802538f,
1.5845383f, 1.5888344f, 1.5931422f, 1.5974616f, 1.6017928f,
1.6061357f, 1.6104903f, 1.6148568f, 1.6192351f, 1.6236253f,
1.6280274f, 1.6324415f, 1.6368675f, 1.6413054f, 1.6457555f,
1.6502176f, 1.6546918f, 1.6591781f, 1.6636766f, 1.6681873f,
1.6727102f, 1.6772454f, 1.6817928f, 1.6863526f, 1.6909248f,
1.6955094f, 1.7001064f, 1.7047158f, 1.7093378f, 1.7139722f,
1.7186193f, 1.7232789f, 1.7279512f, 1.7326362f, 1.7373338f,
1.7420442f, 1.7467674f, 1.7515034f, 1.7562522f, 1.7610138f,
1.7657884f, 1.7705760f, 1.7753765f, 1.7801900f, 1.7850166f,
1.7898563f, 1.7947091f, 1.7995750f, 1.8044542f, 1.8093465f,
1.8142522f, 1.8191711f, 1.8241034f, 1.8290490f, 1.8340081f,
1.8389806f, 1.8439666f, 1.8489661f, 1.8539791f, 1.8590058f,
1.8640460f, 1.8691000f, 1.8741676f, 1.8792490f, 1.8843442f,
1.8894532f, 1.8945760f, 1.8997127f, 1.9048633f, 1.9100280f,
1.9152066f, 1.9203992f, 1.9256059f, 1.9308268f, 1.9360618f,
1.9413110f, 1.9465744f, 1.9518521f, 1.9571441f, 1.9624505f,
1.9677712f, 1.9731064f, 1.9784560f, 1.9838202f, 1.9891988f,
1.9945921f
};


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
	union double_to_i64 ux;
	ux.d = x * RECIP_E_D;
	ux.i -= 0x3ff0000000000000ULL;
	unsigned int array_index = (unsigned int)(ux.i >> 44) & 0b11111111;
	ux.i >>= 52;
	double log = (double)(ux.i);
	log += (double)(LD_ARRAY[array_index]);
	ux.d = log * x;
	int exp = (int)(ux.d);
	double r = ux.d - (double)(exp);
	exp += 1023;
	ux.i = exp;
	ux.i <<= 52;
	r *= 256;
	array_index = (unsigned int)(r);
	return sqrt(TWO_PI_D * x) * ux.d * (double)(P2_ARRAY[array_index]);
}

/*
 * Approximates the gamma function for the single-precision
 * input variable x. The greater the argument is the more
 * precision the result is going to have. The function DOES
 * NOT perform any argument checks! This is done on purpose.
 */

//implement the faster version
FUNC_INLINE float fast_gammaf(float x) {
	x -= 1.0f;
	// standard Stirling formula
	union float_to_i32 ux;
	ux.f = x * RECIP_E_F;
	ux.i -= 0x3f800000UL;
	unsigned int array_index = (unsigned int)(ux.i >> 15) & 0b11111111;
	ux.i >>= 23;
	float log = (float)(ux.i);
	log += LD_ARRAY[array_index];
	ux.f = log * x;
	int exp = (int)(ux.f);
	float r = ux.f - (float)(exp);
	exp += 127;
	ux.i = exp;
	ux.i <<= 23;
	r *= 256;
	array_index = (unsigned int)(r);
	return sqrtf(TWO_PI_F * x) * ux.f * P2_ARRAY[array_index];
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
	ux.i &= 0x3fffffffffffffffULL;
#else
	ux.i -= 0x3ff0000000000000ULL;
	unsigned int array_index = (unsigned int)(ux.i >> 44) & 0b11111111;
#endif
	ux.i >>= 52;
	double log = (double)(ux.i);
#if USE_APPROX_ARRAYS != 0
	log += (double)(LD_ARRAY[array_index]);
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
	ux.i &= 0x3fffffffUL;
#else
	ux.i -= 0x3f800000UL;
	unsigned int array_index = (unsigned int)(ux.i >> 15) & 0b11111111;
#endif
	ux.i >>= 23;
	float log = (float)(ux.i);
#if USE_APPROX_ARRAYS != 0
	log += LD_ARRAY[array_index];
#endif
	log *= LOG_2_F; // a rough estimate of the natural log of x
#if USE_APPROX_ARRAYS == 0
	return x * log;
#else
	return x * (log - 1.0f);
#endif
}

#endif
