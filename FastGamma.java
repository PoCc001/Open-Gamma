/**
* Copyright Johannes Kloimböck 2021.
* Distributed under the Boost Software License, Version 1.0.
* (See accompanying file LICENSE or copy at
* https://www.boost.org/LICENSE_1_0.txt)
*/

public class FastGamma {
	// Some useful private constants that shouldn't have to be recalculated every time.
	private static final double TWOPID = Math.PI + Math.PI;
	private static final float TWOPIF = (float)(TWOPID);
	private static final double LOG_2D = Math.log(2.0);
	private static final float LOG_2F = (float)(LOG_2D);
	private static final double RECED = 1.0 / Math.E;
	private static final float RECEF = (float)(RECED);
	
	/*
	* Approximates the gamma function for the double-precision input variable x.
	* The greater the argument is the more precision the result is going to have.
	* The function DOES NOT perform any argument checks! This is done on purpose.
	*/
	public static double gamma (double x) {
		x = x - 1.0;
		// standard Stirling formula
		return Math.sqrt(TWOPID * x) * Math.pow(x * RECED, x);
	}
	
	/*
	* Approximates the gamma function for the single-precision input variable x.
	* The greater the argument is the more precision the result is going to have.
	* The function DOES NOT perform any argument checks! This is done on purpose.
	*/
	public static float gamma (float x) {
		x = x - 1.0f;
		// standard Stirling formula
		return (float)(Math.sqrt(TWOPIF * x)) * (float)(Math.pow(x * RECEF, x));
	}
	
	/*
	* Approximates the natural logarithm of the gamma function for the double-precision input
	* variable x. The greater the argument is the more precision the result is going to have.
	* The function DOES NOT perform any argument checks! This is done on purpose.
	*/
	public static double logGamma (double x) {
		long bits = Double.doubleToRawLongBits(x);
		bits &= 0x3fffffffffffffffL;
		bits >>= 52;
		double log = bits;
		log *= LOG_2D; // a rough estimate of the natural log of x
		return x * log;
	}
	
	/*
	* Approximates the natural logarithm of the gamma function for the single-precision input
	* variable x. The greater the argument is the more precision the result is going to have.
	* The function DOES NOT perform any argument checks! This is done on purpose.
	*/
	public static float logGamma (float x) {
		int bits = Float.floatToRawIntBits(x);
		bits &= 0x3fffffff;
		bits >>= 23;
		float log = bits;
		log *= LOG_2F; // a rough estimate of the natural log of x
		return x * log;
	}
}
