/**
* Copyright Johannes Kloimböck 2021 - 2022.
* Distributed under the Boost Software License, Version 1.0.
* (See accompanying file LICENSE or copy at
* https://www.boost.org/LICENSE_1_0.txt)
*/

public class OGamma {
	// Three useful private constants that shouldn't have to be recalculated every time.
	private static final double TWOPI = Math.PI + Math.PI;
	private static final double LOGTWOPI = Math.log(TWOPI);
	private static final double RECE = 1.0 / Math.E;
	
	/*
	* If the input is larger than or equal to this number, a faster algorithm for the logGamma function is used
	*/
	private static final double STIRLING_THRESHOLD = Double.longBitsToDouble(0x41a0000000000000L); // 2^27
	
	/*
	* Approximates the gamma function for the input variable x.
	* The greater the argument is the more precision the result is going to have.
	* See https://en.wikipedia.org/wiki/Stirling%27s_approximation#Versions_suitable_for_calculators
	* for more details on this algorithm.
	*/
	private static double gammaApprox (double x) {
		return Math.sqrt(TWOPI / x) * Math.pow(RECE * (x + (1.0 / (12.0 * x - 0.1 / x))), x);
	}
	
	/*
	* Checks, whether x is an integer or not.
	*/
	private static boolean isInteger (double x) {
		return (double)((long)(x)) == x || x > 1E18;
	}
	
	/*
	* Calculates the gamma function of x and returns a double-precision floating-point number.
	*/
	private static double doubleGamma (double x) {
		if (x != x) { // sort out NaN values immediately
			return Double.NaN;
		}
		
		if (x > 0.0) {
			if (isInteger(x)) { // return precalculated values for gamma from an array, when the input is known to be an integer
				return x < 36.0 ? intFactArray[(int)(x) - 1] : Double.POSITIVE_INFINITY;
			} else {
				if (x >= 12.0) { // if x >= 12, the gammaApprox function will return a result that has at least 8 correct decimal places
					return gammaApprox(x);
				} else {
					/*
					* If x < 12, x has to be made larger for the gammaApprox function to return a "correct" result.
					* Then the obtained result is divided a few times. Rearranging the identity gamma(x) * x = gamma(x + 1)
					* we get gamma(x) / (x - 1) = gamma(x - 1). So gamma(9.3) = gamma(12.3) / (11.3 * 10.3 * 9.3).
					*/
					int diff = 12 - (int)(x);
					double y = x + diff;
					
					double r = gammaApprox(y);
					
					for (int i = 1; i <= diff; i++) {
						r /= (y - i);
					}
					
					return r;
				}
			}
		} else { // Here the identity gamma(x) * gamma(1 - x) = pi / sin(x * pi) is used.
			x = -x;
			
			if (isInteger(x)) {
				return Double.NaN;
			}
			
			double product = Math.PI / Math.sin(Math.PI * (x + 1));
			
			return product / doubleGamma(x + 1);
		}
	}
	
	/*
	* Returns gamma(x) as a single-precision floating-point number
	* x is also a single-precision floating-point number
	* If x is 0.0 or a negative integer, NaN is returned.
	* If gamma(x) would be too large for a single-precision float variable positive infinity is returned.
	*/
	public static float gamma (float x) {
		// calls doubleGamma to make sure that there are no rounding errors in the results
		return (float)(doubleGamma(x));
	}
	
	/*
	* Returns x! as a single-precision floating-point number
	* x is a 32-bit integer
	* If x is negative, the function returns NaN.
	* If the result would be too large for a single-precision float variable, positive infinity is returned.
	*/
	public static float factorial (int x) {
		return x >= 0 ? (float)doubleGamma((double)(x + 1)) : Float.NaN;
	}
	
	/*
	* Returns x! as a single-precision floating-point number
	* x is also a single-precision floating-point number
	* If x is negative, the function returns NaN.
	* If the result would be too large for a single-precision float variable, positive infinity is returned.
	*/
	public static float factorial (float x) {
		return x >= 0.0f ? gamma(x + 1.0f) : Float.NaN;
	}
	
	// This array holds the values for all factorials from 0! to 34! as single-precision floating-point numbers.
	private static final float [] intFactArray = { 1.0f, 1.0f, 2.0f, 6.0f, 24.0f, 120.0f, 720.0f, 5040.0f, 40320.0f,
						      362880.0f, 3628800.0f, 3.99168E7f, 4.790016E8f, 6.2270208E9f, 8.7178289E10f,
						      1.30767441E12f, 2.09227906E13f, 3.55687415E14f, 6.4023735E15f, 1.21645105E17f,
						      2.43290202E18f, 5.109094E19f, 1.1240007E21f, 2.5852017E22f, 6.204484E23f,
						      1.551121E25f, 4.0329146E26f, 1.0888869E28f, 3.0488835E29f, 8.841762E30f,
						      2.6525285E32f, 8.2228384E33f, 2.6313083E35f, 8.683318E36f, 2.952328E38f };
	
	/*
	* Approximates the natural logarithm of the gamma function for the input variable x.
	* The greater the argument is the more precision the result is going to have.
	* The advantage of having this function together with gammaApprox, is that logGammaApprox
	* can also handle much larger arguments and return a correct, meaningful result.
	* See https://en.wikipedia.org/wiki/Stirling%27s_approximation#Versions_suitable_for_calculators
	* for more details on this algorithm.
	*/
	private static double logGammaApprox (double x) {
		return 0.5 * (LOGTWOPI - Math.log(x)) + x * (Math.log(x + 1.0 / (12.0 * x - 0.1 / x)) - 1.0);
	}
	
	/*
	* Approximates the natural logarithm of the gamma function for the input variable x.
	* As this is a shortened version of the algorithm used in logGammaApprox, it can only be used
	* for numbers that are quite large. However, it should be noticably faster.
	*/
	private static double logGammaStirling (double x) {
		return x * (Math.log(x) - 1.0);
	}
	
	/*
	* Calculates the natural logarithm of the gamma function of x and returns a double-precision floating-point number.
	* Uses the same kind of algorithm as the doubleGamma function.
	*/
	private static double doubleLogGamma (double x) {
		if (x != x || x <= 0.0) { // only positive inputs are allowed
			return Double.NaN;
		}
		
		if (x == 1.0 || x == 2.0) {
			return 0.0;
		}
		
		if (x > 0.0) {
			if (x >= 12.0) {
				return x >= STIRLING_THRESHOLD ? logGammaStirling(x) : logGammaApprox(x);
			} else {
				int diff = 12 - (int)(x);
				double y = x + diff;
					
				double r = logGammaApprox(y);
				
				for (int i = 1; i <= diff; i++) {
					r -= Math.log(y - i);
				}
					
				return r;
			}
		} else {
			x = -x;
			
			if (isInteger(x)) {
				return Double.NaN;
			}
			
			double product = Math.PI / Math.sin(Math.PI * (x + 1));
			
			return Math.log(product) - doubleLogGamma(x);
		}
	}
	
	/*
	* Returns loggamma(x) as a single-precision floating-point number
	* The base of the logarithm is e (2.71828...)
	* x is also a single-precision floating-point number
	* If x is 0.0 or negative, NaN is returned.
	* If loggamma(x) would be too large for a single-precision float variable positive infinity is returned.
	*/
	public static float logGamma (float x) {
		return (float)(doubleLogGamma(x));
	}
	
	/*
	* Returns log(x!) as a single-precision floating-point number
	* x is a 32-bit integer
	* If x is negative, the function returns NaN.
	* If the result would be too large for a single-precision float variable, positive infinity is returned.
	*/
	public static float logFactorial (int x) {
		if (x == 0x7fffffff) {
			return (float)(doubleLogGamma((double)(0x80000000L)));
		}
		
		return x >= 0 ? (float)doubleLogGamma((double)(x + 1)) : Float.NaN;
	}
	
	/*
	* Returns log(x!) as a single-precision floating-point number
	* x is also a single-precision floating-point number
	* If x is negative, the function returns NaN.
	* If the result would be too large for a single-precision float variable, positive infinity is returned.
	*/
	public static float logFactorial (float x) {
		return x >= 0.0f ? logGamma(x + 1.0f) : Float.NaN;
	}
	
	// Rounds x down to an integer
	private static double downToIntegerValue (double x) {
		return x > 1E18 ? x : (long)(x);
	}
	
	/*
	* Calculates the subfactorial function !x as a double-precision floating-point number
	*/
	private static double doubleSubFactorial (int x) {
		return downToIntegerValue((doubleGamma((double)(x + 1)) + 1.0) * RECE);
	}
	
	/*
	* Returns !x as a single-precision floating-point number
	* x is a 32-bit integer
	* If x is negative, the function returns NaN.
	* If the result would be too large for a single-precision float variable, positive infinity is returned.
	*/
	public static float subFactorial (int x) {
		if (x < 0) {
			return Float.NaN;
		}
		
		return x < 35 ? (float)(doubleSubFactorial(x)) : Float.POSITIVE_INFINITY;
	}
	
	/*
	* Returns log(!x) as a single-precision floating-point number
	* x is a 32-bit integer
	* If x is negative, the function returns NaN.
	* If the result would be too large for a single-precision float variable, positive infinity is returned.
	*/
	public static float logSubFactorial (int x) {
		if (x < 0) {
			return Float.NaN;
		} else if (x == 0x7fffffff) {
			return (float)(doubleLogGamma((double)(0x80000000L)));
		}
		
		if (x < 35) {
			return (float)(Math.log(doubleSubFactorial(x)));
		} else {
			return (float)(doubleLogGamma((double)(x + 1)) - 1.0);
		}
	}
}
