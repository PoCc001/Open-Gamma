/**
* Copyright Johannes Kloimböck 2020.
* Distributed under the Boost Software License, Version 1.0.
* (See accompanying file LICENSE or copy at
* https://www.boost.org/LICENSE_1_0.txt)
*/

public strictfp class OGamma {
	private static double gammaApprox (double x) {
		return Math.sqrt(2.0 * Math.PI / x) * Math.pow(1.0 / Math.E * (x + (1.0 / (12.0 * x - 1.0 / (10.0 * x)))), x);
	}
	
	private static boolean isInteger (double x) {
		return (double)((long)(x)) == x || x > 1E18;
	}
	
	private static double doubleGamma (double x) {
		if (x != x) {
			return Double.NaN;
		} else if (x >= 35.0) {
			return Double.POSITIVE_INFINITY;
		}
		
		if (x > 0.0) {
			if (isInteger(x)) {
				return intFactArray[(int)(x) - 1];
			} else {
				if (x >= 12.0) {
					return gammaApprox(x);
				} else {
					int diff = 12 - (int)(x);
					double y = x + diff;
					
					double r = doubleGamma(y);
					
					for (int i = 0; i < diff; i++) {
						r /= (y - i - 1);
					}
					
					return r;
				}
			}
		} else {
			x = -x;
			
			if (isInteger(x)) {
				return Double.NaN;
			}
			
			double product = Math.PI / Math.sin(Math.PI * (x + 1));
			
			return product / doubleGamma(x + 1);
		}
	}
	
	public static float gamma (float x) {
		return (float)(doubleGamma(x));
	}
	
	public static float factorial (int x) {
		return x >= 0 ? doubleGamma((double)(x + 1)) : Float.NaN;
	}
	
	public static float factorial (float x) {
		return x >= 0.0f ? gamma(x + 1.0f) : Float.NaN;
	}
	
	private static final float [] intFactArray = { 1.0f, 1.0f, 2.0f, 6.0f, 24.0f, 120.0f, 720.0f, 5040.0f, 40320.0f,
						      362880.0f, 3628800.0f, 3.99168E7f, 4.790016E8f, 6.2270208E9f, 8.7178289E10f,
						      1.30767441E12f, 2.09227906E13f, 3.55687415E14f, 6.4023735E15f, 1.21645105E17f,
						      2.43290202E18f, 5.109094E19f, 1.1240007E21f, 2.5852017E22f, 6.204484E23f,
						      1.551121E25f, 4.0329146E26f, 1.0888869E28f, 3.0488835E29f, 8.841762E30f,
						      2.6525285E32f, 8.2228384E33f, 2.6313083E35f, 8.683318E36f, 2.952328E38f };
	
	private static double logGammaApprox (double x) {
		return 0.5 * (Math.log(2.0 * Math.PI) - Math.log(x)) + x * (Math.log(x + 1.0 / (12.0 * x - 1.0 / (10.0 * x))) - 1.0);
	}
	
	private static double doubleLogGamma (double x) {
		if (x != x || x <= 0.0) {
			return Double.NaN;
		}
		
		if (x == 1.0 || x == 2.0) {
			return 0.0;
		}
		
		if (x > 0.0) {
			if (x >= 12.0) {
				return logGammaApprox(x);
			} else {
				int diff = 12 - (int)(x);
				double y = x + diff;
					
				double r = doubleLogGamma(y);
				
				for (int i = 0; i < diff; i++) {
					r -= Math.log(y - i - 1);
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
	
	public static float logGamma (float x) {
		return (float)(doubleLogGamma(x));
	}
	
	public static float logFactorial (int x) {
		if (x == 0x7fffffff) {
			return (float)(doubleLogGamma((double)(0x80000000L)));
		}
		
		return x >= 0 ? doubleLogGamma((double)(x + 1)) : Float.NaN;
	}
	
	public static float logFactorial (float x) {
		return x >= 0.0f ? logGamma(x + 1.0f) : Float.NaN;
	}
	
	private static double downToIntegerValue (double x) {
		return x > 1E18 ? x : (long)(x);
	}
	
	private static double doubleSubFactorial (int x) {
		return downToIntegerValue((doubleGamma((double)(x + 1)) + 1.0) / Math.E)
	}
	
	public static float subFactorial (int x) {
		if (x < 0) {
			return Float.NaN;
		}
		
		return x < 35 ? (float)(doubleSubFactorial(x)) : Float.POSITIVE_INFINITY;
	}
	
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
