# Open-Gamma
## OGamma
Open-Source implementation of the Gamma function (and Factorial) in Java (float datatype).

Use the gamma function in the OGamma class to calculate gamma(x) for x is any real number except for non-positive integers. The factorial function only allows positive real numbers. You can enter any number between about -35 and +35.

In order to be able to evaluate the gamma/factorial function of much larger numbers use the logGamma and logFactorial functions. Note, that they currently only take positive inputs.

## FastGamma
Open-Source implementation of the Gamma and LogGamma function in Java (double and float datatype). Note, that they only return a rough approximation of the actual value! (2-3 digits maximum) For small values, the result may even be completely wrong!

These implementation don't do any argument checks and are only useful for relatively large positive input values.
Therefore, they're quite fast, but at the same time much more inprecise.

## fastgamma
Open-Source implementation of the Gamma and LogGamma function in C/C++ (double and float datatype). Note, that they only return a rough approximation of the actual value (2-3 digits maximum)! For small values, the result may even be completely wrong! You can however get more correct digits by assigning the value 1 to the macros
SUBSTITUTE_LOG_FUNCTION and SUBSTITUTE_POW_FUNCTION (5-6 digits max).

To use this part of the library, just include the fastgamma.h file in your C/C++ project.

These implementation don't do any argument checks and are only useful for relatively large positive input values.
Therefore, they're quite fast, but at the same time much more inprecise.

If, for any reason, you don't want the funcions to be inlined, change the DO_INLINE macro to 0.
