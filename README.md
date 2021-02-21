# Open-Gamma
## OGamma
Open-Source implementation of the Gamma function (and Factorial) for Java (float datatype).

Use the gamma function in the OGamma class to calculate gamma(x) for x is any real number except for non-positive integers. The factorial function only allows positive real numbers. You can enter any number between about -35 and +35.

In order to be able to evaluate the gamma/factorial function of much larger numbers use the logGamma and logFactorial functions. Note, that they currently only take positive inputs.

## FastGamma
Open-Source implementation of the Gamma and LogGamma function for Java (double and float datatype).

These implementation don't do any argument checks and are only useful for relatively large positive input values.
Therefore, they're quite fast, but at the same time much more inprecise.
