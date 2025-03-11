/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard headers and minor functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_STANDARDS
#define HEADERFILE_STANDARDS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <complex>
#include "eigenClasses.hpp"

// Returns the element with the largest absolute value in a vector
double Mod(const Vector &X);

// Computes the frequency of a given mode
double freq(int mode, int N);

// Computes the sign of the input
int sgn(double val);

#endif