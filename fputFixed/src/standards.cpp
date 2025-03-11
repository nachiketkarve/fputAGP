#include "standards.hpp"

// Returns the element with the largest absolute value in a vector
double Mod(const Vector &X)
{
    return ((X.array()).abs()).maxCoeff();
}

// Computes the frequency of a given mode
double freq(int mode, int N)
{
    return 2.0 * sin(pi * double(mode) / (2.0 * double(N + 1)));
}

// Computes the sign of the input
int sgn(double val)
{
    return (0 < val) - (val < 0);
}