////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the AGP and its gradient on a periodic orbit
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_AGP
#define HEADERFILE_AGP

#include "state.hpp"
#include "eigenClasses.hpp"

// Calculate the Floquet matrix of a periodic state
double Floquet(Matrix &Floquet, const State &periodic, int seedMode, double deltaT, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4");

// Calculate the AGP of a periodic state
double AGPperiodic(const State &xPeriodic, int seedMode, double deltaT, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4");

// Calculate the AGP gradient of a periodic orbit starting at the time reversal symmetric point
void AGPGradientPeriodic(Vector &AGPGrad, const State &periodic, double period, int seedMode, Eigen::Vector<std::complex<double>, Eigen::Dynamic> &eigVals, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> &eigVecs, double deltaT, double tmax = 1000.0, int iterations = 1, std::string integrator = "rk4");

#endif