//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define iterations of the newton algorithm to compute periodic orbits under various constraints
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_NEWTON
#define HEADERFILE_NEWTON

#include "state.hpp"

// Perform one iteration of the Newton method on a given state x, to achieve a periodic orbit with a given energy and seed mode
double NewtonReduced(State &x, double Energy, int seedMode, double deltaT, bool log = true, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4");

// Perform one iteration of the Newton method on a given state x, to achieve a periodic orbit with a given seed mode amplitude and momentum
double NewtonFull(State &x, double QSeedMode, double PSeedMode, int seedMode, double deltaT, bool log = true, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4");

// Perform one iteration of the Newton method on a given state x, to achieve a periodic orbit with a given seed mode amplitude
double NewtonFixedAmplitude(State &x, double QSeedMode, int seedMode, double deltaT, bool log = true, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4");

double NewtonConstantEnergy(State &x, double Energy, double PSeedMode, int seedMode, double deltaT, bool log, int iterations, double tmax, std::string integrator);

#endif