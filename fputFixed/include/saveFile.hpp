#ifndef HEADERFILE_SAVEFILE
#define HEADERFILE_SAVEFILE

#include "state.hpp"
#include "eigenClasses.hpp"

// calculate the Floquet multipliers and save the periodic state in a csv file
void savePeriodicOrbitReduced(const State &periodic, int seedMode, double lmd, double deltaT, int iterations = 1, double tmax = 1000.0, std::string saveFolder = "./", std::string integrator = "rk4");
void savePeriodicOrbitFull(const State &periodic, int seedMode, double lmd, double deltaT, int iterations = 1, double tmax = 1000.0, std::string saveFolder = "./", std::string integrator = "rk4");
void savePeriodicOrbitFixedAmplitude(const State &periodic, int seedMode, double lmd, double deltaT, int iterations = 1, double tmax = 1000.0, std::string saveFolder = "./", std::string integrator = "rk4");
void savePeriodicOrbitReducedAGP(const State &periodic, int seedMode, double lmd, double deltaT, int iterations = 1, double tmax = 1000.0, std::string saveFolder = "./", std::string integrator = "rk4");
void savePeriodicOrbitReducedAGPGradient(const State &periodic, int seedMode, double E0init, double lmd, double deltaT, double deltaTAGPGradient, int iterations = 1, double tmax = 1000.0, std::string saveFolder = "./", std::string integrator = "rk4");

#endif