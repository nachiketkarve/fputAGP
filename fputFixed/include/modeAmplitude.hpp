//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the amplitude of the seed mode at a given energy of the system
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_MODEAMPLITUDE
#define HEADERFILE_MODEAMPLITUDE

#include "eigenClasses.hpp"

// Computes the energy of the alpha-FPUT system described by the normal space coordinates Q and P
double EnergyNormalSpaceAlpha(const Vector &Q, const Vector &P, double alpha);

// Computes the energy of the beta-FPUT system described by the normal space coordinates Q and P
double EnergyNormalSpaceBeta(const Vector &Q, const Vector &P, double beta);

// Computes the amplitude of the seed mode such that the alpha-FPUT system has the specified energy
void modeAmplitudeAlpha(Vector &Q, Vector &P, int seedMode, double Energy, double alpha);

// Computes the amplitude of the seed mode such that the beta-FPUT system has the specified energy
void modeAmplitudeBeta(Vector &Q, Vector &P, int seedMode, double E0, double beta);

#endif