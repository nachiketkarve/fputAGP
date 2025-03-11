/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Construct the matrix of tangent vectors under various constraints
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_TANGENT
#define HEADERFILE_TANGENT

#include "eigenClasses.hpp"
#include "state.hpp"

// Constructs the tangent vectors to a given state, projected on a constant energy surface and the mode amplitude space
void ConstructTangentReduced(Matrix &QTangent, State x, int seedMode);

// Constructs the tangent vectors to a given state, projected on a constant energy surface
void ConstructTangentFull(Matrix &QTangent, Matrix &PTangent, State x, int seedMode);

// Constructs the tangent vectors to a given state, with a fixed seed mode amplitude
void ConstructTangentFixedAmplitude(Matrix &QTangent, State x, int seedMode);

void ConstructTangentConstantEnergy(Matrix &QTangent, Matrix &PTangent, State &x, int seedMode);

#endif