//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the cubic and quartic couplings in the mode space of the FPUT system
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_COUPLINGS
#define HEADERFILE_COUPLINGS

// Computes the coupling between the modes i, j, and k in an alpha-FPUT system with size N
double CouplingAlpha(int i, int j, int k, int N);

// Computes the coupling between the modes i, j, k, and l in an beta-FPUT system with size N
double CouplingBeta(int i, int j, int k, int l, int N);

#endif