//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the state class which stores information about the state of the FPUT system
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_STATE
#define HEADERFILE_STATE

#include "eigenClasses.hpp"

// Objects of this class store the state of the system
class State
{
    public:
    
    // Size of the system
    int N;
    // Displacements on the lattice
    Vector q;
    // Momenta on the lattice
    Vector p;
    // Non-linearity parameter
    double nonLin;
    // Type of the FPUT model, possible options are "alpha", "beta", and "toda"
    std::string model;

    // Construct a state given the system size and model
    State(int n, std::string modelName);

    // Duplicate a state
    State(const State &x);
    
    // Return the total energy of the state
    double totalEnergy();

    // Return the interaction potential of the state
    double interactionPotential();

    // Return the phase space gradient of the interaction potential
    void potentialGradient(Vector &grad);

    // Integrate the equations of motion for one time step using the leap frog integrator
    void leapFrog(double deltaT);

    // Integrate the equations of motion of the state and tangent vectors for one time step using the leap frog integrator
    void leapFrogTangent(Matrix &qTangent, Matrix &pTangent, double deltaT);

    // Integrate the equations of motion for one time step using the RK4 integrator
    void rkn4(double deltaT);

    // Integrate the equations of motion of the state and tangent vectors for one time step using the RK4 integrator
    void rknTangent4(Matrix &qTangent, Matrix &pTangent, double deltaT);

    // Integrate the equations of motion for one time step using the RK6 integrator
    void rkn6(double deltaT);

    // Integrate the equations of motion of the state and tangent vectors for one time step using the RK6 integrator
    void rknTangent6(Matrix &qTangent, Matrix &pTangent, double deltaT);

    // Evolve the system for time tmax with a timestep of deltaT, possible integrator options: "leapFrog", "rk4", "rk6"
    void Evolve(double tmax, double deltaT, std::string integrator);

    // Evolve the system until it crosses the Poincare section described by its P_seedMode
    double EvolvePoincare(int seedMode, double deltaT, int iterations, double tmax, std::string integrator);

    // Evolve the system and its tangent vectors until it crosses the Poincare section described by its P_seedMode
    double EvolveTangentPoincare(Matrix &qTangent, Matrix &pTangent, int seedMode, double deltaT, int iterations, double tmax, std::string integrator);
    
    // Initializes the system into a single normal mode
    void initializeNormalMode(double Energy, int mode);


    private:

    // Computes the accelaration at each lattice point
    void eom(Vector &ddq, const Vector &Q);

    // Computes the accelaration at each lattice point and tangent vector
    void eomTangent(Vector &ddq, Matrix &ddqTangent, const Vector &Q, const Matrix &QTangent);

};

#endif