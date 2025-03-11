//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the matrix and vector classes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEADERFILE_EIGENCLASSES
#define HEADERFILE_EIGENCLASSES
#include "Eigen/Dense"
#include "unsupported/Eigen/CXX11/Tensor"
#include "unsupported/Eigen/Polynomials"
#include "Eigen/LU"
#include "Eigen/Eigenvalues"
#include "constants.hpp"

// Define Matrix and Vector classes
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Vector<double, Eigen::Dynamic> Vector;

#endif