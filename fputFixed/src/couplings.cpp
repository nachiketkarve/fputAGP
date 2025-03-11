#include "couplings.hpp"
#include "constants.hpp"
#include "eigenClasses.hpp"
#include "standards.hpp"

// Computes the coupling between the modes i, j, and k in an alpha-FPUT system with size N
double CouplingAlpha(int i, int j, int k, int N)
{
    double C = 0;
    if (i + j + k == 0)
    {
        C = 1.0;
    }
    if (i + j - k == 0)
    {
        C = 1.0;
    }
    if (i - j + k == 0)
    {
        C = 1.0;
    }
    if (i - j - k == 0)
    {
        C = 1.0;
    }

    if (i + j + k == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j - k == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j + k == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j - k == 2 * (N + 1))
    {
        C = -1.0;
    }
    C = C * freq(i, N) * freq(j, N) * freq(k, N) / sqrt(2.0 * double(N + 1));

    return C;
}

// Computes the coupling between the modes i, j, k, and l in an beta-FPUT system with size N
double CouplingBeta(int i, int j, int k, int l, int N)
{

    double C = 0;
    if (i + j + k + l == 0)
    {
        C = 1.0;
    }
    if (i - j + k + l == 0)
    {
        C = 1.0;
    }
    if (i + j - k + l == 0)
    {
        C = 1.0;
    }
    if (i + j + k - l == 0)
    {
        C = 1.0;
    }
    if (i - j - k + l == 0)
    {
        C = 1.0;
    }
    if (i - j + k - l == 0)
    {
        C = 1.0;
    }
    if (i + j - k - l == 0)
    {
        C = 1.0;
    }
    if (i - j - k - l == 0)
    {
        C = 1.0;
    }

    if (i + j + k + l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j + k + l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j - k + l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j + k - l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j - k + l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j + k - l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j - k - l == 2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j - k - l == 2 * (N + 1))
    {
        C = -1.0;
    }

    if (i + j + k + l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j + k + l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j - k + l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j + k - l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j - k + l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j + k - l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i + j - k - l == -2 * (N + 1))
    {
        C = -1.0;
    }
    if (i - j - k - l == -2 * (N + 1))
    {
        C = -1.0;
    }

    C = C * freq(i, N) * freq(j, N) * freq(k, N) * freq(l, N) / (2.0 * double(N + 1));

    return C;
}