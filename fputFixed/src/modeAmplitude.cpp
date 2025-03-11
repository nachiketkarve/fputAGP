#include "modeAmplitude.hpp"
#include "standards.hpp"
#include "couplings.hpp"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Computes the energy of the alpha-FPUT system described by the normal space coordinates Q and P
double EnergyNormalSpaceAlpha(const Vector &Q, const Vector &P, double alpha)
{
    int N = Q.size() - 2;
    double Etot = 0;

    if (P.size() != N + 2)
    {
        std::cout << "Vector length does not match";
        system("pause");
        exit(0);
    }

    for (int i = 0; i < N + 2; i++)
    {
        Etot = Etot + 0.5 * pow(P(i), 2) + 0.5 * pow(freq(i, N) * Q(i), 2);
    }

    for (int j = 0; j < N + 2; j++)
    {
        for (int k = 0; k < N + 2; k++)
        {
            int i = j + k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = j - k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = -j + k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = -j - k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = 2 * (N + 1) + j + k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = 2 * (N + 1) + j - k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = 2 * (N + 1) - j + k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
            i = 2 * (N + 1) - j - k;
            if (0 <= i && i < N + 2)
            {
                Etot = Etot + alpha / 3.0 * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * Q(k);
            }
        }
    }
    return Etot;
}

// Computes the energy of the beta-FPUT system described by the normal space coordinates Q and P
double EnergyNormalSpaceBeta(const Vector &Q, const Vector &P, double beta)
{
    int N = Q.size() - 2;
    double Etot = 0;

    if (P.size() != N + 2)
    {
        std::cout << "Vector length does not match";
        system("pause");
        exit(0);
    }

    for (int i = 0; i < N + 2; i++)
    {
        Etot = Etot + 0.5 * pow(P(i), 2) + 0.5 * pow(freq(i, N) * Q(i), 2);
    }

    for (int j = 0; j < N + 2; j++)
    {
        for (int k = 0; k < N + 2; k++)
        {
            for (int l = 0; l < N + 2; l++)
            {
                int i = j + k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -j + k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = j - k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = j + k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -j - k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -j + k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = j - k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) + j + k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) - j + k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) + j - k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) + j + k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) - j - k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) - j + k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) + j - k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = 2 * (N + 1) - j - k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) + j + k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) - j + k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) + j - k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) + j + k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) - j - k + l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) - j + k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) + j - k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }

                i = -2 * (N + 1) - j - k - l;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + beta / 4.0 * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * Q(k) * Q(l);
                }
            }
        }
    }
    return Etot;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Computes the amplitude of the seed mode such that the alpha-FPUT system has the specified energy
void modeAmplitudeAlpha(Vector &Q, Vector &P, int seedMode, double Energy, double alpha)
{
    int N = Q.size() - 2;
    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Q(seedMode) = 0.0;

    Vector coeff(3);
    coeff.setZero();
    coeff(0) = EnergyNormalSpaceAlpha(Q, P, alpha) - Energy;
    for (int j = 0; j < N + 2; j++)
    {
        int i = j + seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = j - seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = -j + seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = -j - seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = 2 * (N + 1) + j + seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = 2 * (N + 1) + j - seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = 2 * (N + 1) - j + seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
        i = 2 * (N + 1) - j - seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(1) = coeff(1) + alpha * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
        }
    }
    coeff(2) = coeff(2) + 0.5 * pow(freq(seedMode, N), 2);

    int i = seedMode + seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(2) = coeff(2) + alpha * CouplingAlpha(i, seedMode, seedMode, N) * Q(i);
    }
    i = 2 * (N + 1) - seedMode - seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(2) = coeff(2) + alpha * CouplingAlpha(i, seedMode, seedMode, N) * Q(i);
    }

    Eigen::PolynomialSolver<double, 2> solver;
    solver.compute(coeff);
    bool realEx;
    Q(seedMode) = solver.greatestRealRoot(realEx, 0.000001);
}

// Computes the amplitude of the seed mode such that the beta-FPUT system has the specified energy
void modeAmplitudeBeta(Vector &Q, Vector &P, int seedMode, double E0, double beta)
{
    int N = Q.size() - 2;
    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector coeff(5);
    coeff.setZero();
    Q(seedMode) = 0.0;
    coeff(0) = EnergyNormalSpaceBeta(Q, P, beta) - E0;
    for (int j = 0; j < N + 2; j++)
    {
        for (int l = 0; l < N + 2; l++)
        {
            int i = j + seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -j + seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = j - seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = j + seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -j - seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -j + seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = j - seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -j - seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) + j + seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) - j + seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) + j - seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) + j + seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) - j - seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) - j + seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) + j - seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = 2 * (N + 1) - j - seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) + j + seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) - j + seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) + j - seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) + j + seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) - j - seedMode + l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) - j + seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) + j - seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }

            i = -2 * (N + 1) - j - seedMode - l;
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + beta * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
            }
        }
    }
    coeff(2) = coeff(2) + 0.5 * pow(freq(seedMode, N), 2);

    for (int j = 0; j < N + 2; j++)
    {
        int i = j + 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -j + 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = j;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 12.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -j;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 12.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = j - 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -j - 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = 2 * (N + 1) + j + 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = 2 * (N + 1) - j + 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = 2 * (N + 1) + j;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 12.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = 2 * (N + 1) - j;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 12.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = 2 * (N + 1) + j - 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = 2 * (N + 1) - j - 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -2 * (N + 1) + j + 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -2 * (N + 1) - j + 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -2 * (N + 1) + j;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 12.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -2 * (N + 1) - j;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 12.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -2 * (N + 1) + j - 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }

        i = -2 * (N + 1) - j - 2 * seedMode;
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + 6.0 / 4.0 * beta * CouplingBeta(i, j, seedMode, seedMode, N) * Q(i) * Q(j);
        }
    }

    int i = 3 * seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(3) = coeff(3) + beta * CouplingBeta(i, seedMode, seedMode, seedMode, N) * Q(i);
    }

    i = seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(3) = coeff(3) + 3.0 * beta * CouplingBeta(i, seedMode, seedMode, seedMode, N) * Q(i);
    }

    i = 2 * (N + 1) - seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(3) = coeff(3) + 3.0 * beta * CouplingBeta(i, seedMode, seedMode, seedMode, N) * Q(i);
    }

    i = 2 * (N + 1) - 3 * seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(3) = coeff(3) + beta * CouplingBeta(i, seedMode, seedMode, seedMode, N) * Q(i);
    }

    i = -2 * (N + 1) + 3 * seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(3) = coeff(3) + beta * CouplingBeta(i, seedMode, seedMode, seedMode, N) * Q(i);
    }

    i = -2 * (N + 1) + seedMode;
    if (0 <= i && i < N + 2)
    {
        coeff(3) = coeff(3) + 3.0 * beta * CouplingBeta(i, seedMode, seedMode, seedMode, N) * Q(i);
    }

    coeff(4) = coeff(4) + 3.0 / 4.0 * beta * CouplingBeta(seedMode, seedMode, seedMode, seedMode, N);

    Eigen::PolynomialSolver<double, 4> solver;
    solver.compute(coeff);
    bool realEx;
    Q(seedMode) = solver.greatestRealRoot(realEx, 0.000001);
}