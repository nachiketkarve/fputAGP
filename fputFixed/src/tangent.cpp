#include "tangent.hpp"
#include "couplings.hpp"
#include "standards.hpp"

// Constructs the tangent vectors to a given state, projected on a constant energy surface and the mode amplitude space
void ConstructTangentReduced(Matrix &QTangent, State x, int seedMode)
{
    int N = x.N;
    if (QTangent.rows() != N + 2 || QTangent.cols() != N - 1)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    QTangent.setZero();
    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q(N + 2);
    Q = FourierComponents * x.q;

    int mode = 1;
    int colNum = 0;

    while (colNum < N - 1)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        QTangent(mode, colNum) = 1.0;
        Vector coeff(2);
        coeff.setZero();

        if (x.model == "alpha" || x.model == "toda")
        {
            coeff(0) = coeff(0) + freq(seedMode, N) * freq(seedMode, N) * Q(seedMode);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    coeff(0) = coeff(0) + x.nonLin * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
                }
            }
            for (int i = 0; i < N + 2; i++)
            {
                coeff(1) = coeff(1) + pow(freq(i, N), 2) * Q(i) * QTangent(i, colNum);
            }
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int k = 0; k < N + 2; k++)
                    {
                        coeff(1) = coeff(1) + x.nonLin * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * QTangent(k, colNum);
                    }
                }
            }
        }
        else if (x.model == "beta")
        {
            coeff(0) = coeff(0) + freq(seedMode, N) * freq(seedMode, N) * Q(seedMode);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int l = 0; l < N + 2; l++)
                    {
                        coeff(0) = coeff(0) + x.nonLin * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
                    }
                }
            }

            for (int i = 0; i < N + 2; i++)
            {
                coeff(1) = coeff(1) + pow(freq(i, N), 2) * Q(i) * QTangent(i, colNum);
            }
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int k = 0; k < N + 2; k++)
                    {
                        for (int l = 0; l < N + 2; l++)
                        {
                            coeff(1) = coeff(1) + x.nonLin * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * QTangent(k, colNum) * Q(l);
                        }
                    }
                }
            }
        }

        QTangent(seedMode, colNum) = -coeff(1) / coeff(0);

        mode = mode + 1;
        colNum = colNum + 1;
    }
}


// Constructs the tangent vectors to a given state, projected on a constant energy surface
void ConstructTangentFull(Matrix &QTangent, Matrix &PTangent, State x, int seedMode)
{
    int N = x.N;
    if (QTangent.rows() != N + 2 || QTangent.cols() != 2*N - 2 || PTangent.rows() != N + 2 || PTangent.cols() != 2*N - 2)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    QTangent.setZero();
    PTangent.setZero();

    int mode = 1;
    int colNum = 0;

    while (colNum < N - 1)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }

        QTangent(mode,colNum) = 1.0;
        mode = mode + 1;
        colNum = colNum + 1;
    }

    mode = 1;
    while (colNum < 2*N - 2)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }

        PTangent(mode,colNum) = 1.0;
        mode = mode + 1;
        colNum = colNum + 1;
    }
}

// Constructs the tangent vectors to a given state, with a fixed seed mode amplitude
void ConstructTangentFixedAmplitude(Matrix &QTangent, State x, int seedMode)
{
    int N = x.N;
    if (QTangent.rows() != N + 2 || QTangent.cols() != N - 1)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    QTangent.setZero();

    int mode = 1;
    int colNum = 0;

    while (colNum < N - 1)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        QTangent(mode, colNum) = 1.0;
        mode = mode + 1;
        colNum = colNum + 1;
    }
}

void ConstructTangentConstantEnergy(Matrix &QTangent, Matrix &PTangent, State &x, int seedMode)
{
    int N = x.N;
    if (QTangent.rows() != N + 2 || QTangent.cols() != 2*N - 2 || PTangent.rows() != N + 2 || PTangent.cols() != 2*N - 2)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    QTangent.setZero();
    PTangent.setZero();
    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q(N + 2);
    Vector P(N + 2);
    Q = FourierComponents * x.q;
    P = FourierComponents * x.p;

    int mode = 1;
    int colNum = 0;

    while (colNum < N - 1)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        QTangent(mode, colNum) = 1.0;
        Vector coeff(2);
        coeff.setZero();

        if (x.model == "alpha" || x.model == "toda")
        {
            coeff(0) = coeff(0) + freq(seedMode, N) * freq(seedMode, N) * Q(seedMode);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    coeff(0) = coeff(0) + x.nonLin * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
                }
            }
            for (int i = 0; i < N + 2; i++)
            {
                coeff(1) = coeff(1) + pow(freq(i, N), 2) * Q(i) * QTangent(i, colNum);
            }
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int k = 0; k < N + 2; k++)
                    {
                        coeff(1) = coeff(1) + x.nonLin * CouplingAlpha(i, j, k, N) * Q(i) * Q(j) * QTangent(k, colNum);
                    }
                }
            }
        }
        else if (x.model == "beta")
        {
            coeff(0) = coeff(0) + freq(seedMode, N) * freq(seedMode, N) * Q(seedMode);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int l = 0; l < N + 2; l++)
                    {
                        coeff(0) = coeff(0) + x.nonLin * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
                    }
                }
            }

            for (int i = 0; i < N + 2; i++)
            {
                coeff(1) = coeff(1) + pow(freq(i, N), 2) * Q(i) * QTangent(i, colNum);
            }
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int k = 0; k < N + 2; k++)
                    {
                        for (int l = 0; l < N + 2; l++)
                        {
                            coeff(1) = coeff(1) + x.nonLin * CouplingBeta(i, j, k, l, N) * Q(i) * Q(j) * QTangent(k, colNum) * Q(l);
                        }
                    }
                }
            }
        }

        QTangent(seedMode, colNum) = -coeff(1) / coeff(0);

        mode = mode + 1;
        colNum = colNum + 1;
    }

    mode = 1;

    while (colNum < 2*N - 2)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        PTangent(mode, colNum) = 1.0;
        Vector coeff(2);
        coeff.setZero();

        if (x.model == "alpha" || x.model == "toda")
        {
            coeff(0) = coeff(0) + freq(seedMode, N) * freq(seedMode, N) * Q(seedMode);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    coeff(0) = coeff(0) + x.nonLin * CouplingAlpha(i, j, seedMode, N) * Q(i) * Q(j);
                }
            }
            for (int i = 0; i < N + 2; i++)
            {
                coeff(1) = coeff(1) + P(i) * PTangent(i, colNum);
            }
        }
        else if (x.model == "beta")
        {
            coeff(0) = coeff(0) + freq(seedMode, N) * freq(seedMode, N) * Q(seedMode);
            for (int i = 0; i < N + 2; i++)
            {
                for (int j = 0; j < N + 2; j++)
                {
                    for (int l = 0; l < N + 2; l++)
                    {
                        coeff(0) = coeff(0) + x.nonLin * CouplingBeta(i, j, seedMode, l, N) * Q(i) * Q(j) * Q(l);
                    }
                }
            }

            for (int i = 0; i < N + 2; i++)
            {
                coeff(1) = coeff(1) + P(i) * PTangent(i, colNum);
            }
        }

        QTangent(seedMode, colNum) = -coeff(1) / coeff(0);

        mode = mode + 1;
        colNum = colNum + 1;
    }




}