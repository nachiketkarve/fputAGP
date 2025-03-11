#include "newton.hpp"
#include "standards.hpp"
#include "tangent.hpp"
#include "modeAmplitude.hpp"

// Perform one iteration of the Newton method on a given state x, to achieve a periodic orbit with a given energy and seed mode
double NewtonReduced(State &x, double Energy, int seedMode, double deltaT, bool log, int iterations, double tmax, std::string integrator)
{
    int N = x.N;

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q0(N + 2);
    Q0 = FourierComponents * x.q;

    Vector Z0(N - 1);
    int mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z0(i) = Q0(mode);
        mode = mode + 1;
    }

    Matrix QTangent(N + 2, N - 1);
    Matrix qTangent(N + 2, N - 1);
    Matrix pTangent(N + 2, N - 1);
    ConstructTangentReduced(QTangent, x, seedMode);
    qTangent = FourierComponents * QTangent;
    pTangent.setZero();
    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    Vector Q(N + 2);
    Vector P(N + 2);
    Q = FourierComponents * x.q;

    Vector Z(N - 1);
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z(i) = Q(mode);
        mode = mode + 1;
    }

    QTangent = FourierComponents * qTangent;

    Matrix Gradient(N - 1, N - 1);
    Gradient.setZero();
    for (int colNum = 0; colNum < N - 1; colNum++)
    {
        int rowNum = 0;
        mode = 1;
        while (rowNum < N - 1)
        {
            if (mode == seedMode)
            {
                mode = mode + 1;
            }
            Gradient(rowNum, colNum) = QTangent(mode, colNum);
            rowNum = rowNum + 1;
            mode = mode + 1;
        }
    }

    Matrix Id(N - 1, N - 1);
    Id.setZero();
    for (int i = 0; i < N - 1; i++)
    {
        Id(i, i) = 1.0;
    }

    Matrix Mat(N - 1, N - 1);
    Mat.setZero();
    Mat = Gradient - Id;

    Z = Z0 - (Mat.inverse()) * (Z - Z0);
    Q.setZero();
    P.setZero();
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Q(mode) = Z(i);
        mode = mode + 1;
    }

    if (x.model == "alpha" || x.model == "toda")
    {
        modeAmplitudeAlpha(Q, P, seedMode, Energy, x.nonLin);
    }
    else if (x.model == "beta")
    {
        modeAmplitudeBeta(Q, P, seedMode, Energy, x.nonLin);
    }

    x.q = FourierComponents * Q;
    x.p.setZero();

    State x1 = x;
    period = x1.EvolvePoincare(seedMode, deltaT, iterations, tmax, integrator);
    Vector Q1(N + 2);
    Vector P1(N + 2);
    Q1 = FourierComponents * x1.q;
    P1.setZero();
    if (x.model == "alpha" || x.model == "toda")
    {
        modeAmplitudeAlpha(Q1, P1, seedMode, Energy, x1.nonLin);
    }
    else if (x.model == "beta")
    {
        modeAmplitudeBeta(Q1, P1, seedMode, Energy, x1.nonLin);
    }

    double lmd = Mod(Q1 - Q) / Mod(Q);

    if (log)
    {
        std::cout << "Period = " << period << "\n"
                    << "Energy = " << x.totalEnergy() << " lmd = " << lmd << "\n";
    }

    return lmd;
}

// Perform one iteration of the Newton method on a given state x, to achieve a periodic orbit with a given seed mode amplitude
// and momentum
double NewtonFull(State &x, double QSeedMode, double PSeedMode, int seedMode, double deltaT, bool log, int iterations, double tmax, std::string integrator)
{
    int N = x.N;

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q0(N + 2);
    Vector P0(N + 2);
    Q0 = FourierComponents * x.q;
    P0 = FourierComponents * x.p;

    Vector Z0(2*N - 2);
    int mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z0(i) = Q0(mode);
        mode = mode + 1;
    }

    mode = 1;
    for (int i = N - 1; i < 2*N - 2; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z0(i) = P0(mode);
        mode = mode + 1;
    }

    Matrix QTangent(N + 2, 2*N - 2);
    Matrix PTangent(N + 2, 2*N - 2);
    Matrix qTangent(N + 2, 2*N - 2);
    Matrix pTangent(N + 2, 2*N - 2);
    ConstructTangentFull(QTangent, PTangent, x, seedMode);
    qTangent = FourierComponents * QTangent;
    pTangent = FourierComponents * PTangent;
    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    Vector Q(N + 2);
    Vector P(N + 2);
    Q = FourierComponents * x.q;
    P = FourierComponents * x.p;

    Vector Z(2*N - 2);
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z(i) = Q(mode);
        mode = mode + 1;
    }
    mode = 1;
    for (int i = N - 1; i < 2*N - 2; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z(i) = P(mode);
        mode = mode + 1;
    }

    QTangent = FourierComponents * qTangent;
    PTangent = FourierComponents * pTangent;

    Matrix Gradient(2*N - 2, 2*N - 2);
    Gradient.setZero();
    for (int colNum = 0; colNum < 2*N - 2; colNum++)
    {
        int rowNum = 0;
        mode = 1;
        while (rowNum < N - 1)
        {
            if (mode == seedMode)
            {
                mode = mode + 1;
            }
            Gradient(rowNum, colNum) = QTangent(mode, colNum);
            rowNum = rowNum + 1;
            mode = mode + 1;
        }
        mode = 1;
        while (rowNum < 2*N - 2)
        {
            if (mode == seedMode)
            {
                mode = mode + 1;
            }
            Gradient(rowNum, colNum) = PTangent(mode, colNum);
            rowNum = rowNum + 1;
            mode = mode + 1;
        }
    }

    Matrix Id(2*N - 2, 2*N - 2);
    Id.setZero();
    for (int i = 0; i < 2*N - 2; i++)
    {
        Id(i, i) = 1.0;
    }

    Matrix Mat(2*N - 2, 2*N - 2);
    Mat.setZero();
    Mat = Gradient - Id;

    Z = Z0 - (Mat.inverse()) * (Z - Z0);
    Q.setZero();
    P.setZero();
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Q(mode) = Z(i);
        mode = mode + 1;
    }
    mode = 1;
    for (int i = N - 1; i < 2*N - 2; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        P(mode) = Z(i);
        mode = mode + 1;
    }

    Q(seedMode) = QSeedMode;
    P(seedMode) = PSeedMode;

    x.q = FourierComponents * Q;
    x.p = FourierComponents * P;

    State x1 = x;
    period = x1.EvolvePoincare(seedMode, deltaT, iterations, tmax, integrator);
    Vector Q1(N + 2);
    Vector P1(N + 2);
    Q1 = FourierComponents * x1.q;
    P1 = FourierComponents * x1.p;
    Q1(seedMode) = QSeedMode;
    P1(seedMode) = PSeedMode;

    double lmd = Mod(Q1 - Q) / Mod(Q);

    if (log)
    {
        std::cout << "Period = " << period << "\n"
                    << "Energy = " << x.totalEnergy() << " lmd = " << lmd << "\n";
    }

    return lmd;
}

// Perform one iteration of the Newton method on a given state x, to achieve a periodic orbit with a given seed mode amplitude
double NewtonFixedAmplitude(State &x, double QSeedMode, int seedMode, double deltaT, bool log, int iterations, double tmax, std::string integrator)
{
    int N = x.N;

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q0(N + 2);
    Q0 = FourierComponents * x.q;

    Vector Z0(N - 1);
    int mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z0(i) = Q0(mode);
        mode = mode + 1;
    }

    Matrix QTangent(N + 2, N - 1);
    Matrix qTangent(N + 2, N - 1);
    Matrix pTangent(N + 2, N - 1);
    ConstructTangentFixedAmplitude(QTangent, x, seedMode);
    qTangent = FourierComponents * QTangent;
    pTangent.setZero();
    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    Vector Q(N + 2);
    Q = FourierComponents * x.q;

    Vector Z(N - 1);
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z(i) = Q(mode);
        mode = mode + 1;
    }

    QTangent = FourierComponents * qTangent;

    Matrix Gradient(N - 1, N - 1);
    Gradient.setZero();
    for (int colNum = 0; colNum < N - 1; colNum++)
    {
        int rowNum = 0;
        mode = 1;
        while (rowNum < N - 1)
        {
            if (mode == seedMode)
            {
                mode = mode + 1;
            }
            Gradient(rowNum, colNum) = QTangent(mode, colNum);
            rowNum = rowNum + 1;
            mode = mode + 1;
        }
    }

    Matrix Id(N - 1, N - 1);
    Id.setZero();
    for (int i = 0; i < N - 1; i++)
    {
        Id(i, i) = 1.0;
    }

    Matrix Mat(N - 1, N - 1);
    Mat.setZero();
    Mat = Gradient - Id;

    Z = Z0 - (Mat.inverse()) * (Z - Z0);
    Q.setZero();
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Q(mode) = Z(i);
        mode = mode + 1;
    }

    Q(seedMode) = QSeedMode;

    x.q = FourierComponents * Q;
    x.p.setZero();

    State x1 = x;
    period = x1.EvolvePoincare(seedMode, deltaT, iterations, tmax, integrator);
    Vector Q1(N + 2);
    Q1 = FourierComponents * x1.q;
    Q1(seedMode) = QSeedMode;

    double lmd = Mod(Q1 - Q) / Mod(Q);

    if (log)
    {
        std::cout << "Period = " << period << "\n"
                    << "Energy = " << x.totalEnergy() << " lmd = " << lmd << "\n";
    }

    return lmd;
}


double NewtonConstantEnergy(State &x, double Energy, double PSeedMode, int seedMode, double deltaT, bool log, int iterations, double tmax, std::string integrator)
{
    int N = x.N;

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q0(N + 2);
    Vector P0(N + 2);
    Q0 = FourierComponents * x.q;
    P0 = FourierComponents * x.p;

    Vector Z0(2*N - 2);
    int mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z0(i) = Q0(mode);
        mode = mode + 1;
    }

    mode = 1;
    for (int i = N - 1; i < 2*N - 2; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z0(i) = P0(mode);
        mode = mode + 1;
    }

    Matrix QTangent(N + 2, 2*N - 2);
    Matrix PTangent(N + 2, 2*N - 2);
    Matrix qTangent(N + 2, 2*N - 2);
    Matrix pTangent(N + 2, 2*N - 2);
    ConstructTangentConstantEnergy(QTangent, PTangent, x, seedMode);
    qTangent = FourierComponents * QTangent;
    pTangent = FourierComponents * PTangent;
    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    Vector Q(N + 2);
    Vector P(N + 2);
    Q = FourierComponents * x.q;
    P = FourierComponents * x.p;

    Vector Z(2*N - 2);
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z(i) = Q(mode);
        mode = mode + 1;
    }
    mode = 1;
    for (int i = N - 1; i < 2*N - 2; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Z(i) = P(mode);
        mode = mode + 1;
    }

    QTangent = FourierComponents * qTangent;
    PTangent = FourierComponents * pTangent;

    Matrix Gradient(2*N - 2, 2*N - 2);
    Gradient.setZero();
    for (int colNum = 0; colNum < 2*N - 2; colNum++)
    {
        int rowNum = 0;
        mode = 1;
        while (rowNum < N - 1)
        {
            if (mode == seedMode)
            {
                mode = mode + 1;
            }
            Gradient(rowNum, colNum) = QTangent(mode, colNum);
            rowNum = rowNum + 1;
            mode = mode + 1;
        }
        mode = 1;
        while (rowNum < 2*N - 2)
        {
            if (mode == seedMode)
            {
                mode = mode + 1;
            }
            Gradient(rowNum, colNum) = PTangent(mode, colNum);
            rowNum = rowNum + 1;
            mode = mode + 1;
        }
    }

    Matrix Id(2*N - 2, 2*N - 2);
    Id.setZero();
    for (int i = 0; i < 2*N - 2; i++)
    {
        Id(i, i) = 1.0;
    }

    Matrix Mat(2*N - 2, 2*N - 2);
    Mat.setZero();
    Mat = Gradient - Id;

    Z = Z0 - (Mat.inverse()) * (Z - Z0);
    Q.setZero();
    P.setZero();
    mode = 1;
    for (int i = 0; i < N - 1; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        Q(mode) = Z(i);
        mode = mode + 1;
    }
    mode = 1;
    for (int i = N - 1; i < 2*N - 2; i++)
    {
        if (mode == seedMode)
        {
            mode = mode + 1;
        }
        P(mode) = Z(i);
        mode = mode + 1;
    }

    P(seedMode) = PSeedMode;
    if(x.model == "alpha" || x.model == "toda")
    {
        modeAmplitudeAlpha(Q, P, seedMode, Energy, x.nonLin);
    } else if (x.model == "beta")
    {
        modeAmplitudeBeta(Q, P, seedMode, Energy, x.nonLin);
    }

    x.q = FourierComponents * Q;
    x.p = FourierComponents * P;

    State x1 = x;
    period = x1.EvolvePoincare(seedMode, deltaT, iterations, tmax, integrator);
    Vector Q1(N + 2);
    Vector P1(N + 2);
    Q1 = FourierComponents * x1.q;
    P1 = FourierComponents * x1.p;
    P1(seedMode) = PSeedMode;
    if(x.model == "alpha" || x.model == "toda")
    {
        modeAmplitudeAlpha(Q1, P1, seedMode, Energy, x.nonLin);
    } else if (x.model == "beta")
    {
        modeAmplitudeBeta(Q1, P1, seedMode, Energy, x.nonLin);
    }

    double lmd = Mod(Q1 - Q) / Mod(Q);

    if (log)
    {
        std::cout << "Period = " << period << "\n"
                    << "Energy = " << x.totalEnergy() << " lmd = " << lmd << "\n";
    }

    return lmd;
}