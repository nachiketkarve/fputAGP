#include "saveFile.hpp"
#include "standards.hpp"
#include "agp.hpp"

// calculate the Floquet multipliers and save the periodic state in a csv file
void savePeriodicOrbitReduced(const State &periodic, int seedMode, double lmd, double deltaT, int iterations, double tmax, std::string saveFolder, std::string integrator)
{
    State x = periodic;
    int N = periodic.N;
    Matrix Phi(2 * N, 2 * N);
    Matrix QTangent(N + 2, 2 * N);
    Matrix PTangent(N + 2, 2 * N);
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
    double E0 = x.totalEnergy();
    double nonLin = x.nonLin;

    Phi.setZero();
    for (int i = 0; i < 2 * N; i++)
    {
        Phi(i, i) = 1;
    }

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            QTangent(i, j) = Phi(i - 1, j);
            PTangent(i, j) = Phi(i + N - 1, j);
        }
    }

    Matrix qTangent(N + 2, 2 * N);
    Matrix pTangent(N + 2, 2 * N);
    qTangent = FourierComponents * QTangent;
    pTangent = FourierComponents * PTangent;

    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    QTangent = FourierComponents * qTangent;
    PTangent = FourierComponents * pTangent;

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            Phi(i - 1, j) = QTangent(i, j);
            Phi(i + N - 1, j) = PTangent(i, j);
        }
    }

    Eigen::EigenSolver<Matrix> eigensolver(Phi);
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> eigVals = eigensolver.eigenvalues();

    std::ofstream file;
    std::string FileName;
    if (periodic.model == "alpha" || periodic.model == "toda")
    {
        FileName = saveFolder + periodic.model + "Br-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (periodic.model == "beta")
    {
        FileName = saveFolder + periodic.model + "Br-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0) + "-B" + std::to_string(nonLin) + ".csv";
    }
    file.open(FileName);
    file << std::setprecision(15);
    file << "Q" << "," << "P" << "," << "EigValReal" << "," << "EigValImag" << "," << "N" << "," << x.model << "," << "k" << "," << "E0" << "," << "lmd" << "," << "period" << "\n";
    file << Q[0] << "," << P[0] << "," << std::real(eigVals(0)) << "," << std::imag(eigVals(0)) << "," << N << "," << nonLin << "," << seedMode << "," << E0 << "," << lmd << "," << period << "\n";
    for (int i = 1; i < eigVals.size(); i++)
    {
        if (i < Q.size())
        {
            file << Q(i) << "," << P(i);
        }
        else
        {
            file << "" << "," << "";
        }

        if (i < eigVals.size())
        {
            file << "," << std::real(eigVals(i)) << "," << std::imag(eigVals(i)) << "\n";
        }
        else
        {
            file << "\n";
        }
    }
}

void savePeriodicOrbitFull(const State &periodic, int seedMode, double lmd, double deltaT, int iterations, double tmax, std::string saveFolder, std::string integrator)
{
    State x = periodic;
    int N = periodic.N;
    Matrix Phi(2 * N, 2 * N);
    Matrix QTangent(N + 2, 2 * N);
    Matrix PTangent(N + 2, 2 * N);
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
    double E0 = x.totalEnergy();
    double nonLin = x.nonLin;
    double QSeedMode = Q(seedMode);
    double PSeedMode = P(seedMode);

    Phi.setZero();
    for (int i = 0; i < 2 * N; i++)
    {
        Phi(i, i) = 1;
    }

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            QTangent(i, j) = Phi(i - 1, j);
            PTangent(i, j) = Phi(i + N - 1, j);
        }
    }

    Matrix qTangent(N + 2, 2 * N);
    Matrix pTangent(N + 2, 2 * N);
    qTangent = FourierComponents * QTangent;
    pTangent = FourierComponents * PTangent;

    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    QTangent = FourierComponents * qTangent;
    PTangent = FourierComponents * pTangent;

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            Phi(i - 1, j) = QTangent(i, j);
            Phi(i + N - 1, j) = PTangent(i, j);
        }
    }

    Eigen::EigenSolver<Matrix> eigensolver(Phi);
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> eigVals = eigensolver.eigenvalues();

    std::ofstream file;
    std::string FileName = saveFolder + "Fixed-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-Q" + std::to_string(QSeedMode) + "-P" + std::to_string(PSeedMode) + "-A" + std::to_string(nonLin) + ".csv";
    file.open(FileName);
    file << std::setprecision(15);
    file << "Q" << "," << "P" << "," << "EigValReal" << "," << "EigValImag" << "," << "N" << "," << x.model << "," << "k" << "," << "Qk0" << "," << "Pk0" << "," << "lmd" << "," << "period" << "\n";
    file << Q[0] << "," << P[0] << "," << std::real(eigVals(0)) << "," << std::imag(eigVals(0)) << "," << N << "," << nonLin << "," << seedMode << "," << QSeedMode << "," << PSeedMode << "," << lmd << "," << period << "\n";
    for (int i = 1; i < eigVals.size(); i++)
    {
        if (i < Q.size())
        {
            file << Q(i) << "," << P(i);
        }
        else
        {
            file << "" << "," << "";
        }

        if (i < eigVals.size())
        {
            file << "," << std::real(eigVals(i)) << "," << std::imag(eigVals(i)) << "\n";
        }
        else
        {
            file << "\n";
        }
    }
}

void savePeriodicOrbitFixedAmplitude(const State &periodic, int seedMode, double lmd, double deltaT, int iterations, double tmax, std::string saveFolder, std::string integrator)
{
    State x = periodic;
    int N = periodic.N;
    Matrix Phi(2 * N, 2 * N);
    Matrix QTangent(N + 2, 2 * N);
    Matrix PTangent(N + 2, 2 * N);
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
    double Qk0 = Q(seedMode);
    double alpha = x.nonLin;

    Phi.setZero();
    for (int i = 0; i < 2 * N; i++)
    {
        Phi(i, i) = 1;
    }

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            QTangent(i, j) = Phi(i - 1, j);
            PTangent(i, j) = Phi(i + N - 1, j);
        }
    }

    Matrix qTangent(N + 2, 2 * N);
    Matrix pTangent(N + 2, 2 * N);
    qTangent = FourierComponents * QTangent;
    pTangent = FourierComponents * PTangent;

    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    QTangent = FourierComponents * qTangent;
    PTangent = FourierComponents * pTangent;

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            Phi(i - 1, j) = QTangent(i, j);
            Phi(i + N - 1, j) = PTangent(i, j);
        }
    }

    Eigen::EigenSolver<Matrix> eigensolver(Phi);
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> eigVals = eigensolver.eigenvalues();

    std::ofstream file;
    std::string FileName = saveFolder + "Fixed-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-Q" + std::to_string(Qk0) + "-A" + std::to_string(alpha) + ".csv";
    file.open(FileName);
    file << std::setprecision(15);
    file << "Q" << "," << "P" << "," << "EigValReal" << "," << "EigValImag" << "," << "N" << "," << x.model << "," << "k" << "," << "Qk0" << "," << "lmd" << "," << "period" << "\n";
    file << Q[0] << "," << P[0] << "," << std::real(eigVals(0)) << "," << std::imag(eigVals(0)) << "," << N << "," << alpha << "," << seedMode << "," << Qk0 << "," << lmd << "," << period << "\n";
    for (int i = 1; i < eigVals.size(); i++)
    {
        if (i < Q.size())
        {
            file << Q(i) << "," << P(i);
        }
        else
        {
            file << "" << "," << "";
        }

        if (i < eigVals.size())
        {
            file << "," << std::real(eigVals(i)) << "," << std::imag(eigVals(i)) << "\n";
        }
        else
        {
            file << "\n";
        }
    }
}


void savePeriodicOrbitReducedAGP(const State &periodic, int seedMode, double lmd, double deltaT, int iterations, double tmax, std::string saveFolder, std::string integrator)
{
    State x = periodic;
    int N = periodic.N;
    Matrix Phi(2 * N, 2 * N);
    Matrix QTangent(N + 2, 2 * N);
    Matrix PTangent(N + 2, 2 * N);
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
    double E0 = x.totalEnergy();
    if (x.model == "toda")
    {
        State xE = x;
        xE.model = "alpha";
        E0 = xE.totalEnergy();
    }
    
    double alpha = x.nonLin;

    Phi.setZero();
    for (int i = 0; i < 2 * N; i++)
    {
        Phi(i, i) = 1;
    }

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            QTangent(i, j) = Phi(i - 1, j);
            PTangent(i, j) = Phi(i + N - 1, j);
        }
    }

    Matrix qTangent(N + 2, 2 * N);
    Matrix pTangent(N + 2, 2 * N);
    qTangent = FourierComponents * QTangent;
    pTangent = FourierComponents * PTangent;

    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    QTangent = FourierComponents * qTangent;
    PTangent = FourierComponents * pTangent;

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            Phi(i - 1, j) = QTangent(i, j);
            Phi(i + N - 1, j) = PTangent(i, j);
        }
    }

    Eigen::EigenSolver<Matrix> eigensolver(Phi);
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> eigVals = eigensolver.eigenvalues();

    State xAGP = periodic;
    int numPts = int(period/deltaT);
    Vector AGPs(numPts);
    Vector times(numPts);
    AGPs.setZero();
    times.setZero();

    double currT = 0;
    double prevT = 0;
    double dH = xAGP.interactionPotential();
    double dHPrev = xAGP.interactionPotential();
    double AGP = 0;

    for (int i = 0; i < int(period/deltaT); i++)
    {
    
        dH = xAGP.interactionPotential();
        AGP = AGP - 0.5*(dH + dHPrev)*(currT-prevT);
        AGPs(i) = AGP;
        times(i) = currT;

        dHPrev = dH;
        prevT = currT;

        if (integrator == "rk6")
        {
            xAGP.rkn6(deltaT);
        } else
        {
            xAGP.rkn4(deltaT);
        }
        
        currT = currT + deltaT;
    }
    AGPs = AGPs - (times)*AGP/period;

    std::ofstream file;
    std::string FileName;
    if (x.model == "alpha")
    {
        FileName = saveFolder + "alphaBr-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0) + "-A" + std::to_string(alpha) + ".csv";
    } else if (x.model == "toda")
    {
        FileName = saveFolder + "todaBr-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0) + "-A" + std::to_string(alpha) + ".csv";
    } else if (x.model == "beta")
    {
        FileName = saveFolder + "betaBr-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0) + "-B" + std::to_string(alpha) + ".csv";
    }


    file.open(FileName);
    file << std::setprecision(15);

    double AGPmean = 0;
    double AGPnorm = 0;
    for (int i = 0; i < AGPs.size(); i++)
    {
        AGPnorm = AGPnorm + AGPs(i)*AGPs(i)/AGPs.size();
        AGPmean = AGPmean + AGPs(i)/AGPs.size();
    }

    AGPnorm = AGPnorm - AGPmean*AGPmean;

    file << "EigValReal";
    for (int i = 0; i < eigVals.size(); i++)
    {
        file << "," << std::real(eigVals(i));
    }
    file << "\n";

    file << "EigValImag";
    for (int i = 0; i < eigVals.size(); i++)
    {
        file << "," << std::imag(eigVals(i));
    }
    file << "\n";

    file << "Q";
    for (int i = 0; i < Q.size(); i++)
    {
        file << "," << Q(i);
    }
    file << "\n";

    file << "P";
    for (int i = 0; i < P.size(); i++)
    {
        file << "," << P(i);
    }
    file << "\n";
    
    file << "N" << "," << N << "\n";
    file << x.model << "," << x.nonLin << "\n";
    file << "k" << "," << seedMode << "\n";
    file << "E0" << "," << E0 << "\n";
    file << "lmd" << "," << lmd << "\n";
    file << "period" << "," << period << "\n";
    file << "AGPnorm";
    file << "," << AGPnorm;
    file << "\n";
    file.close();

}

void savePeriodicOrbitReducedAGPGradient(const State &periodic, int seedMode, double E0init, double lmd, double deltaT, double deltaTAGPGradient, int iterations, double tmax, std::string saveFolder, std::string integrator)
{
    State x = periodic;
    int N = periodic.N;
    Matrix Phi(2 * N, 2 * N);
    Matrix qTangent(N + 2, 2 * N);
    Matrix pTangent(N + 2, 2 * N);
    qTangent.setZero();
    pTangent.setZero();

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
    double E0 = x.totalEnergy();
    double alpha = x.nonLin;
    if (x.model == "toda")
    {
        State xE = x;
        xE.model = "alpha";
        E0 = xE.totalEnergy();
    }

    Phi.setZero();
    for (int i = 0; i < 2 * N; i++)
    {
        Phi(i, i) = 1;
    }

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            qTangent(i, j) = Phi(i - 1, j);
            pTangent(i, j) = Phi(i + N - 1, j);
        }
    }

    double period = x.EvolveTangentPoincare(qTangent, pTangent, seedMode, deltaT, iterations, tmax, integrator);

    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            Phi(i - 1, j) = qTangent(i, j);
            Phi(i + N - 1, j) = pTangent(i, j);
        }
    }

    Eigen::EigenSolver<Matrix> eigensolver(Phi);
    Eigen::Vector<std::complex<double>, Eigen::Dynamic> eigVals = eigensolver.eigenvalues();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigVecs = eigensolver.eigenvectors();
    

    Vector AGPgradient(2*N);
    AGPgradient.setZero();
    AGPGradientPeriodic(AGPgradient, periodic, period, seedMode, eigVals, eigVecs, deltaTAGPGradient, tmax, 1, integrator);

    // Vector QGrad(N+2);
    // QGrad.setZero();
    // for (int i = 1; i < N+1; i ++)
    // {
    //     QGrad(i) = AGPgradient(i-1);
    // }

    // QGrad = FourierComponents * QGrad;

    //std::cout << AGPgradient << "\n";
    //std::cout << QGrad << "\n";

    std::ofstream file;
    std::string FileName;

    if (x.model == "alpha")
    {
        FileName = saveFolder + "alphaBr-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0init) + "-A" + std::to_string(alpha) + ".csv";
    } else if (x.model == "toda")
    {
        FileName = saveFolder + "todaBr-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0init) + "-A" + std::to_string(alpha) + ".csv";
    } else if (x.model == "beta")
    {
        FileName = saveFolder + "betaBr-N" + std::to_string(N) + "-K" + std::to_string(int(seedMode)) + "-E" + std::to_string(E0init) + "-B" + std::to_string(alpha) + ".csv";
    }

    file.open(FileName);
    file << std::setprecision(15);
    file << "Q" << "," << "P" << "," << "EigValReal" << "," << "EigValImag" << "," << "AGPgrad" << "," << "N" << "," << x.model << "," << "k" << "," << "E0" << "," << "lmd" << "," << "period" << "\n";
    file << Q[0] << "," << P[0] << "," << std::real(eigVals(0)) << "," << std::imag(eigVals(0)) << "," << AGPgradient(0) << "," << N << "," << alpha << "," << seedMode << "," << E0 << "," << lmd << "," << period << "\n";
    for (int i = 1; i < eigVals.size(); i++)
    {
        if (i < Q.size())
        {
            file << Q(i) << "," << P(i);
        }
        else
        {
            file << "" << "," << "";
        }

        if (i < eigVals.size())
        {
            file << "," << std::real(eigVals(i)) << "," << std::imag(eigVals(i)) << "," << AGPgradient(i) <<"\n";
        }
        else
        {
            file << "\n";
        }
    }
}