#include "fputFixed.hpp"
#include "nlohmann/json.hpp"
#include "agp.hpp"

using json = nlohmann::json;

int main()
{
    std::ifstream file("_initBreather.json");
    json initData = json::parse(file);

    std::string model = initData["model"];
    int N = initData["N"];
    double E0 = initData["Energy"];
    int k0 = initData["seedMode"];
    double nonLinStart = initData["nonLinStart"];
    double nonLinEnd = initData["nonLinEnd"];
    double deltaNonLin = initData["deltaNonLin"];
    double deltaT = initData["deltaT"];
    double deltaTAGPGradient = initData["deltaTAGPGradient"];
    double tmax = initData["tmax"];
    int iterations = initData["iterations"];
    bool log = initData["log"];
    std::string integrator = initData["integrator"];
    std::string saveFolder = initData["saveFolder"];

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    double nonLin = nonLinStart;
    int count = 0;
    double lmd = 1.0;
    double E0init = E0;

    State x(N, model);
    x.nonLin = nonLin;
    x.initializeNormalMode(E0, k0);

    Matrix fq(2*N,2*N);

    while (nonLin <= nonLinEnd)
    {
        lmd = 1.0;
        count = 0;
        if (log && (model == "alpha" || model == "toda"))
        {
            std::cout << "Alpha = " << nonLin << "\n";
        } else if (log && (model == "beta"))
        {
            std::cout << "Beta = " << nonLin << "\n";
        }
        
        while (lmd > pow(10, -12))
        {
            lmd = NewtonReduced(x, E0, k0, deltaT, log, iterations, tmax, integrator);
            count = count + 1;
            if (lmd > pow(10, -1) || count > 15 || std::isnan(lmd) || std::isnan(x.totalEnergy()))
            {
                break;
            }
        }

        lmd = 1;

        while (lmd > pow(10,-12))
        {
            lmd = NewtonConstantEnergy(x,E0,0,k0,deltaT,log,iterations,tmax,integrator);
        }

        int nonLinE3 = round(nonLin * 1000000);

        if (lmd <= pow(10,-12))
        {
            double period = Floquet(fq,x,k0,deltaT,iterations,tmax,integrator);
            Eigen::EigenSolver<Matrix> eigensolver(fq);
            Eigen::Vector<std::complex<double>, Eigen::Dynamic> eigVals = eigensolver.eigenvalues();
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigVecs = eigensolver.eigenvectors();

            Vector AGPgradient(2*N);
            AGPgradient.setZero();
            AGPGradientPeriodic(AGPgradient, x, period, k0, eigVals, eigVecs, deltaTAGPGradient, tmax, 1, integrator);
            State xN = x;
            for (int i = 0; i < N; i++)
            {
                xN.q(i) = xN.q(i) + AGPgradient(N+i)*deltaNonLin;
                xN.p(i) = xN.p(i) - AGPgradient(i)*deltaNonLin;
            }

            State x2 = xN;
            x2.nonLin = nonLin + deltaNonLin;
            x2.EvolvePoincare(k0, deltaT, iterations, tmax, integrator);
            Vector Q2(N + 2);
            Vector Q1(N + 2);
            Q2 = FourierComponents * x2.q;
            Q1 = FourierComponents * xN.q;

            double lmd2 = Mod(Q2 - Q1) / Mod(Q1);
            std::cout << "Energy = " << xN.totalEnergy() << " lmd = " << lmd2 << "\n";

            E0 = xN.totalEnergy();
        }

        if (lmd <= pow(10, -12) && nonLinE3 % 10000 == 0)
        {
            savePeriodicOrbitReducedAGPGradient(x, k0, E0init, lmd, deltaT, deltaTAGPGradient, iterations, tmax, saveFolder, integrator);
        }

        nonLin = nonLin + deltaNonLin;
        x.nonLin = nonLin;
    }

    system("pause");
    // return 1;
}