#include "fputFixed.hpp"
#include "nlohmann/json.hpp"

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

    State xb = x;

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
            //lmd = NewtonConstantEnergy(x,E0,0,k0,deltaT,log,iterations,tmax,integrator);
            lmd = NewtonReduced(x, E0, k0, deltaT, log, iterations, tmax, integrator);
            count = count + 1;
            if (lmd > pow(10, -1) || count > 15 || std::isnan(lmd) || std::isnan(x.totalEnergy()))
            {
                x = xb;
                break;
            }
        }
        /*
        if (lmd < pow(10,-12))
        {
            NewtonConstantEnergy(x,E0,0,k0,deltaT,log,iterations,tmax,integrator);
        }
        */

        int nonLinE3 = round(nonLin * 1000000);
        if (lmd <= pow(10, -12) && nonLinE3 % 10000 == 0)
        {
            xb = x;
            savePeriodicOrbitReducedAGPGradient(x, k0, E0init, lmd, deltaT, deltaTAGPGradient, iterations, tmax, saveFolder, integrator);
        }

        nonLin = nonLin + deltaNonLin;
        x.nonLin = nonLin;
    }

    system("pause");
    // return 1;
}