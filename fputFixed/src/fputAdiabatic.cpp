#include "fputFixed.hpp"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

using json = nlohmann::json;

int main()
{
    std::ifstream file("_initAdiabatic.json");
    json initData = json::parse(file);

    std::string model = initData["model"];
    int N = initData["N"];
    double E0 = initData["Energy"];
    int k0 = initData["seedMode"];
    double nonLin = 0.0;
    double deltaT = initData["deltaT"];
    double tmax = initData["tmax"];
    double nonLinRate = initData["nonLinRate"];
    bool log = initData["log"];
    std::string integrator = initData["integrator"];
    std::string saveFolder = initData["saveFolder"];


    std::ofstream saveFile;
    std::string FileName;
    if (model == "alpha")
    {
        FileName = saveFolder + "alphaAd-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + ".csv";
    } else if (model == "toda")
    {
        FileName = saveFolder + "todaAd-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + ".csv";
    } else if (model == "beta")
    {
        FileName = saveFolder + "betaAd-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + ".csv";
    }

    saveFile.open(FileName);
    saveFile << std::setprecision(15);

    saveFile << "Time" << "," << "NonLin"<< "," << "Energy" << "\n";

    State x(N, model);
    x.nonLin = nonLin;
    x.initializeNormalMode(E0, k0);

    int iSkip = int(tmax/(100*deltaT));


    int i = 0;
    double tCurr = 0;

    while (tCurr <= tmax)
    {
        if (i%iSkip == 0)
        {
            saveFile << tCurr << "," << nonLin << "," << x.totalEnergy() << "\n";
        }
        if (integrator == "leapFrog")
        {
            x.leapFrog(deltaT);
        } else if (integrator == "rk6")
        {
            x.rkn6(deltaT);
        } else
        {
            x.rkn4(deltaT);
        }

        tCurr = tCurr + deltaT;
        nonLin = nonLin + nonLinRate * deltaT;
        x.nonLin = nonLin;
        i = i + 1;

        
    }

    std::cout << x.q;

    saveFile.close();
    FileName = saveFolder + "endstate.csv";
    saveFile.open(FileName);
    saveFile << std::setprecision(15);
    saveFile << "q" << "," << "p" << "\n";
    for (int i = 0; i < N+2; i++)
    {
        saveFile << x.q(i) << "," << x.p(i) << "\n";    
    }

    saveFile.close();
    system("pause");
    // return 1;
}