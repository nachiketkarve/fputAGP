#include "fputFixed.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main()
{
    std::ifstream file("_initBreatherAmplitude.json");
    json initData = json::parse(file);

    std::string model = initData["model"];
    int N = initData["N"];
    int k0 = initData["seedMode"];
    double nonLin = initData["nonLin"];
    double amplitudeStart = initData["amplitudeStart"];
    double amplitudeEnd = initData["amplitudeEnd"];
    double deltaAmplitude = initData["deltaAmplitude"];
    double deltaT = initData["deltaT"];
    double tmax = initData["tmax"];
    int iterations = initData["iterations"];
    bool log = initData["log"];
    std::string integrator = initData["integrator"];
    std::string saveFolder = initData["saveFolder"];

    double amplitude = amplitudeStart;
    int count = 0;
    double lmd = 1.0;

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q(N+2);
    Q.setZero();
    Q(k0) = amplitude;

    State x(N, model);
    x.nonLin = nonLin;
    x.q = FourierComponents * Q;

    State xb = x;

    while (amplitude <= amplitudeEnd)
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
            lmd = NewtonFixedAmplitude(x,amplitude,k0,deltaT,log,iterations,tmax,integrator);
            count = count + 1;
            if (lmd > pow(10, -1) || count > 15 || std::isnan(lmd) || std::isnan(x.totalEnergy()))
            {
                x = xb;
                break;
            }
        }

        int nonLinE3 = round(nonLin * 1000000);
        if (lmd <= pow(10, -12))
        {
            savePeriodicOrbitFixedAmplitude(x,k0,lmd,deltaT,iterations,tmax,saveFolder,integrator);
            xb = x;
        }

        amplitude = amplitude + deltaAmplitude;
    }

    system("pause");
    // return 1;
}
