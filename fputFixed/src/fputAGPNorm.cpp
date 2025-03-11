#include "fputFixed.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    int N = 8;
    double E0 = 1;
    int k0 = 1;
    double nonLin = 0.0;

    if (argc != 5)
    {
        std::cout << argc << "\n";
        return -1;
    }

    N = int(std::atof(argv[1]));
    E0 = std::atof(argv[2]);
    k0 = std::atof(argv[3]);
    nonLin = std::atof(argv[4]);
    
    std::ifstream dataFile("_params.json");
    json initData = json::parse(dataFile);

    std::string model = initData["model"];
    double deltaT = initData["deltaT"];
    double tmax = initData["tmax"];
    int dataPoints = initData["dataPoints"];
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

    State x(N, model);
    x.nonLin = nonLin;
    x.initializeNormalMode(E0, k0);

    int iStep = int(tmax/(deltaT*dataPoints));

    std::string FileName;

    if (model == "alpha")
    {
        FileName = "alphaInit-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "toda")
    {
        FileName = "todaInit-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "beta")
    {
        FileName = "betaInit-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-B" + std::to_string(nonLin) + ".csv";
    }

    std::ofstream file;
    file.open(FileName);
    file << std::setprecision(15);
    file << "Time" << "," << "AGPnorm" << "," << "Energy";

    double currT = 0;
    double dHInt = 0;
    double dH = 0;
    double dHPrev = 0;

    Vector dHIntArray(int(tmax/deltaT));
    Vector timeArray(int(tmax/deltaT));
    Vector AGPnorm(int(tmax/deltaT));
    Vector AGP(int(tmax/deltaT));

    double sumX = 0;
    double sumY = 0;
    double sumX2 = 0;
    double sumXY = 0;
    double slope = 0;

    for (int i = 0; i < int(tmax/deltaT); i++)
    {
        dH = x.interactionPotential();
        dHInt = dHInt + 0.5*deltaT*(dH+dHPrev);
        dHPrev = dH;
        dHIntArray(i) = dHInt;
        timeArray(i) = currT;

        sumX = sumX + currT;
        sumX2 = sumX2 + currT*currT;
        sumY = sumY + dHInt;
        sumXY = sumXY + currT*dHInt;
        if (i > 0)
        {
            slope = ((i+1)*sumXY - sumX*sumY)/((i+1)*sumX2 - sumX*sumX);
        }

        if (i % iStep == 0)
        {
            double AGPmean = 0;
            double AGP2mean = 0;
            int count = 0;
            for (int j = 0; j < i+1; j = j + iStep)
            {
                count = count + 1;
            }

            for (int j = 0; j < i+1; j = j + iStep)
            {
                AGP(j) = -dHIntArray(j) + timeArray(j)*slope;
                AGPmean = AGPmean + AGP(j)/double(int(i/iStep + 1));
                AGP2mean = AGP2mean + AGP(j)*AGP(j)/double(int(i/iStep + 1));
            }

            AGPnorm(i) = AGP2mean - AGPmean * AGPmean;

            file << "\n" << currT << "," << AGPnorm(i) << "," << x.totalEnergy();
        }
        if (integrator == "leapFrog")
        {
            x.leapFrog(deltaT);
        } else if (integrator == "rk4")
        {
            x.rkn4(deltaT);
        } else if (integrator == "rk6")
        {
            x.rkn6(deltaT);
        }
        currT = currT + deltaT;
        if (x.model == "beta")
        {
            Vector Q(N+2);
            Vector P(N+2);
            Q = FourierComponents * x.q;
            P = FourierComponents * x.p;

            for(int k = k0%2+1; k < N+2; k=k+2)
            {
                Q(k) = 0;
                P(k) = 0;
            }

            x.q = FourierComponents * Q;
            x.p = FourierComponents * P;
        }
    }

    system("pause");
    // return 1;
}