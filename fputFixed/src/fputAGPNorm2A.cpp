#include "fputFixed.hpp"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    int N = 8;
    double Q1 = 1;
    double Q2 = 1;
    double deltaQ1 = 0.01;
    double deltaQ2 = 0.01;
    double nonLin = 0.0;

    if (argc != 7)
    {
        std::cout << argc << "\n";
        return -1;
    }

    N = int(std::atof(argv[1]));
    Q1 = std::atof(argv[2]);
    deltaQ1 = std::atof(argv[3]);
    Q2 = std::atof(argv[4]);
    deltaQ2 = std::atof(argv[5]);
    nonLin = std::atof(argv[6]);

    Vector Q1s(int(std::abs(Q1)/deltaQ1));
    Vector Q2s(int(std::abs(Q2)/deltaQ2));
    Q1s.setZero();
    Q2s.setZero();

    for (int i = 1; i < Q1s.size(); i++)
    {
        Q1s(i) = Q1s(i-1) + deltaQ1;
    }

    for (int i = 1; i < Q2s.size(); i++)
    {
        Q2s(i) = Q2s(i-1) + deltaQ2;
    }

    Matrix AGPnorms(int(std::abs(Q1)/deltaQ1),int(std::abs(Q2)/deltaQ2));
    AGPnorms.setZero();
    
    std::ifstream dataFile("_params.json");
    json initData = json::parse(dataFile);

    std::string model = initData["model"];
    double deltaT = initData["deltaT"];
    double tmax = initData["tmax"];
    int dataPoints = initData["dataPoints"];
    std::string integrator = initData["integrator"];
    std::string saveFolder = initData["saveFolder"];

    std::string FileName;

    if (model == "alpha")
    {
        FileName = "alphaInit-N" + std::to_string(N) + "-Q1" + std::to_string(Q1) + "-Q2" + std::to_string(Q2) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "toda")
    {
        FileName = "todaInit-N" + std::to_string(N) + "-Q1" + std::to_string(Q1) + "-Q2" + std::to_string(Q2) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "beta")
    {
        FileName = "betaInit-N" + std::to_string(N) + "-Q1" + std::to_string(Q1) + "-Q2" + std::to_string(Q2) + "-B" + std::to_string(nonLin) + ".csv";
    }

    std::ofstream file;
    file.open(FileName);
    file << std::setprecision(15);

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }


    for (int i1 = 0; i1 < Q1s.size(); i1++)
    {
        for (int i2 = 0; i2 < Q2s.size(); i2++)
        {
            State x(N, model);
            x.nonLin = nonLin;
            Vector Q(N+2);
            Q.setZero();
            Q(1) = Q1s(i1);
            Q(2) = Q2s(i2);
            x.q = FourierComponents * Q;

            int iStep = int(tmax/(deltaT*dataPoints));
            
            double currT = 0;
            double dHInt = 0;
            double dH = 0;
            double dHPrev = 0;

            double sumX = 0;
            double sumY = 0;
            double sumX2 = 0;
            double sumXY = 0;
            double slope = 0;

            Vector AGP(int(tmax/deltaT));
            Vector dHIntArray(int(tmax/deltaT));
            Vector timeArray(int(tmax/deltaT));

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
            }

            slope = (round(tmax/deltaT)*sumXY - sumX*sumY)/(round(tmax/deltaT)*sumX2 - sumX*sumX);
            double AGPmean = 0;
            double AGP2mean = 0;
            for (int j = 0; j < AGP.size(); j = j + 1)
            {
                AGP(j) = -dHIntArray(j) + timeArray(j)*slope;
                AGPmean = AGPmean + AGP(j)/double(AGP.size());
                AGP2mean = AGP2mean + AGP(j)*AGP(j)/double(AGP.size());
            }
            AGPnorms(i1,i2) = AGP2mean - AGPmean*AGPmean;
        }
    }

    for (int i1 = 0; i1 < Q1s.size(); i1++)
    {
        for (int i2 = 0; i2 < Q2s.size(); i2++)
        {
            file << AGPnorms(i1,i2);
            if (i2 < Q2s.size() - 1)
            {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();

    system("pause");
    // return 1;
}