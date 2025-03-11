#include "fputFixed.hpp"
#include "nlohmann/json.hpp"
#include <random>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    int N = 8;
    double E0 = 1;
    int k0 = 1;
    double nonLin = 0.0;
    int seed = 1;

    if (argc != 7)
    {
        std::cout << argc << "\n";
        return -1;
    }

    N = int(std::atof(argv[1]));
    E0 = std::atof(argv[2]);
    k0 = std::atof(argv[3]);
    double nonLinNum = std::atof(argv[4]);
    double nonLinDen = std::atof(argv[5]);
    seed = int(std::atof(argv[6]));
    nonLin = nonLinNum/nonLinDen;

    std::mt19937                        generator(seed);
    std::uniform_real_distribution<double>  distr(0,2.0*pi);
    
    std::ifstream dataFile("_params.json");
    json initData = json::parse(dataFile);

    std::string model = initData["model"];
    double deltaT = initData["deltaT"];
    double tmax = initData["tmax"];
    int dataPoints = initData["dataPoints"];
    std::string integrator = initData["integrator"];
    std::string saveFolder = initData["saveFolder"];
    int averages = initData["averages"];

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

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
    file << "Time" << "," << "AGPnorm" << "," << "Energy" << "," << "Entropy";

    Vector dHIntArray(int(tmax/deltaT));
    Vector timeArray(int(tmax/deltaT));
    Vector AGPnorm(int(tmax/deltaT));
    Vector AGP(int(tmax/deltaT));
    Vector En(int(tmax/deltaT));
    Vector Entropy(int(tmax/deltaT));

    AGPnorm.setZero();
    En.setZero();
    Entropy.setZero();

    for (int average = 0; average < averages; average++)
    {
        State x(N, model);
        x.nonLin = nonLin;

        double theta = distr(generator);

        Vector Q(N+2);
        Vector P(N+2);
        Q.setZero();
        P.setZero();
        Q(k0) = std::sqrt(2.0*E0)*std::cos(theta)/freq(k0,N);
        P(k0) = std::sqrt(2.0*E0)*std::sin(theta);
        if (x.model == "beta" && nonLin*CouplingBeta(k0,k0,k0,k0,N)*Q(k0)*Q(k0)*Q(k0)*Q(k0) != 0)
        {
            Vector coeff(5);
            coeff.setZero();
            coeff(0) = -E0;
            coeff(2) = 0.5*P(k0)*P(k0) + 0.5*freq(k0,N)*freq(k0,N)*Q(k0)*Q(k0);
            coeff(4) = 3.0/4.0*nonLin*CouplingBeta(k0,k0,k0,k0,N)*Q(k0)*Q(k0)*Q(k0)*Q(k0);
            Eigen::PolynomialSolver<double, 4> solver;
            solver.compute(coeff);
            bool realEx;
            double scale = solver.greatestRealRoot(realEx, 0.000001);
            Q(k0) = scale*Q(k0);
            P(k0) = scale*P(k0);
        }
        x.q = FourierComponents * Q;
        x.p = FourierComponents * P;

        double currT = 0;
        double dHInt = 0;
        double dH = 0;
        double dHPrev = 0;

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

                AGPnorm(i) = AGPnorm(i) + (AGP2mean - AGPmean * AGPmean)/averages;
                En(i) = En(i) + x.totalEnergy()/averages;

                Vector Q(N+2);
                Vector P(N+2);
                Q = FourierComponents * x.q;
                P = FourierComponents * x.p;
                
                double entr = 0;
                for (int l = 1; l < N+1; l++)
                {
                    double enMode = 0.5*P(l)*P(l) + 0.5*freq(l,N)*freq(l,N)*Q(l)*Q(l);
                    if (enMode > pow(10,-20))
                    {
                        entr = entr - enMode/E0 * std::log(enMode/E0)/std::log(N);
                    }
                }

                Entropy(i) = Entropy(i) + entr/averages;


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
    }

    for (int i = 0; i < int(tmax/deltaT); i = i + iStep)
    {
        file << "\n" << timeArray(i) << "," << AGPnorm(i) << "," << En(i) << "," << Entropy(i);
    }

    file.close();

    system("pause");
    // return 1;
}