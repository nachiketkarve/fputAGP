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
        FileName = "alphaInitRev-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "toda")
    {
        FileName = "todaInitRev-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "beta")
    {
        FileName = "betaInitRev-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-B" + std::to_string(nonLin) + ".csv";
    }

    FileName = saveFolder + FileName;

    std::ofstream file;
    file.open(FileName);
    file << std::setprecision(15);
    file << "EInit" << "," << "EFinal";

    Vector EInit(N+2);
    Vector EFinal(N+2);

    EInit.setZero();
    EFinal.setZero();


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

        for (int i = 1; i < N+1; i++)
        {
            EInit(i) = EInit(i) + (0.5*P(i)*P(i) + 0.5*freq(i,N)*freq(i,N)*Q(i)*Q(i))/averages;
        }


        for (int i = 0; i < int(tmax/deltaT); i++)
        {

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

        x.p = -x.p;

        for (int i = 0; i < int(tmax/deltaT); i++)
        {

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

        Q = FourierComponents * x.q;
        P = FourierComponents * x.p;

        for (int i = 1; i < N+1; i++)
        {
            EFinal(i) = EFinal(i) + (0.5*P(i)*P(i) + 0.5*freq(i,N)*freq(i,N)*Q(i)*Q(i))/averages;
        }
    }

    for (int i = 0; i < N+2; i++)
    {
        file << "\n" << EInit(i) << "," << EFinal(i);
    }

    file.close();

    //system("pause");
    // return 1;
}