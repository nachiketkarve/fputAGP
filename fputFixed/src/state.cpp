#include "state.hpp"
#include "couplings.hpp"
#include "standards.hpp"
#include "modeAmplitude.hpp"

State::State(int n, std::string modelName = "alpha") : N(n), q(n + 2), p(n + 2), nonLin(0), model(modelName)
{
    q.setZero();
    p.setZero();
    if (model != "alpha" && model != "beta" && model != "toda")
    {
        throw std::invalid_argument("Invalid Model");
    }
}

State::State(const State &x)
{
    q = x.q;
    p = x.p;
    N = x.N;
    nonLin = x.nonLin;
    model = x.model;
}

// Computes the acceleration at each site on the lattice
void State::eom(Vector &ddq, const Vector &Q)
{
    ddq.setZero();
    if (model == "alpha")
    {
        for (int i = 1; i <= N; i++)
        {
            ddq(i) = Q(i + 1) + Q(i - 1) - 2.0 * Q(i) + nonLin * ((Q(i + 1) - Q(i)) * (Q(i + 1) - Q(i)) - (Q(i) - Q(i - 1)) * (Q(i) - Q(i - 1)));
        }
    }
    else if (model == "beta")
    {
        for (int i = 1; i <= N; i++)
        {
            ddq(i) = Q(i + 1) + Q(i - 1) - 2.0 * Q(i) + nonLin * ((Q(i + 1) - Q(i)) * (Q(i + 1) - Q(i)) * (Q(i + 1) - Q(i)) - (Q(i) - Q(i - 1)) * (Q(i) - Q(i - 1)) * (Q(i) - Q(i - 1)));
        }
    }
    else if (model == "toda")
    {
        double V0 = 1/(4.0*nonLin*nonLin);
        double lmd0 = 2.0*nonLin;
        for (int i = 1; i <= N; i++)
        {
            ddq(i) = V0*lmd0*std::exp(lmd0*(Q(i+1)-Q(i))) - V0*lmd0*std::exp(lmd0*(Q(i)-Q(i-1)));
        }
    }
}

// Computes the acceleration at each site on the lattice and of the tangent vectors
void State::eomTangent(Vector &ddq, Matrix &ddqTangent, const Vector &Q, const Matrix &QTangent)
{
    ddq.setZero();
    ddqTangent.setZero();
    for (int i = 1; i <= N; i++)
    {
        if (model == "alpha")
        {
            ddq(i) = Q(i + 1) + Q(i - 1) - 2.0 * Q(i) + nonLin * (pow(Q(i + 1) - Q(i), 2) - pow(Q(i) - Q(i - 1), 2));
            for (int j = 0; j < QTangent.cols(); j++)
            {
                ddqTangent(i, j) = QTangent(i + 1, j) + QTangent(i - 1, j) - 2.0 * QTangent(i, j) + 2.0 * nonLin * ((Q(i + 1) - Q(i)) * (QTangent(i + 1, j) - QTangent(i, j)) - (Q(i) - Q(i - 1)) * (QTangent(i, j) - QTangent(i - 1, j)));
            }
        }
        else if (model == "beta")
        {
            ddq(i) = Q(i + 1) + Q(i - 1) - 2.0 * Q(i) + nonLin * (pow(Q(i + 1) - Q(i), 3) - pow(Q(i) - Q(i - 1), 3));
            for (int j = 0; j < QTangent.cols(); j++)
            {
                ddqTangent(i, j) = QTangent(i + 1, j) + QTangent(i - 1, j) - 2.0 * QTangent(i, j) + 3.0 * nonLin * (pow(Q(i + 1) - Q(i), 2) * (QTangent(i + 1, j) - QTangent(i, j)) - pow(Q(i) - Q(i - 1), 2) * (QTangent(i, j) - QTangent(i - 1, j)));
            }
        }
        else if (model == "toda")
        {
            double V0 = 1/(4.0*nonLin*nonLin);
            double lmd0 = 2.0*nonLin;
            ddq(i) = V0*lmd0*std::exp(lmd0*(Q(i+1)-Q(i))) - V0*lmd0*std::exp(lmd0*(Q(i)-Q(i-1)));
            for (int j = 0; j < QTangent.cols(); j++)
            {
                ddqTangent(i, j) = std::exp(lmd0*(Q(i+1)-Q(i)))*(QTangent(i+1,j)-QTangent(i,j)) - std::exp(lmd0*(Q(i)-Q(i-1)))*(QTangent(i,j)-QTangent(i-1,j));
            }
        }
    }
}

// Computes the total energy of the system
double State::totalEnergy()
{
    double Energy = 0.0;
    if (model == "alpha")
    {
        for (int i = 0; i < N + 2; i++)
        {
            Energy = Energy + 0.5 * p(i) * p(i);
            if (i < N + 1)
            {
                Energy = Energy + 0.5 * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) + nonLin / 3.0 * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) * (q(i + 1) - q(i));
            }
        }
    }
    else if (model == "beta")
    {
        for (int i = 0; i < N + 2; i++)
        {
            Energy = Energy + 0.5 * p(i) * p(i);
            if (i < N + 1)
            {
                Energy = Energy + 0.5 * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) + nonLin / 4.0 * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) * (q(i + 1) - q(i));
            }
        }
    }

    else if (model == "toda")
    {
        double V0 = 1/(4.0*nonLin*nonLin);
        double lmd0 = 2.0*nonLin;
        for (int i = 0; i < N + 2; i++)
        {
            Energy = Energy + 0.5 * pow(p(i), 2);
            if (i < N + 1)
            {
                Energy = Energy + V0*(std::exp(lmd0*(q(i+1)-q(i))) - 1.0 - lmd0*(q(i+1)-q(i)));
            }
        }
    }

    return Energy;
};

// Computes the interaction potential relevant to the AGP of the system
double State::interactionPotential()
{

    double Energy = 0.0;
    if (model == "alpha")
    {
        for (int i = 0; i < N + 2; i++)
        {
            if (i < N + 1)
            {
                Energy = Energy + 1.0 / 3.0 * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) * (q(i + 1) - q(i));
            }
        }
    }
    else if (model == "beta")
    {
        for (int i = 0; i < N + 2; i++)
        {
            if (i < N + 1)
            {
                Energy = Energy + 1.0 / 4.0 * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) * (q(i + 1) - q(i)) * (q(i + 1) - q(i));
            }
        }
    }

    else if (model == "toda")
    {
        double V0 = 1/(4.0*nonLin*nonLin);
        double lmd0 = 2.0*nonLin;
        for (int i = 0; i < N + 2; i++)
        {
            if (i < N + 1)
            {
                Energy = Energy + 1.0/(2.0*nonLin*nonLin*nonLin) * (1.0 - std::exp(lmd0*(q(i+1)-q(i)))) + (q(i+1)-q(i))/(2.0*nonLin*nonLin) * (1.0 + std::exp(lmd0*(q(i+1)-q(i))));
            }
        }
    }

    return Energy;
};

// Calculates the phase space gradient of the interaction potential and stores it in the argument
void State::potentialGradient(Vector &grad)
{
    if (grad.size() != 2 * N)
    {
        throw std::invalid_argument("Incompatible Gradient Vector");
    }

    grad.setZero();

    for (int i = 1; i < N+1; i++)
    {
        if (model == "alpha")
        {
            grad(i-1) =  - (q(i+1) - q(i)) * (q(i+1) - q(i)) + (q(i) - q(i-1)) * (q(i) - q(i-1));
        } else if (model == "beta")
        {
            grad(i-1) =  - (q(i+1) - q(i)) * (q(i+1) - q(i)) * (q(i+1) - q(i)) + (q(i) - q(i-1)) * (q(i) - q(i-1)) * (q(i) - q(i-1));
        } else if (model == "toda")
        {
            grad(i-1) = std::exp(2.0*nonLin*(q(i)-q(i-1)))*((q(i)-q(i-1))/nonLin - 1.0/(2.0*nonLin*nonLin)) - std::exp(2.0*nonLin*(q(i+1)-q(i)))*((q(i+1)-q(i))/nonLin - 1.0/(2.0*nonLin*nonLin));
        }
    }
};

// Integrates the equations of motion over one time step using the leap frog integrator
void State::leapFrog(double deltaT)
{
    Vector a1(N + 2);
    Vector a2(N + 2);
    eom(a1,q);
    q = q + deltaT*p + 0.5*deltaT*deltaT*a1;
    eom(a2,q);
    p = p + 0.5*deltaT*(a1 + a2);
}

// Integrates the equations of motion of the lattice and tangent vectors for one time step using leap frog
void State::leapFrogTangent(Matrix &qTangent, Matrix &pTangent, double deltaT)
{
    if (qTangent.rows() != N + 2 || pTangent.rows() != N + 2)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    Vector a1(N + 2);
    Vector a2(N + 2);

    Matrix a1Tangent(qTangent.rows(), qTangent.cols());
    Matrix a2Tangent(qTangent.rows(), qTangent.cols());

    eomTangent(a1, a1Tangent, q, qTangent);
    q = q + deltaT*p + 0.5*deltaT*deltaT*a1;
    qTangent = qTangent + deltaT*pTangent + 0.5*deltaT*deltaT*a1Tangent;
    eomTangent(a2, a2Tangent, q, qTangent);
    p = p + 0.5*deltaT*(a1 + a2);
    pTangent = pTangent + 0.5*deltaT*(a1Tangent + a2Tangent);

}

// Integrates the equations of motion over one time step using the 4th order symplectic Runge-Kutta integrator
void State::rkn4(double deltaT)
{
    Vector pPrev(N + 2);
    Vector EOM1(N + 2);
    Vector EOM2(N + 2);
    Vector EOM3(N + 2);
    Vector EOM4(N + 2);
    Vector EOM5(N + 2);

    Vector Q1(N + 2);
    Vector Q2(N + 2);
    Vector Q3(N + 2);
    Vector Q4(N + 2);
    Vector Q5(N + 2);

    Q1 = q + deltaT * g1 * p;
    eom(EOM1, Q1);
    Q2 = q + deltaT * g2 * p + pow(deltaT, 2) * a21 * EOM1;
    eom(EOM2, Q2);
    Q3 = q + deltaT * g3 * p + pow(deltaT, 2) * a31 * EOM1 + pow(deltaT, 2) * a32 * EOM2;
    eom(EOM3, Q3);
    Q4 = q + deltaT * g4 * p + pow(deltaT, 2) * a41 * EOM1 + pow(deltaT, 2) * a42 * EOM2 + pow(deltaT, 2) * a43 * EOM3;
    eom(EOM4, Q4);
    Q5 = q + deltaT * g5 * p + pow(deltaT, 2) * a51 * EOM1 + pow(deltaT, 2) * a52 * EOM2 + pow(deltaT, 2) * a53 * EOM3 + pow(deltaT, 2) * a54 * EOM4;
    eom(EOM5, Q5);

    pPrev = p;
    p = p + deltaT * b1 * EOM1 + +deltaT * b2 * EOM2 + deltaT * b3 * EOM3 + deltaT * b4 * EOM4 + deltaT * b5 * EOM5;
    q = q + deltaT * pPrev + pow(deltaT, 2) * B1 * EOM1 + pow(deltaT, 2) * B2 * EOM2 + pow(deltaT, 2) * B3 * EOM3 + pow(deltaT, 2) * B4 * EOM4 + pow(deltaT, 2) * B5 * EOM5;
};

// Integrates the equations of motion of the lattice and tangent vectors for one time step using RK4
void State::rknTangent4(Matrix &qTangent, Matrix &pTangent, double deltaT)
{

    if (qTangent.rows() != N + 2 || pTangent.rows() != N + 2)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    Vector pPrev(N + 2);
    Vector EOM1(N + 2);
    Vector EOM2(N + 2);
    Vector EOM3(N + 2);
    Vector EOM4(N + 2);
    Vector EOM5(N + 2);

    Vector Q1(N + 2);
    Vector Q2(N + 2);
    Vector Q3(N + 2);
    Vector Q4(N + 2);
    Vector Q5(N + 2);

    Matrix pPreve(qTangent.rows(), qTangent.cols());
    Matrix EOM1e(qTangent.rows(), qTangent.cols());
    Matrix EOM2e(qTangent.rows(), qTangent.cols());
    Matrix EOM3e(qTangent.rows(), qTangent.cols());
    Matrix EOM4e(qTangent.rows(), qTangent.cols());
    Matrix EOM5e(qTangent.rows(), qTangent.cols());

    Matrix Q1e(qTangent.rows(), qTangent.cols());
    Matrix Q2e(qTangent.rows(), qTangent.cols());
    Matrix Q3e(qTangent.rows(), qTangent.cols());
    Matrix Q4e(qTangent.rows(), qTangent.cols());
    Matrix Q5e(qTangent.rows(), qTangent.cols());

    Q1 = q + deltaT * g1 * p;
    Q1e = qTangent + deltaT * g1 * pTangent;
    eomTangent(EOM1, EOM1e, Q1, Q1e);
    Q2 = q + deltaT * g2 * p + pow(deltaT, 2) * a21 * EOM1;
    Q2e = qTangent + deltaT * g2 * pTangent + pow(deltaT, 2) * a21 * EOM1e;
    eomTangent(EOM2, EOM2e, Q2, Q2e);
    Q3 = q + deltaT * g3 * p + pow(deltaT, 2) * a31 * EOM1 + pow(deltaT, 2) * a32 * EOM2;
    Q3e = qTangent + deltaT * g3 * pTangent + pow(deltaT, 2) * a31 * EOM1e + pow(deltaT, 2) * a32 * EOM2e;
    eomTangent(EOM3, EOM3e, Q3, Q3e);
    Q4 = q + deltaT * g4 * p + pow(deltaT, 2) * a41 * EOM1 + pow(deltaT, 2) * a42 * EOM2 + pow(deltaT, 2) * a43 * EOM3;
    Q4e = qTangent + deltaT * g4 * pTangent + pow(deltaT, 2) * a41 * EOM1e + pow(deltaT, 2) * a42 * EOM2e + pow(deltaT, 2) * a43 * EOM3e;
    eomTangent(EOM4, EOM4e, Q4, Q4e);
    Q5 = q + deltaT * g5 * p + pow(deltaT, 2) * a51 * EOM1 + pow(deltaT, 2) * a52 * EOM2 + pow(deltaT, 2) * a53 * EOM3 + pow(deltaT, 2) * a54 * EOM4;
    Q5e = qTangent + deltaT * g5 * pTangent + pow(deltaT, 2) * a51 * EOM1e + pow(deltaT, 2) * a52 * EOM2e + pow(deltaT, 2) * a53 * EOM3e + pow(deltaT, 2) * a54 * EOM4e;
    eomTangent(EOM5, EOM5e, Q5, Q5e);

    pPrev = p;
    pPreve = pTangent;
    p = p + deltaT * b1 * EOM1 + deltaT * b2 * EOM2 + deltaT * b3 * EOM3 + deltaT * b4 * EOM4 + deltaT * b5 * EOM5;
    pTangent = pTangent + deltaT * b1 * EOM1e + deltaT * b2 * EOM2e + deltaT * b3 * EOM3e + deltaT * b4 * EOM4e + deltaT * b5 * EOM5e;
    q = q + deltaT * pPrev + pow(deltaT, 2) * B1 * EOM1 + pow(deltaT, 2) * B2 * EOM2 + pow(deltaT, 2) * B3 * EOM3 + pow(deltaT, 2) * B4 * EOM4 + pow(deltaT, 2) * B5 * EOM5;
    qTangent = qTangent + deltaT * pPreve + pow(deltaT, 2) * B1 * EOM1e + pow(deltaT, 2) * B2 * EOM2e + pow(deltaT, 2) * B3 * EOM3e + pow(deltaT, 2) * B4 * EOM4e + pow(deltaT, 2) * B5 * EOM5e;
};

// Integrates the equations of motion over one time step using the 4th order symplectic Runge-Kutta integrator
void State::rkn6(double deltaT)
{
    Vector pPrev(N + 2);
    Vector EOM1(N + 2);
    Vector EOM2(N + 2);
    Vector EOM3(N + 2);
    Vector EOM4(N + 2);
    Vector EOM5(N + 2);
    Vector EOM6(N + 2);
    Vector EOM7(N + 2);

    Vector Q1(N + 2);
    Vector Q2(N + 2);
    Vector Q3(N + 2);
    Vector Q4(N + 2);
    Vector Q5(N + 2);
    Vector Q6(N + 2);
    Vector Q7(N + 2);

    Q1 = q + deltaT * g61 * p;
    eom(EOM1, Q1);
    Q2 = q + deltaT * g62 * p + pow(deltaT, 2) * a621 * EOM1;
    eom(EOM2, Q2);
    Q3 = q + deltaT * g63 * p + pow(deltaT, 2) * a631 * EOM1 + pow(deltaT, 2) * a632 * EOM2;
    eom(EOM3, Q3);
    Q4 = q + deltaT * g64 * p + pow(deltaT, 2) * a641 * EOM1 + pow(deltaT, 2) * a642 * EOM2 + pow(deltaT, 2) * a643 * EOM3;
    eom(EOM4, Q4);
    Q5 = q + deltaT * g65 * p + pow(deltaT, 2) * a651 * EOM1 + pow(deltaT, 2) * a652 * EOM2 + pow(deltaT, 2) * a653 * EOM3 + pow(deltaT, 2) * a654 * EOM4;
    eom(EOM5, Q5);
    Q6 = q + deltaT * g66 * p + pow(deltaT, 2) * a661 * EOM1 + pow(deltaT, 2) * a662 * EOM2 + pow(deltaT, 2) * a663 * EOM3 + pow(deltaT, 2) * a664 * EOM4 + pow(deltaT, 2) * a665 * EOM5;
    eom(EOM6, Q6);
    Q7 = q + deltaT * g67 * p + pow(deltaT, 2) * a671 * EOM1 + pow(deltaT, 2) * a672 * EOM2 + pow(deltaT, 2) * a673 * EOM3 + pow(deltaT, 2) * a674 * EOM4 + pow(deltaT, 2) * a675 * EOM5 + pow(deltaT, 2) * a676 * EOM6;
    eom(EOM7, Q7);

    pPrev = p;
    p = p + deltaT * b61 * EOM1 + +deltaT * b62 * EOM2 + deltaT * b63 * EOM3 + deltaT * b64 * EOM4 + deltaT * b65 * EOM5 + deltaT * b66 * EOM6 + deltaT * b67 * EOM7;
    q = q + deltaT * pPrev + pow(deltaT, 2) * B61 * EOM1 + pow(deltaT, 2) * B62 * EOM2 + pow(deltaT, 2) * B63 * EOM3 + pow(deltaT, 2) * B64 * EOM4 + pow(deltaT, 2) * B65 * EOM5 + pow(deltaT, 2) * B66 * EOM6 + pow(deltaT, 2) * B67 * EOM7;
};

// Integrates the equations of motion of the lattice and tangent vectors for one time step using RK6
void State::rknTangent6(Matrix &qTangent, Matrix &pTangent, double deltaT)
{

    if (qTangent.rows() != N + 2 || pTangent.rows() != N + 2)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    Vector pPrev(N + 2);
    Vector EOM1(N + 2);
    Vector EOM2(N + 2);
    Vector EOM3(N + 2);
    Vector EOM4(N + 2);
    Vector EOM5(N + 2);
    Vector EOM6(N + 2);
    Vector EOM7(N + 2);

    Vector Q1(N + 2);
    Vector Q2(N + 2);
    Vector Q3(N + 2);
    Vector Q4(N + 2);
    Vector Q5(N + 2);
    Vector Q6(N + 2);
    Vector Q7(N + 2);

    Matrix pPreve(qTangent.rows(), qTangent.cols());
    Matrix EOM1e(qTangent.rows(), qTangent.cols());
    Matrix EOM2e(qTangent.rows(), qTangent.cols());
    Matrix EOM3e(qTangent.rows(), qTangent.cols());
    Matrix EOM4e(qTangent.rows(), qTangent.cols());
    Matrix EOM5e(qTangent.rows(), qTangent.cols());
    Matrix EOM6e(qTangent.rows(), qTangent.cols());
    Matrix EOM7e(qTangent.rows(), qTangent.cols());

    Matrix Q1e(qTangent.rows(), qTangent.cols());
    Matrix Q2e(qTangent.rows(), qTangent.cols());
    Matrix Q3e(qTangent.rows(), qTangent.cols());
    Matrix Q4e(qTangent.rows(), qTangent.cols());
    Matrix Q5e(qTangent.rows(), qTangent.cols());
    Matrix Q6e(qTangent.rows(), qTangent.cols());
    Matrix Q7e(qTangent.rows(), qTangent.cols());

    Q1 = q + deltaT * g61 * p;
    Q1e = qTangent + deltaT * g61 * pTangent;
    eomTangent(EOM1, EOM1e, Q1, Q1e);
    Q2 = q + deltaT * g62 * p + pow(deltaT, 2) * a621 * EOM1;
    Q2e = qTangent + deltaT * g62 * pTangent + pow(deltaT, 2) * a621 * EOM1e;
    eomTangent(EOM2, EOM2e, Q2, Q2e);
    Q3 = q + deltaT * g63 * p + pow(deltaT, 2) * a631 * EOM1 + pow(deltaT, 2) * a632 * EOM2;
    Q3e = qTangent + deltaT * g63 * pTangent + pow(deltaT, 2) * a631 * EOM1e + pow(deltaT, 2) * a632 * EOM2e;
    eomTangent(EOM3, EOM3e, Q3, Q3e);
    Q4 = q + deltaT * g64 * p + pow(deltaT, 2) * a641 * EOM1 + pow(deltaT, 2) * a642 * EOM2 + pow(deltaT, 2) * a643 * EOM3;
    Q4e = qTangent + deltaT * g64 * pTangent + pow(deltaT, 2) * a641 * EOM1e + pow(deltaT, 2) * a642 * EOM2e + pow(deltaT, 2) * a643 * EOM3e;
    eomTangent(EOM4, EOM4e, Q4, Q4e);
    Q5 = q + deltaT * g65 * p + pow(deltaT, 2) * a651 * EOM1 + pow(deltaT, 2) * a652 * EOM2 + pow(deltaT, 2) * a653 * EOM3 + pow(deltaT, 2) * a654 * EOM4;
    Q5e = qTangent + deltaT * g65 * pTangent + pow(deltaT, 2) * a651 * EOM1e + pow(deltaT, 2) * a652 * EOM2e + pow(deltaT, 2) * a653 * EOM3e + pow(deltaT, 2) * a654 * EOM4e;
    eomTangent(EOM5, EOM5e, Q5, Q5e);
    Q6 = q + deltaT * g66 * p + pow(deltaT, 2) * a661 * EOM1 + pow(deltaT, 2) * a662 * EOM2 + pow(deltaT, 2) * a663 * EOM3 + pow(deltaT, 2) * a664 * EOM4 + pow(deltaT, 2) * a665 * EOM5;
    Q6e = qTangent + deltaT * g66 * pTangent + pow(deltaT, 2) * a661 * EOM1e + pow(deltaT, 2) * a662 * EOM2e + pow(deltaT, 2) * a663 * EOM3e + pow(deltaT, 2) * a664 * EOM4e + pow(deltaT, 2) * a665 * EOM5e;
    eomTangent(EOM6, EOM6e, Q6, Q6e);
    Q7 = q + deltaT * g67 * p + pow(deltaT, 2) * a671 * EOM1 + pow(deltaT, 2) * a672 * EOM2 + pow(deltaT, 2) * a673 * EOM3 + pow(deltaT, 2) * a674 * EOM4 + pow(deltaT, 2) * a675 * EOM5 + pow(deltaT, 2) * a676 * EOM6;
    Q7e = qTangent + deltaT * g67 * pTangent + pow(deltaT, 2) * a671 * EOM1e + pow(deltaT, 2) * a672 * EOM2e + pow(deltaT, 2) * a673 * EOM3e + pow(deltaT, 2) * a674 * EOM4e + pow(deltaT, 2) * a675 * EOM5e + pow(deltaT, 2) * a676 * EOM6e;
    eomTangent(EOM7, EOM7e, Q7, Q7e);

    pPrev = p;
    pPreve = pTangent;
    p = p + deltaT * b61 * EOM1 + deltaT * b62 * EOM2 + deltaT * b63 * EOM3 + deltaT * b64 * EOM4 + deltaT * b65 * EOM5 + deltaT * b66 * EOM6 + deltaT * b67 * EOM7;
    pTangent = pTangent + deltaT * b61 * EOM1e + deltaT * b62 * EOM2e + deltaT * b63 * EOM3e + deltaT * b64 * EOM4e + deltaT * b65 * EOM5e + deltaT * b66 * EOM6e + deltaT * b67 * EOM7e;
    q = q + deltaT * pPrev + pow(deltaT, 2) * B61 * EOM1 + pow(deltaT, 2) * B62 * EOM2 + pow(deltaT, 2) * B63 * EOM3 + pow(deltaT, 2) * B64 * EOM4 + pow(deltaT, 2) * B65 * EOM5 + pow(deltaT, 2) * B66 * EOM6 + pow(deltaT, 2) * B67 * EOM7;
    qTangent = qTangent + deltaT * pPreve + pow(deltaT, 2) * B61 * EOM1e + pow(deltaT, 2) * B62 * EOM2e + pow(deltaT, 2) * B63 * EOM3e + pow(deltaT, 2) * B64 * EOM4e + pow(deltaT, 2) * B65 * EOM5e + pow(deltaT, 2) * B66 * EOM6e + pow(deltaT, 2) * B67 * EOM7e;
};

// Evolve the system over a given period of time
void State::Evolve(double tmax, double deltaT, std::string integrator = "rk4")
{
    for (int i = 0; i < int(tmax / deltaT); i++)
    {
        if (integrator == "leapFrog")
        {
            leapFrog(deltaT);
        } else if (integrator == "rk6")
        {
            rkn6(deltaT);
        } else 
        {
            rkn4(deltaT);
        }
    }
};

// Evolve the system until it crosses the Poincare section described by the normal mode momentum of the seed mode
double State::EvolvePoincare(int seedMode, double deltaT, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4")
{
    double period = 0.0;
    int count = 0;
    Matrix FourierComponentsSeedMode(1, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        FourierComponentsSeedMode(0, i) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(seedMode) * pi / double(N + 1));
    }

    double QSeedModeInit = (FourierComponentsSeedMode * q)(0, 0);
    double PSeedModeInit = (FourierComponentsSeedMode * p)(0, 0);

    double QSeedMode = QSeedModeInit;
    double PSeedMode = PSeedModeInit;
    double PSeedModePrev = PSeedModeInit;

    for (int i = 0; i < int(tmax / deltaT); i++)
    {
        PSeedModePrev = PSeedMode;
        if (integrator == "leapFrog")
        {
            leapFrog(deltaT);
        } else if (integrator == "rk6")
        {
            rkn6(deltaT);
        } else 
        {
            rkn4(deltaT);
        }
        period = period + deltaT;
        QSeedMode = (FourierComponentsSeedMode * q)(0, 0);
        PSeedMode = (FourierComponentsSeedMode * p)(0, 0);

        if (i > 1 && sgn(QSeedMode) == sgn(QSeedModeInit) && sgn(PSeedMode - PSeedModeInit) != sgn(PSeedModePrev - PSeedModeInit))
        {
            count = count + 1;
            /*
            double tBack = -(PSeedMode - PSeedModeInit)/(PSeedMode - PSeedModePrev)*deltaT;
            period = period + tBack;
            if (integrator == "leapFrog")
            {
                leapFrog(tBack);
            } else if (integrator == "rk6")
            {
                rkn6(tBack);
            } else 
            {
                rkn4(tBack);
            }
            */
        }
        if (count >= iterations)
        {
            break;
        }
    }

    return period;
};

// Evolve the system and the tangent vectors until it crosses the Poincare section described by the normal mode momentum
// of the seed mode
double State::EvolveTangentPoincare(Matrix &qTangent, Matrix &pTangent, int seedMode, double deltaT, int iterations = 1, double tmax = 1000.0, std::string integrator = "rk4")
{
    if (qTangent.rows() != N + 2 || pTangent.rows() != N + 2)
    {
        throw std::invalid_argument("Incompatible Tangent Vectors");
    }

    double period = 0.0;
    int count = 0;
    Matrix FourierComponentsSeedMode(1, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        FourierComponentsSeedMode(0, i) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(seedMode) * pi / double(N + 1));
    }

    double QSeedModeInit = (FourierComponentsSeedMode * q)(0, 0);
    double PSeedModeInit = (FourierComponentsSeedMode * p)(0, 0);

    double QSeedMode = QSeedModeInit;
    double PSeedMode = PSeedModeInit;
    double PSeedModePrev = PSeedModeInit;

    for (int i = 0; i < int(tmax / deltaT); i++)
    {
        PSeedModePrev = PSeedMode;
        if (integrator == "leapFrog")
        {
            leapFrogTangent(qTangent,pTangent,deltaT);
        } else if (integrator == "rk6")
        {
            rknTangent6(qTangent,pTangent,deltaT);
        } else 
        {
            rknTangent4(qTangent,pTangent,deltaT);
        }
        period = period + deltaT;
        QSeedMode = (FourierComponentsSeedMode * q)(0, 0);
        PSeedMode = (FourierComponentsSeedMode * p)(0, 0);

        if (i > 1 && sgn(QSeedMode) == sgn(QSeedModeInit) && sgn(PSeedMode - PSeedModeInit) != sgn(PSeedModePrev - PSeedModeInit))
        {
            count = count + 1;
            /*
            double tBack = -(PSeedMode - PSeedModeInit)/(PSeedMode - PSeedModePrev)*deltaT;
            period = period + tBack;
            if (integrator == "leapFrog")
            {
                leapFrog(tBack);
            } else if (integrator == "rk6")
            {
                rkn6(tBack);
            } else 
            {
                rkn4(tBack);
            }
            */
        }
        if (count >= iterations)
        {
            break;
        }
    }

    return period;
};

// Initializes the system into a single normal mode with a given energy
void State::initializeNormalMode(double Energy, int mode)
{
    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Q(N+2);
    Vector P(N+2);
    Q.setZero();
    P.setZero();
    if (model == "alpha" || model == "toda")
    {
        modeAmplitudeAlpha(Q,P,mode,Energy,nonLin);
    } else if (model == "beta")
    {
        modeAmplitudeBeta(Q,P,mode,Energy,nonLin);
    }

    p.setZero();
    q = FourierComponents * Q;
};