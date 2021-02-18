#include <cmath>
#include <vector>
#include "gsl/gsl_sf_gamma.h"

#include "../include/basis.h"
#include "../include/Integral.h"

#define COEFF_IGNORE_LIMIT 1e-20

double SIntegral(double ra[], double rb[], int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta)
{

    if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0)
    {
        return 0;
    }
    double zeta = alpha + beta;
    if (ax > 0)
    {
        double pi = (alpha * ra[0] + beta * rb[0]) / zeta;
        return (pi - ra[0]) * SIntegral(ra, rb, ax - 1, ay, az, bx, by, bz, alpha, beta) + 0.5 * bx / zeta * SIntegral(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha, beta) + 0.5 * (ax - 1) / zeta * SIntegral(ra, rb, ax - 2, ay, az, bx, by, bz, alpha, beta);
    }
    else if (ay > 0)
    {
        double pi = (alpha * ra[1] + beta * rb[1]) / zeta;
        return (pi - ra[1]) * SIntegral(ra, rb, ax, ay - 1, az, bx, by, bz, alpha, beta) + 0.5 * by / zeta * SIntegral(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha, beta) + 0.5 * (ay - 1) / zeta * SIntegral(ra, rb, ax, ay - 2, az, bx, by, bz, alpha, beta);
    }
    else if (az > 0)
    {
        double pi = (alpha * ra[2] + beta * rb[2]) / zeta;
        return (pi - ra[2]) * SIntegral(ra, rb, ax, ay, az - 1, bx, by, bz, alpha, beta) + 0.5 * bz / zeta * SIntegral(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha, beta) + 0.5 * (az - 1) / zeta * SIntegral(ra, rb, ax, ay, az - 2, bx, by, bz, alpha, beta);
    }
    else if (bx > 0)
    {
        double pi = (alpha * ra[0] + beta * rb[0]) / zeta;
        return (pi - rb[0]) * SIntegral(ra, rb, ax, ay, az, bx - 1, by, bz, alpha, beta) + 0.5 * ax / zeta * SIntegral(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha, beta) + 0.5 * (bx - 1) / zeta * SIntegral(ra, rb, ax, ay, az, bx - 2, by, bz, alpha, beta);
    }
    else if (by > 0)
    {
        double pi = (alpha * ra[1] + beta * rb[1]) / zeta;
        return (pi - rb[1]) * SIntegral(ra, rb, ax, ay, az, bx, by - 1, bz, alpha, beta) + 0.5 * ay / zeta * SIntegral(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha, beta) + 0.5 * (by - 1) / zeta * SIntegral(ra, rb, ax, ay, az, bx, by - 2, bz, alpha, beta);
    }
    else if (bz > 0)
    {
        double pi = (alpha * ra[2] + beta * rb[2]) / zeta;
        return (pi - rb[2]) * SIntegral(ra, rb, ax, ay, az, bx, by, bz - 1, alpha, beta) + 0.5 * az / zeta * SIntegral(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha, beta) + 0.5 * (bz - 1) / zeta * SIntegral(ra, rb, ax, ay, az, bx, by, bz - 2, alpha, beta);
    }
    else
    {
        // all zero
        double eta = alpha * beta / zeta;
        return pow(M_PI / zeta, 1.5) * exp(-eta * (pow(ra[0] - rb[0], 2) + pow(ra[1] - rb[1], 2) + pow(ra[2] - rb[2], 2)));
    }
}

double GTO_SIntegral(GTO &gto1, GTO &gto2)
{
    return gto1.coefficient * gto2.coefficient * SIntegral(gto1.cartesian, gto2.cartesian, gto1.ang[0], gto1.ang[1], gto1.ang[2], gto2.ang[0], gto2.ang[1], gto2.ang[2], gto1.orbital_exponent, gto2.orbital_exponent);
}

double Orbital_SIntegral(Orbital &ob1, Orbital &ob2)
{
    double result = 0;
    for (GTO &gtoa : ob1.component)
    {
        for (GTO &gtob : ob2.component)
        {
            result += GTO_SIntegral(gtoa, gtob);
        }
    }
    return result;
}

// derivative in direction dir for one gto
std::vector<GTO> GTO_Derivative(GTO &gto, int dir)
{
    std::vector<GTO> result;
    if (gto.ang[dir] == 0)
    {
        // only one part
        GTO g1 = gto;
        g1.ang[dir] += 1;
        g1.coefficient = g1.coefficient * -2 * g1.orbital_exponent;
        result.push_back(g1);
    }
    else
    {
        // two part
        GTO g1 = gto, g2 = gto;
        g1.ang[dir] += 1;
        g1.coefficient = g1.coefficient * -2.0 * g1.orbital_exponent;
        g2.coefficient *= g2.ang[dir];
        g2.ang[dir] -= 1;
        result.push_back(g1);
        result.push_back(g2);
    }
    return result;
}

/*
do second order derivative
0-x 1-y 2-z
*/
std::vector<GTO> GTO_SecondOrder_Derivative(GTO &gto, int dir)
{
    std::vector<GTO> result;
    auto d = GTO_Derivative(gto, dir);
    for (auto &dgto : d)
    {
        auto dd = GTO_Derivative(dgto, dir);
        result.insert(result.end(), dd.begin(), dd.end());
    }
    return result;
}

double GTO_Kinetic_Integral(GTO &LGTO, GTO &RGTO)
{
    double result = 0;
    for (int i = 0; i < 3; i++)
    {
        // for each direction
        // make second derivative to Right GTO
        std::vector<GTO> ddgto = GTO_SecondOrder_Derivative(RGTO, i);
        for (GTO &gto : ddgto)
        {
            // for each second derivative , calculate SIntegral
            result += GTO_SIntegral(LGTO, gto);
        }
    }
    // result need to multipy by -0.5
    result *= -0.5;
    return result;
}

double Orbital_Kinetic_Integral(Orbital &LOrbital, Orbital &ROrbital)
{
    double result = 0;
    for (auto &LGto : LOrbital.component)
    {
        for (auto &RGto : ROrbital.component)
        {
            result += GTO_Kinetic_Integral(LGto, RGto);
        }
    }
    return result;
}

// boys function
double Boys_Function(double x, int m)
{
    if (fabs(x) < 1e-8)
        return 1.0 / (1.0 + 2.0 * m);
    return 0.5 * pow(x, -0.5 - m) * (gsl_sf_gamma(m + 0.5) - gsl_sf_gamma_inc(0.5 + m, x));
}

double ZIntegral(const double ra[], const double rc[], const double rb[], const int anga[], const int angb[], double alpha, double beta, int m)
{
    for (int i = 0; i < 3; i++)
    {
        if (anga[i] < 0 || angb[i] < 0)
            return 0;
    }

    double zeta = alpha + beta;
    double xi = alpha * beta / zeta;
    double P[3];
    // AB and PC squared
    double AB = 0, PC = 0;

    for (int i = 0; i < 3; i++)
        P[i] = (alpha * ra[i] + beta * rb[i]) / zeta;

    for (int i = 0; i < 3; i++)
    {
        // cal AB squared and PC squared
        AB += pow(ra[i] - rb[i], 2);
        PC += pow(P[i] - rc[i], 2);
    }

    for (int i = 0; i < 3; i++)
    {
        if (anga[i] > 0)
        {
            // a minus 1 in one direction
            int am[3]{anga[0], anga[1], anga[2]};
            // a minus 2 in one direction
            int amm[3]{anga[0], anga[1], anga[2]};
            // b minus 1 in one direction
            int bm[3]{angb[0], angb[1], angb[2]};
            am[i] -= 1;
            amm[i] -= 2;
            bm[i] -= 1;
            return beta / zeta * (rb[i] - ra[i]) * ZIntegral(ra, rc, rb, am, angb, alpha, beta, m) + (rc[i] - ra[i] - beta / zeta * (rb[i] - ra[i])) * ZIntegral(ra, rc, rb, am, angb, alpha, beta, m + 1) + 0.5 * angb[i] / zeta * (ZIntegral(ra, rc, rb, am, bm, alpha, beta, m) - ZIntegral(ra, rc, rb, am, bm, alpha, beta, m + 1)) + 0.5 * am[i] / zeta * (ZIntegral(ra, rc, rb, amm, angb, alpha, beta, m) - ZIntegral(ra, rc, rb, amm, angb, alpha, beta, m + 1));
        }
        else if (angb[i] > 0)
        {
            // a minus 1 in one direction
            int am[3]{anga[0], anga[1], anga[2]};
            // b minus 2 in one direction
            int bmm[3]{angb[0], angb[1], angb[2]};
            // b minus 1 in one direction
            int bm[3]{angb[0], angb[1], angb[2]};
            am[i] -= 1;
            bm[i] -= 1;
            bmm[i] -= 2;

            return alpha / zeta * (ra[i] - rb[i]) * ZIntegral(ra, rc, rb, anga, bm, alpha, beta, m) + (rc[i] - rb[i] - alpha / zeta * (ra[i] - rb[i])) * ZIntegral(ra, rc, rb, anga, bm, alpha, beta, m + 1) + 0.5 * anga[i] / zeta * (ZIntegral(ra, rc, rb, am, bm, alpha, beta, m) - ZIntegral(ra, rc, rb, am, bm, alpha, beta, m + 1)) + 0.5 * bm[i] / zeta * (ZIntegral(ra, rc, rb, anga, bmm, alpha, beta, m) - ZIntegral(ra, rc, rb, anga, bmm, alpha, beta, m + 1));
        }
    }
    // all equal to zero, here must return init value
    return 2.0 * M_PI / zeta * exp(-xi * AB) * Boys_Function(zeta * PC, m);
}

double GTO_ZIntegral(GTO &LGTO, GTO &RGTO, const double rz[])
{
    return LGTO.coefficient * RGTO.coefficient * ZIntegral(LGTO.cartesian, rz, RGTO.cartesian, LGTO.ang, RGTO.ang, LGTO.orbital_exponent, RGTO.orbital_exponent, 0);
}

double Orbital_ZIntegral(Orbital &LOrbital, Orbital &ROrbital, const double rz[])
{
    double result = 0;
    for (auto &gto1 : LOrbital.component)
        for (auto &gto2 : ROrbital.component)
            result += GTO_ZIntegral(gto1, gto2, rz);
    return result;
}

double JIntegral(const double ra[], const double rb[], const int anga[], const int angb[], double alpha, double beta, int m)
{
    // if any of anga or angb low than 0 return 0;
    for (int i = 0; i < 3; i++)
        if (anga[i] < 0 || angb[i] < 0)
            return 0;
    double zeta = alpha + beta;
    double xi = alpha * beta / zeta;
    // AB squared
    double AB = 0;
    double P[3];
    for (int i = 0; i < 3; i++)
    {
        AB += pow(ra[i] - rb[i], 2);
        P[i] = (alpha * ra[i] + beta * rb[i]) / zeta;
    }

    for (int i = 0; i < 3; i++)
    {
        if (anga[i] > 0)
        {
            // a minus 1 in one direction
            int am[3]{anga[0], anga[1], anga[2]};
            // a minus 2 in one direction
            int amm[3]{anga[0], anga[1], anga[2]};
            // b minus 1 in one direction
            int bm[3]{angb[0], angb[1], angb[2]};
            am[i] -= 1;
            amm[i] -= 2;
            bm[i] -= 1;
            return (P[i] - ra[i]) * JIntegral(ra, rb, am, angb, alpha, beta, m + 1) + 0.5 * am[i] / alpha * (JIntegral(ra, rb, amm, angb, alpha, beta, m) - beta / zeta * JIntegral(ra, rb, amm, angb, alpha, beta, m + 1)) + 0.5 * angb[i] / zeta * JIntegral(ra, rb, am, bm, alpha, beta, m + 1);
        }
        if (angb[i] > 0)
        {
            // a minus 1 in one direction
            int am[3]{anga[0], anga[1], anga[2]};
            // b minus 2 in one direction
            int bmm[3]{angb[0], angb[1], angb[2]};
            // b minus 1 in one direction
            int bm[3]{angb[0], angb[1], angb[2]};
            am[i] -= 1;
            bm[i] -= 1;
            bmm[i] -= 2;

            return (P[i] - rb[i]) * JIntegral(ra, rb, anga, bm, alpha, beta, m + 1) + 0.5 * bm[i] / beta * (JIntegral(ra, rb, anga, bmm, alpha, beta, m) - alpha / zeta * JIntegral(ra, rb, anga, bmm, alpha, beta, m + 1)) + 0.5 * anga[i] / zeta * JIntegral(ra, rb, am, bm, alpha, beta, m + 1);
        }
    }
    //all zero
    return 2 * pow(M_PI, 2.5) / alpha / beta / sqrt(zeta) * Boys_Function(xi * AB, m);
}

double GTO_JIntegral(GTO &LGTO, GTO &RGTO)
{
    return LGTO.coefficient * RGTO.coefficient * JIntegral(LGTO.cartesian, RGTO.cartesian, LGTO.ang, RGTO.ang, LGTO.orbital_exponent, RGTO.orbital_exponent, 0);
}

double Orbital_JIntegral(Orbital &LOb, Orbital &ROb)
{
    double result = 0;
    for (auto &LGTO : LOb.component)
        for (auto &RGTO : ROb.component)
            result += GTO_JIntegral(LGTO, RGTO);
    return result;
}

std::vector<GTO> *Two_Electron_Transform(const Orbital &Oba, const Orbital &Obb)
{
    std::vector<GTO> *gtos = new std::vector<GTO>;

    double AB = 0;
    double zeta, xi, alpha, beta, P[3], coeff;
    int ang[3];
    for (int i = 0; i < 3; i++)
    {
        AB += pow(Oba.cartesian[i] - Obb.cartesian[i], 2);
    }

    for (auto &gtoa : Oba.component)
    {
        for (auto &gtob : Obb.component)
        {
            alpha = gtoa.orbital_exponent;
            beta = gtob.orbital_exponent;
            zeta = alpha + beta;
            xi = alpha * beta / zeta;
            for (int i = 0; i < 3; i++)
            {
                P[i] = (alpha * gtoa.cartesian[i] + beta * gtob.cartesian[i]) / zeta;
            }

            for (ang[0] = 0; ang[0] <= gtoa.ang[0] + gtob.ang[0]; ang[0]++)
                for (ang[1] = 0; ang[1] <= gtoa.ang[1] + gtob.ang[1]; ang[1]++)
                    for (ang[2] = 0; ang[2] <= gtoa.ang[2] + gtob.ang[2]; ang[2]++)
                    {
                        GTO gto;
                        // calculate coefficient
                        // first , coeff of gtoa and gtob themselves
                        // while reading basis gto is born with normalized coefficient
                        coeff = gtoa.coefficient * gtob.coefficient;
                        // second calculate exp coeff
                        coeff *= exp(-xi * AB);
                        // third calculate accumalated bionomal expasion coefficient
                        for (int i = 0; i < 3; i++)
                        {
                            coeff *= Expand_coeff(ang[i], gtoa.ang[i], gtob.ang[i], P[i] - gtoa.cartesian[i], P[i] - gtob.cartesian[i]);
                        }

                        // if coeff equal to zero then ignore it
                        if (fabs(coeff) > COEFF_IGNORE_LIMIT)
                        { // set exponent , ang , cartesian and coeff
                            gto.orbital_exponent = zeta;
                            for (int i = 0; i < 3; i++)
                            {
                                gto.ang[i] = ang[i];
                                gto.cartesian[i] = P[i];
                            }
                            gto.coefficient = coeff;

                            // put it into gtos to return
                            gtos->push_back(gto);
                        }
                    }
        }
    }

    return gtos;
}

double Two_Electron_JIntegral(const Orbital &l1, const Orbital &l2, const Orbital &r1, const Orbital &r2)
{
    double result = 0;
    std::vector<GTO> *left, *right;
    left = Two_Electron_Transform(l1, l2);
    right = Two_Electron_Transform(r1, r2);

    for (auto &LGTO : *left)
        for (auto &RGTO : *right)
        {
            result += GTO_JIntegral(LGTO, RGTO);
        }
    delete left;
    delete right;

    return result;
}
double Expand_coeff(int m, int anga, int angb, double PA, double PB)
{
    double result = 0;
    for (int i = 0; i <= m; i++)
    {
        if (anga < i || angb - m + i < 0)
            continue;
        result += Binomials(anga, i) * Binomials(angb, m - i) * pow(PA, anga - i) * pow(PB, angb - m + i);
    }
    return result;
}

double Binomials(unsigned int n, unsigned int k)
{
    double result = 1;
    if (n < k)
        throw "k should lower than n";
    if (n == k || k == 0)
        return 1;
    for (int i = n; i > n - k; i--)
    {
        result *= i;
    }
    for (int i = 1; i <= k; i++)
    {
        result /= i;
    }
    return result;
}