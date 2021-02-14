#include <cmath>

#include "../include/basis.h"
#include "../include/Intergal.h"

double SIntergal(double ra[], double rb[], int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta)
{

    if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0)
    {
        return 0;
    }
    double zeta = alpha + beta;
    if (ax > 0)
    {
        double pi = (alpha * ra[0] + beta * rb[0]) / zeta;
        return (pi - ra[0]) * SIntergal(ra, rb, ax - 1, ay, az, bx, by, bz, alpha, beta) + 0.5 * bx / zeta * SIntergal(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha, beta) + 0.5 * (ax - 1) / alpha * SIntergal(ra, rb, ax - 2, ay, az, bx, by, bz, alpha, beta);
    }
    else if (ay > 0)
    {
        double pi = (alpha * ra[1] + beta * rb[1]) / zeta;
        return (pi - ra[1]) * SIntergal(ra, rb, ax, ay - 1, az, bx, by, bz, alpha, beta) + 0.5 * by / zeta * SIntergal(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha, beta) + 0.5 * (ay - 1) / alpha * SIntergal(ra, rb, ax, ay - 2, az, bx, by, bz, alpha, beta);
    }
    else if (az > 0)
    {
        double pi = (alpha * ra[2] + beta * rb[2]) / zeta;
        return (pi - ra[2]) * SIntergal(ra, rb, ax, ay, az - 1, bx, by, bz, alpha, beta) + 0.5 * bz / zeta * SIntergal(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha, beta) + 0.5 * (az - 1) / alpha * SIntergal(ra, rb, ax, ay, az - 2, bx, by, bz, alpha, beta);
    }
    else if (bx > 0)
    {
        double pi = (alpha * ra[0] + beta * rb[0]) / zeta;
        return (pi - rb[0]) * SIntergal(ra, rb, ax, ay, az, bx - 1, by, bz, alpha, beta) + 0.5 * ax / zeta * SIntergal(ra, rb, ax - 1, ay, az, bx - 1, by, bz, alpha, beta) + 0.5 * (bx - 1) / beta * SIntergal(ra, rb, ax, ay, az, bx - 2, by, bz, alpha, beta);
    }
    else if (by > 0)
    {
        double pi = (alpha * ra[1] + beta * rb[1]) / zeta;
        return (pi - rb[1]) * SIntergal(ra, rb, ax, ay, az, bx, by - 1, bz, alpha, beta) + 0.5 * ay / zeta * SIntergal(ra, rb, ax, ay - 1, az, bx, by - 1, bz, alpha, beta) + 0.5 * (by - 1) / beta * SIntergal(ra, rb, ax, ay, az, bx, by - 2, bz, alpha, beta);
    }
    else if (bz > 0)
    {
        double pi = (alpha * ra[2] + beta * rb[2]) / zeta;
        return (pi - rb[2]) * SIntergal(ra, rb, ax, ay, az, bx, by, bz - 1, alpha, beta) + 0.5 * az / zeta * SIntergal(ra, rb, ax, ay, az - 1, bx, by, bz - 1, alpha, beta) + 0.5 * (bz - 1) / beta * SIntergal(ra, rb, ax, ay, az, bx, by, bz - 2, alpha, beta);
    }
    else
    {
        // all zero
        double eta = alpha * beta / zeta;
        return pow(M_PI / zeta, 1.5) * exp(-eta * (pow(ra[0] - rb[0], 2) + pow(ra[1] - rb[1], 2) + pow(ra[2] - rb[2], 2)));
    }
}

double Orbital_SIntergal(Orbital &ob1, Orbital &ob2)
{
    double result=0;
    for(GTO &gtoa:ob1.component){
        for (GTO &gtob:ob2.component){
            result += GTO_SIntergal(gtoa,gtob);
        }
    }
    return result;
}

double GTO_SIntergal(GTO &gto1, GTO &gto2)
{
    return gto1.coefficient * gto2.coefficient * SIntergal(gto1.cartesian, gto2.cartesian, gto1.ang[0], gto1.ang[1], gto1.ang[2], gto2.ang[0], gto2.ang[1], gto2.ang[2], gto1.orbital_exponent, gto2.orbital_exponent);
}