_Pragma("once")

#include "../include/basis.h"

double SIntegral(double *ra, double *rb, int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta);

double Orbital_SIntegral(Orbital &ob1, Orbital &ob2);

double GTO_SIntegral(GTO &gto1, GTO &gto2);

std::vector<GTO> GTO_Derivative(GTO &gto,int dir);

std::vector<GTO> GTO_SecondOrder_Derivative(GTO &gto,int dir);

double GTO_Kinetic_Integral(GTO &LGTO, GTO &RGTO);

double Orbital_Kinetic_Integral(Orbital &LOrbital, Orbital &ROrbital);
double Boys_Function(double x, int m);
double ZIntegral(const double *ra, const double *rc, const double *rb, const int *anga, const int *angb, double alpha, double beta, int m);
double GTO_ZIntegral(GTO &LGTO, GTO &RGTO, const double *rz);
double Orbital_ZIntegral(Orbital &LOrbital, Orbital &ROrbital, const double *rz);

double JIntegral(const double ra[], const double rb[], const int anga[], const int angb[], double alpha, double beta, int m);
double GTO_JIntegral(GTO &LGTO, GTO &RGTO);
double Orbital_JIntegral(Orbital &LOb, Orbital &ROb);
double Two_Electron_JIntegral(const Orbital &l1, const Orbital &l2, const Orbital &r1, const Orbital &r2);

std::vector<GTO> *Two_Electron_Transform(const Orbital &Oba, const Orbital &Obb);
double Expand_coeff(int m, int anga, int angb, double PA, double PB);
double Binomials(unsigned int n, unsigned int k);
