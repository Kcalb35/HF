_Pragma("once")

#include "../include/basis.h"

double SIntergal(double ra[], double rb[], int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta);

double Orbital_SIntergal(Orbital &ob1, Orbital &ob2);

double GTO_SIntergal(GTO &gto1, GTO &gto2);

std::vector<GTO> GTO_Derivative(GTO &gto,int dir);

std::vector<GTO> GTO_SecondOrder_Derivative(GTO &gto,int dir);

double GTO_Kinetic_Intergal(GTO &LGTO, GTO &RGTO);

double Orbital_Kinetic_Intergal(Orbital &LOrbital, Orbital &ROrbital);
double Boys_Function(double x, int m);
double ZIntergal(const double ra[], const double rc[], const double rb[], const int anga[], const int angb[], double alpha, double beta, int m);
double GTO_ZIntergal(GTO &LGTO, GTO &RGTO, const double rz[]);
double Orbital_ZIntergal(Orbital &LOrbital, Orbital &ROrbital,const double rz[]);