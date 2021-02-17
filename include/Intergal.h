#ifndef INTERGAL_H
#define INTERGAL_H
#include "../include/basis.h"

double SIntergal(double ra[], double rb[], int ax, int ay, int az, int bx, int by, int bz, double alpha, double beta);

double Orbital_SIntergal(Orbital &ob1, Orbital &ob2);

double GTO_SIntergal(GTO &gto1, GTO &gto2);

std::vector<GTO> GTO_Derivative(GTO &gto,int dir);

std::vector<GTO> GTO_SecondOrder_Derivative(GTO &gto,int dir);

double GTO_Kinetic_Intergal(GTO &LGTO, GTO &RGTO);

double Orbital_Kinetic_Intergal(Orbital &LOrbital, Orbital &ROrbital);
#endif