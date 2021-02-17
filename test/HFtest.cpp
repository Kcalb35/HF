#include <iostream>
#include <gtest/gtest.h>
#include <cmath>
#include "../include/Intergal.h"
#include "../include/basis.h"

#define MAX_TEST_ERR 1e-5

using namespace std;

void Init_GTO(GTO &gto)
{
    InitGTOAng(gto, 0, 0, 0);
    InitGTOCartisian(gto, 0, 0, 0);
    gto.orbital_exponent = 3.2;
}

void NormalizeGTO(GTO &gto)
{
    gto.coefficient = Normalize(gto.orbital_exponent, gto.ang[0], gto.ang[1], gto.ang[2]);
}

TEST(test_SIntergal, gto)
{
    GTO gto;
    Init_GTO(gto);
    NormalizeGTO(gto);
    EXPECT_LE(fabs(GTO_SIntergal(gto, gto) - 1), MAX_TEST_ERR) << "SIntergal for same gto is not 1." << gto.ang[0] << gto.ang[1] << gto.ang[2];

    for (size_t i = 0; i < 3; i++)
    {
        Init_GTO(gto);
        gto.ang[i] = 2;
        gto.cartesian[i] = 3.5;
        NormalizeGTO(gto);
        EXPECT_LE(fabs(GTO_SIntergal(gto, gto) - 1), MAX_TEST_ERR) << "SIntergal for same gto is not 1." << gto.ang[0] << gto.ang[1] << gto.ang[2];
    }
}

TEST(test_KineticIntergal, gto)
{
    GTO gto;
    InitGTOAng(gto, 0, 1, 2);
    InitGTOCartisian(gto,3, 4, 5);
    gto.orbital_exponent = 3.1;
    NormalizeGTO(gto);

    EXPECT_LE(fabs(GTO_Kinetic_Intergal(gto, gto) - 9.81667), MAX_TEST_ERR);
}