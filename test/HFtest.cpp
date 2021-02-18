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
    InitGTOCartisian(gto, 3, 4, 5);
    gto.orbital_exponent = 3.1;
    NormalizeGTO(gto);

    EXPECT_LE(fabs(GTO_Kinetic_Intergal(gto, gto) - 9.81667), MAX_TEST_ERR);
}

TEST(test_ZIntergal, boysfunction)
{
    EXPECT_LE(fabs(Boys_Function(1.5, 2) - 0.0723635), MAX_TEST_ERR);
    EXPECT_LE(fabs(Boys_Function(0.5, 2) -0.140751), MAX_TEST_ERR);
    EXPECT_LE(fabs(Boys_Function(0.005, 2) -0.199287), MAX_TEST_ERR);
}

TEST(test_ZIntergal, gto)
{
    double ra[]{1,1,1};
    double rz[3]{0, 0, 0};
    int ang[]{1,2,3};

    double err = ZIntergal(ra,rz,ra,ang,ang,3.1,3.2,0);
    EXPECT_LE(fabs(err - 2.37528e-6), MAX_TEST_ERR);
}