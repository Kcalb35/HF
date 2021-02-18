#include <iostream>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../include/Integral.h"
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
void PrintGTO(GTO &gto)
{
    cout << gto.ang[0] << gto.ang[1] << gto.ang[2] << " " << gto.cartesian[0] << " " << gto.cartesian[1] << " " << gto.cartesian[2] << " " << gto.orbital_exponent << " " << gto.coefficient << endl;
}

TEST(SIntegral, gto)
{
    GTO gto;
    Init_GTO(gto);
    NormalizeGTO(gto);
    EXPECT_LE(fabs(GTO_SIntegral(gto, gto) - 1), MAX_TEST_ERR) << "SIntegral for same gto is not 1." << gto.ang[0] << gto.ang[1] << gto.ang[2];

    for (size_t i = 0; i < 3; i++)
    {
        Init_GTO(gto);
        gto.ang[i] = 2;
        gto.cartesian[i] = 3.5;
        NormalizeGTO(gto);
        EXPECT_LE(fabs(GTO_SIntegral(gto, gto) - 1), MAX_TEST_ERR) << "SIntegral for same gto is not 1." << gto.ang[0] << gto.ang[1] << gto.ang[2];
    }
}

TEST(KineticIntegral, gto)
{
    GTO gto;
    InitGTOAng(gto, 0, 1, 2);
    InitGTOCartisian(gto, 3, 4, 5);
    gto.orbital_exponent = 3.1;
    NormalizeGTO(gto);

    EXPECT_LE(fabs(GTO_Kinetic_Integral(gto, gto) - 9.81667), MAX_TEST_ERR);
}

TEST(ZIntegral, boysfunction)
{
    EXPECT_LE(fabs(Boys_Function(1.5, 2) - 0.0723635), MAX_TEST_ERR);
    EXPECT_LE(fabs(Boys_Function(0.5, 2) - 0.140751), MAX_TEST_ERR);
    EXPECT_LE(fabs(Boys_Function(0.005, 2) - 0.199287), MAX_TEST_ERR);
}

TEST(ZIntegral, gto)
{
    double ra[]{1, 1, 1};
    double rz[3]{0, 0, 0};
    int ang[]{1, 2, 3};

    double err = ZIntegral(ra, rz, ra, ang, ang, 3.1, 3.2, 0);
    EXPECT_LE(fabs(err - 2.37528e-6), MAX_TEST_ERR);
}

TEST(JIntegral, allzero)
{
    double r0[]{0, 0, 0};
    int ang0[]{0, 0, 0};
    double err = JIntegral(r0, r0, ang0, ang0, 3.1, 2.1, 0);
    EXPECT_LE(fabs(err - 2.3568), MAX_TEST_ERR);
}
TEST(JIntegral, gto)
{
    double ra[]{1, 2, 3};
    int ang[]{2, 2, 2};
    double result = JIntegral(ra, ra, ang, ang, 3.1, 2.1, 0);
    EXPECT_LE(fabs(result - 7.64505e-5), MAX_TEST_ERR);
}

TEST(JIntegral, Binomials)
{
    EXPECT_EQ((int)Binomials(6, 2), 15);
    EXPECT_EQ((int)Binomials(10, 4), 210);
    EXPECT_EQ((int)Binomials(5, 5), 1);
    EXPECT_EQ((int)Binomials(5, 0), 1);
}

// {{0.0047396, 0.00378586, -0.0196865}
//  {0.00378586, 0.00302403, -0.015725}
//  {-0.0196865, -0.015725, 0.0817699}}
TEST(Two_Electron_Transform, coeff)
{
    GTO gto1{{1, 1, 1}, {1, 2, 3}, 3.1, 1};
    GTO gto2{{1, 1, 1}, {2, 3, 3}, 2.1, 1};
    Orbital o1, o2;
    o1.component.push_back(gto1);
    o2.component.push_back(gto2);
    o1.cartesian[0] = 1;
    o1.cartesian[1] = 2;
    o1.cartesian[2] = 3;
    o2.cartesian[0] = 2;
    o2.cartesian[1] = 3;
    o2.cartesian[2] = 3;
    std::vector<GTO> *gtos = Two_Electron_Transform(o1, o2);

    for (auto &g : (*gtos))
    {
        PrintGTO(g);
    }

    cout << "--all : " << (*gtos).size() << endl;
    delete gtos;
}