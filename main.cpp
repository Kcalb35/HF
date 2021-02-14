#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
double double_factorial(int n)
{
    if (n == 0 || n == 1)
        return 1;
    else if (n == -1)
        return 1;
    else
        return n * double_factorial(n - 2);
}

//Calculate the normalization coefficient for the gaussian orbital
double normalize(double alpha, int ax, int ay, int az)
{
    return pow(2 * alpha / M_PI, 0.75) * pow(4 * alpha, (double)(ax + ay + az) / 2.0) / sqrt(double_factorial(2 * ax - 1) * double_factorial(2 * ay - 1) * double_factorial(2 * az - 1));
}
double Normalize(double exponent, int i, int j, int k)
{
    return pow(2 * exponent / M_PI, 0.75) * pow((pow(8 * exponent, i + j + k) * std::tgamma(i + 1) * std::tgamma(j + 1) * std::tgamma(k + 1)) / (std::tgamma(2 * i + 1) * std::tgamma(2 * j + 1) * std::tgamma(2 * k + 1)), 0.5);
}

int main(int argc, char const *argv[])
{
    cout << normalize(3.4252, 0, 0, 0) << endl;
    cout << Normalize(3.4252, 0, 0, 0) <<endl;
}
