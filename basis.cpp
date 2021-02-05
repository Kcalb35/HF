#include <vector>
#include <string>

struct GTO
{
    int ang[3];
    double cartesian[3];
    double orbital_exponent;
    double coefficient;
};

struct Orbital
{
    // quantum number
    int n;
    // angular number
    int l;
    // magentic number
    int m;

    std::string name;

    // an orbital consist of GTO s
    std::vector<GTO> component;

    double cartesian[3];
};

struct Atom
{
    // atom number
    int n;
    // atom name
    std::string name;
    // an atom consists of orbitals
    std::vector<Orbital> Orbitals;
    // cartesian position
    double cartesian[3];
};

// TODO 双电子积分

std::vector<Atom> basis_scanf(FILE *basis)
{
    // 先实现STO-3G
    std::vector<Atom> atoms;


    return atoms;
}