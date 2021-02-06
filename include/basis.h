#include <vector>
#include <string>
#include <fstream>

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

    // obital name
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

std::vector<Atom> ReadBasis(std::fstream &ifile);
Orbital SetOrbital(int m, Orbital ob);