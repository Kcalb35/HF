_Pragma("once")

#include <fstream>
#include <vector>
#include <string>

#include "../include/basis.h"

struct HFoption
{
    // basis name
    std::string BASIS_NAME;

    // max energy difference error
    double MAX_ERR = 1e-6;
    // maximum SCF iteration
    int SCF_MAX = 50;
    // repeat times to meet converage
    int counter = 3;
    // whether to print initial guess
    bool SCF_INITIAL_PRINT = false;
    // whether to print fock matrix during SCF
    bool SCF_FOCK_PRINT = false;
    // whether to print coefficient matrix during SCF
    bool SCF_COEF_PRINT = false;

    // number of electrons
    int ELECTRON_NUMBER = 0;
    // number of all orbitals
    int ORBITAL_NUMBER = 0;
};

void ReadInputFile(std::fstream &ifile, HFoption &option, std::vector<Atom> &atom_list, std::vector<Atom> &AtomListByOrder,std::vector<Orbital> & AllOrbitals);