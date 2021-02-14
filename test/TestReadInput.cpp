#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "../include/basis.h"
#include "../include/input.h"

int main(int argc, char const *argv[])
{
    using namespace std;
    std::vector<Atom> basis;
    std::vector<Atom> atoms;
    HFoption option;
    std::fstream InputFile;
    std::vector<Orbital> all_orbitals;

    cout << "start read file" <<endl;
    InputFile.open("input_example/H2O.in",std::ios::in);
    ReadInputFile(InputFile,option,atoms,basis,all_orbitals);

    cout << "read file done" <<endl;
    // test option
    cout <<option.BASIS_NAME <<endl;
    cout <<"number of electrons:" <<option.ELECTRON_NUMBER <<endl;
    cout <<"numebr of orbitals:" << option.ORBITAL_NUMBER <<endl;

    // test basis
    for (auto &atom : atoms){
        cout << atom.name << " ";
        for (size_t i = 0; i < 3; i++)
        {
            cout <<atom.cartesian[i] <<" ";
        }
        cout<<endl;
        for(auto &ob:atom.Orbitals){
           cout <<ob.name<<" ";
           for (size_t i = 0; i < 3; i++)
           {
               cout << ob.cartesian[i]<<" ";
           }
           cout <<endl;
            
        }        
    }
}
