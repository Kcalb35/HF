#include <iostream>
#include <fstream>
#include "../include/basis.h"

int main(int argc, char const *argv[])
{
    using namespace std;

    fstream ifile;
    ifile.open("basis/STO-3G.txt",ios::in);

    auto atoms = ReadBasis(ifile);

    Atom Li = atoms[1];
    cout << Li.name <<endl;
    for(auto &orbital : Li.Orbitals){
        cout << orbital.name <<" " << orbital.n << orbital.l<<orbital.m<<endl;
        for(auto gto : orbital.component){
            cout <<gto.orbital_exponent << " " << gto.coefficient <<endl;
            cout << gto.ang[0]<<" "<<gto.ang[1]<< " "<<gto.ang[2]<<endl;
        }
    }

    ifile.close();

    
    return 0;
}
