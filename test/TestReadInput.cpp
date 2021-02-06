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

    InputFile.open("input_example\\H2O.in",std::ios::in);
    ReadInputFile(InputFile,option,atoms,basis);

    // test option
    cout <<option.BASIS_NAME <<endl;
    cout <<"number of electrons:" <<option.ELECTRON_NUMBER <<endl;


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
