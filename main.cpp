#include <cstring>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <iomanip>

#include "include/basis.h"
#include "include/input.h"
#include "include/RHF.h"



int main(int argc, char const *argv[])
{
    using namespace std;
    string input;

    HFoption option;
    fstream ifstream;
    vector<Atom> atomByOrder,allAtoms;
    vector<Orbital> allOrbitals;

    gsl_vector * energy;
    gsl_matrix * coeff;
    double total_energy;

    // reading command line params
    for (int i = 0; i < argc; ++i) {
       if (strcmp(argv[i],"-f")==0){
          input.append(argv[i+1]);
       }
       else if(strcmp(argv[i],"-h")==0){
          cout <<  "-h :print help message\n";
          cout << "-f <filename> : specify the input file\n";
          return 0;
       }
    }

    // reading file into option
    ifstream.open(input,ios_base::in);
    if (!ifstream.is_open()){
        cout << "ERROR: NO SUCH INPUT FILE\n";
        return 1;
    }
    ReadInputFile(ifstream,option,allAtoms,atomByOrder,allOrbitals);


    // initialize
    energy = gsl_vector_calloc(option.ORBITAL_NUMBER);
    coeff = gsl_matrix_calloc(option.ORBITAL_NUMBER,option.ORBITAL_NUMBER);

    // start RHF
    RHF_SCF(&total_energy,energy,coeff,allAtoms,allOrbitals,option);

    // print atom orbital labels
    cout <<"Atom orbital labels : [";
    bool first=true;
    for(auto &ob:allOrbitals){
        if(first) {
            cout << ob.name;
            first = false;
        }
        else cout <<", "<<ob.name;
    }
    cout <<"]\n";

    // print energy
    cout <<"Total energy : "<<setprecision(12)<<total_energy + nuclei_repulsion(allAtoms)<<endl;


    cout.precision(6);
    // print MOs
    cout <<"MOs"<<endl;
    for (int i = 0; i < option.ORBITAL_NUMBER; ++i) {
        cout << i<<" energy="<<setw(10)<<gsl_vector_get(energy,i)<<", occ=";
        if (2*i<option.ELECTRON_NUMBER) cout<<2<<endl;
        else cout <<0<<endl;
    }
    return 0;
}
