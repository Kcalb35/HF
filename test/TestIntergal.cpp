#include <iostream>
#include <fstream>

#include "../include/basis.h"
#include "../include/Intergal.h"

using namespace std;

void PrintOb(Orbital &ob)
{
    cout << ob.name<< " " <<ob.n<<ob.l<<ob.m <<endl;
    // for(GTO &gto : ob.component){
    //     cout << "coeef "<<  gto.coefficient <<endl;
    //     cout << "exponent " <<gto.orbital_exponent<<endl;
    // }
}

int main(int argc, char const *argv[])
{
    fstream fs;
    fs.open("basis/STO-3G.txt");

    auto atoms = ReadBasis(fs);
    Atom H = atoms[0];
    Orbital ob = H.Orbitals[0];
    cout <<"--------"<<endl;
    PrintOb(ob);
    cout << "SIntergal: "<<Orbital_SIntergal(ob,ob)<<endl;


    cout <<"--------"<<endl;
    Atom Li = atoms[2];
    for(auto &ob:Li.Orbitals){
        PrintOb(ob);
        cout << "SIntergal: "<<Orbital_SIntergal(ob,ob)<<endl;
    }

    cout <<"----gto intergal test----"<<endl;
    GTO gto;
    gto.orbital_exponent = 2;
    gto.ang[0] = gto.ang[1] = gto.ang[2] = 0;
    gto.cartesian[0] = gto.cartesian[1] = gto.cartesian[2] = 0;
    gto.coefficient = Normalize(2,0,0,0);
    cout << "GTO intergal:"<< GTO_SIntergal(gto,gto) << endl;
    return 0;
}
