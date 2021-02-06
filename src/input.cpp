#include <fstream>
#include <vector>
#include <string>

#include "../include/input.h"

void ReadInputFile(std::fstream &ifile, HFoption &option, std::vector<Atom> &atom_list, std::vector<Atom> &AtomListByOrder)
{
    std::string inputLine;
    std::map<std::string, int> index;
    while (!ifile.eof())
    {
        ifile >> inputLine;
        if (inputLine == "#BASIS")
        {
            ifile >> option.BASIS_NAME;
            std::string basisFilePath = "basis\\" + option.BASIS_NAME + ".txt";

            // start read in basis
            std::fstream basisFstream;
            basisFstream.open(basisFilePath, std::ios::in);
            // TODO file not found situation
            // read basis
            AtomListByOrder = ReadBasis(basisFstream);
            // get index map
            index = GetIndexMap(AtomListByOrder);
        }
        else if (inputLine == "SCF_BEGIN")
        {
            // read SCF options
            while (!ifile.eof())
            {
                ifile >> inputLine;
                if (inputLine == "#ERR")
                    ifile >> option.MAX_ERR;
                else if (inputLine == "#MAX")
                    ifile >> option.SCF_MAX;
                else if (inputLine == "#INITIAL")
                    option.SCF_INITIAL_PRINT = true;
                else if (inputLine == "#FOCK")
                    option.SCF_FOCK_PRINT = true;
                else if (inputLine == "#COEFF")
                    option.SCF_COEF_PRINT = true;
                else if (inputLine == "#SCF_END")
                    break;
            }
        }
        else if (inputLine == "#COORD_BEGIN")
        {
            // read coordinates
            while (!ifile.eof())
            {
                ifile >> inputLine;
                if (inputLine == "#COORD_END")
                    break;
                else
                {
                    int n = index[inputLine];
                    // add up electron numebrs
                    option.ELECTRON_NUMBER += n;
                    Atom atom = AtomListByOrder[n-1];
                    double tmpCOORD;
                    for (int i = 0; i < 3; i++)
                    {
                       ifile >> tmpCOORD;
                       atom.cartesian[i] = tmpCOORD/ 0.529177210903;
                    }
                    // sync coord to every orbital and GTO
                    SyncCoordAndName(atom);
                    atom_list.push_back(atom);                     
                }
            }
        }
    }
}
