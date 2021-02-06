#include <vector>
#include <string>
#include <fstream>

#include "../include/basis.h"

// TODO 双电子积分

std::vector<Atom> ReadBasis(std::fstream &ifile)
{
    // 先实现STO-3G
    std::vector<Atom> atomsByOrder;
    int NumOfAtoms;

    ifile >> NumOfAtoms;
    for (int i = 0; i < NumOfAtoms; i++)
    {
        // read atom one by one
        int orbitalNum;
        Atom atom;
        ifile >> atom.n >> atom.name >> orbitalNum;
        for (int j = 0; j < orbitalNum; j++)
        {
            // read orbital one by one
            int GTONum;
            Orbital orbital;
            // get n and l of an orbital,also with GTO number
            ifile >> orbital.n >> orbital.l >> GTONum;
            for (int k = 0; k < GTONum; k++)
            {
                GTO gto;
                ifile >> gto.orbital_exponent >> gto.coefficient;
                orbital.component.push_back(gto);
            }
            switch (orbital.l)
            {
            case 0:
                atom.Orbitals.push_back(SetOrbital(0, orbital));
                break;
            case 1:
                // p orbital
                atom.Orbitals.push_back(SetOrbital(0, orbital));
                atom.Orbitals.push_back(SetOrbital(1, orbital));
                atom.Orbitals.push_back(SetOrbital(-1, orbital));
                break;

            case 2:
                // d orbital
                throw "Not Implement";
                break;
            }
        }

        // label and set magentic finish
        atomsByOrder.push_back(atom);
    }
    return atomsByOrder;
}

Orbital SetOrbital(int m, Orbital ob)
{
    ob.m = m;
    switch (ob.l)
    {
    case 0:
        // s orbital
        for (auto &gto : ob.component)
        {
            gto.ang[0] = 0;
            gto.ang[1] = 0;
            gto.ang[2] = 0;
        }
        ob.name = std::to_string(ob.n) + "s";
        break;
    case 1:
        // p orbital
        switch (ob.m)
        {
        case -1:
            // py
            for (auto &gto : ob.component)
            {
                gto.ang[0] = 0;
                gto.ang[1] = 1;
                gto.ang[0] = 0;
            }
            ob.name = std::to_string(ob.n) + "py";
            break;
        case 0:
            // pz
            for (auto &gto : ob.component)
            {
                gto.ang[0] = 0;
                gto.ang[1] = 0;
                gto.ang[2] = 1;
            }
            ob.name = std::to_string(ob.n) + "pz";
            break;
        case 1:
            for (auto &gto : ob.component)
            {
                gto.ang[0] = 1;
                gto.ang[1] = 0;
                gto.ang[2] = 0;
            }
            ob.name = std::to_string(ob.n) + "px";
            break;
        }
        break;
    case 2:
        // d orbital
        throw "Not Implement";
        break;
    case 3:
        // f orbital
        throw "Not Implement";
        break;
    }
    return ob;
}