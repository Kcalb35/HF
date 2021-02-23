//
// Created by grasscube on 2021/2/21.
//

#include "../include/RHF.h"
#include "../include/gslextra.h"
#include "../include/basis.h"
#include "../include/Integral.h"
#include "../include/input.h"
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <cmath>
#include <iostream>
#include <iomanip>

/// single electron hamiltonian element
/// \param LOb left orbital
/// \param ROb right orbital
/// \param atoms list of atoms to calculate nuclear attraction integral
/// \return integral value
double single_electron_hamiltonian_element(Orbital &LOb, Orbital &ROb, std::vector<Atom> &atoms){
    double result = 0;
    // calculate kinetic integral
    result += Orbital_Kinetic_Integral(LOb,ROb);
    result += nuclear_attraction_energy(LOb,ROb,atoms);
    return result;
}

/// set core hamiltonian matrix to dest matrix
/// \param dest destiny matrix, should be pre-calloced
/// \param obs a list of all orbitals
/// \param atoms a list of all atoms
void core_hamiltonian_matrix_set(gsl_matrix * dest, std::vector<Orbital> &obs, std::vector<Atom> &atoms){
    // todo optimize with symmetric matrix
    int length = obs.size();
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            gsl_matrix_set(dest,i,j,single_electron_hamiltonian_element(obs[i],obs[j],atoms));
        }
    }
}

/// make initial guess
/// \param coeff_matrix
/// \param core_hamiltonian_matrix
/// \param overlap_matrix
void initial_guess(gsl_matrix * coeff_matrix,gsl_matrix * core_hamiltonian_matrix,gsl_matrix * overlap_matrix){

    gsl_vector * eigen = gsl_vector_calloc(core_hamiltonian_matrix->size1);

    gsl_eigen_Lowdin_diag(core_hamiltonian_matrix,overlap_matrix,eigen,coeff_matrix);

    gsl_vector_free(eigen);
}

/// calculate density matrix
/// \param density_matrix destiny density matrix
/// \param coeff_matrix
/// \param electron_number
void density_matrix_set(gsl_matrix * density_matrix,gsl_matrix * coeff_matrix,int electron_number){
    int length = coeff_matrix->size1;
    gsl_vector * v = gsl_vector_calloc(length);
    gsl_matrix * m1=  gsl_matrix_calloc(length,length);

    for (int i = 0; i < electron_number/2; ++i) {
        gsl_matrix_get_col(v,coeff_matrix,i);
        gsl_matrix_set_col(m1,i,v);
    }

    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,m1,m1,0,density_matrix);

    gsl_matrix_free(m1);
    gsl_vector_free(v);
}

/// fock matrix element calculation
/// \param tensor
/// \param density_matrix
/// \param h_core_matrix
/// \param i
/// \param j
/// \return
double Fock_matrix_element(gsl_quad_tensor * tensor,gsl_matrix * density_matrix,gsl_matrix * h_core_matrix,int i,int j){
    double result =0;
    int length = density_matrix->size1;
    result += gsl_matrix_get(h_core_matrix,i,j);
    for (int k = 0; k < length; ++k) {
        for (int l = 0; l < length; ++l) {
            // coulomb integral
           result += 2 * gsl_matrix_get(density_matrix,k,l) * gsl_quad_tensor_get(tensor,i,j,k,l);
           result -= gsl_matrix_get(density_matrix,k,l) * gsl_quad_tensor_get(tensor,i,l,k,j);
        } 
    }
    return result;
}


/// fock matrix set
/// \param fock result fock matrix
/// \param v
/// \param density
/// \param h_core
void Fock_matrix(gsl_matrix * fock,gsl_quad_tensor * v,gsl_matrix * density,gsl_matrix * h_core){
    int length = density->size1;
    for (int i = 0; i <length ; ++i) {
        for (int j = 0; j < length; ++j) {
            gsl_matrix_set(fock,i,j,Fock_matrix_element(v,density,h_core,i,j));
        }
    }
}

/// calculate nuclei repulsion energy Vn
/// \param atoms a list of all atoms
/// \return Vnn
double nuclei_repulsion(const std::vector<Atom> &atoms){
    if (atoms.size() <=1){
        return 0;
    }
    double result = 0;
    double distance;
    for (auto &a:atoms){
        for(auto &b:atoms){
            if (&a !=&b){
                distance = sqrt(pow(a.cartesian[0]-b.cartesian[0],2)+pow(a.cartesian[1]-b.cartesian[1],2)+pow(a.cartesian[2]-b.cartesian[2],2));
                result += a.n * b.n /distance ;
            }
        }
    }
    return result /2.0;
}

/// HF energy without Vnn
/// \param v quad tensor
/// \param density density matrix
/// \param h_core h core matrix
/// \return HF energy without Vnn
double HF_energy(gsl_quad_tensor *v,gsl_matrix * density,gsl_matrix * h_core){
    double result = 0;
    double pij;
    int length = density->size1;
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            pij = gsl_matrix_get(density,i,j);
            result += gsl_matrix_get(h_core,i,j) *pij ;
            for (int k = 0; k < length; ++k) {
                for (int l = 0; l < length; ++l) {
                    result += pij * gsl_matrix_get(density,k,l) * (gsl_quad_tensor_get(v,i,j,k,l)-0.5 *gsl_quad_tensor_get(v,i,l,k,j));
                }
            }
        }
    }
    return 2 * result;
}

void two_electron_quad_tensor_set(gsl_quad_tensor * dest,std::vector<Orbital> &obs){
    int len = obs.size();
    // todo optimize
    for (int i = 0; i <len; ++i) {
        for (int j = 0; j < len; ++j) {
            for (int k = 0; k < len; ++k) {
                for (int l = 0; l < len; ++l) {
                   gsl_quad_tensor_set(dest,i,j,k,l,Two_Electron_JIntegral(obs[i],obs[j],obs[k],obs[l]));
                }
            }
        }
    }
}


/// set kinetic integral matrix
/// \param m output matrix ,which should be pre-allocated
/// \param obs orbitals
void kinetic_matrix_set(gsl_matrix *m, std::vector<Orbital> &obs) {
    int len = obs.size();
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            gsl_matrix_set(m,i,j,Orbital_Kinetic_Integral(obs[i],obs[j]));
        }
    }
}

double nuclear_attraction_energy(Orbital &LOb, Orbital &ROb, std::vector<Atom> &atoms) {
   double result = 0;
   for(auto &atom:atoms){
       result -= atom.n * Orbital_ZIntegral(LOb,ROb,atom.cartesian);
   }
    return result;
}

/// set nuclear attraction energy matrix,
/// \param m output matrix should be pre-allocated
/// \param obs all orbitals
/// \param atoms all atoms
void nuclear_attraction_energy_matrix_set(gsl_matrix *m, std::vector<Orbital> obs,std::vector<Atom> &atoms) {
   int len = obs.size();
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            gsl_matrix_set(m, i, j, nuclear_attraction_energy(obs[i], obs[j], atoms));
        }
    }
}


/// RHF SCF method
/// \param total_energy output total energyLevel
/// \param energyLevel all energyLevel levels
/// \param coeff coeff matrix
/// \param allAtoms all atoms in the system
/// \param allOrbitals all orbitals of atoms
/// \param option HF option
void RHF_SCF(double * total_energy, gsl_vector * energyLevel, gsl_matrix * coeff, std::vector<Atom> &allAtoms, std::vector<Orbital> &allOrbitals, const HFoption &option){
    int ORBITAL_NUMBER = option.ORBITAL_NUMBER;
    double fock_energy,tmp_energy;

    // pre allocate
    gsl_matrix * H_core = gsl_matrix_calloc(ORBITAL_NUMBER,ORBITAL_NUMBER);
    gsl_matrix * S_overlap = gsl_matrix_calloc(ORBITAL_NUMBER,ORBITAL_NUMBER);
    gsl_quad_tensor * quadTensor = gsl_quad_tensor_calloc(ORBITAL_NUMBER,ORBITAL_NUMBER,ORBITAL_NUMBER,ORBITAL_NUMBER);
    gsl_matrix * P_density = gsl_matrix_calloc(ORBITAL_NUMBER,ORBITAL_NUMBER);
    gsl_matrix * Fock = gsl_matrix_calloc(ORBITAL_NUMBER,ORBITAL_NUMBER);
    gsl_matrix * debug_info = gsl_matrix_calloc(ORBITAL_NUMBER,ORBITAL_NUMBER);


    // preparing data, calculate initial guess
    using namespace std;
    cout << "Electron numbers:"<<option.ELECTRON_NUMBER<<endl;
    cout <<"Nuclear repulsion:"<< nuclei_repulsion(allAtoms)<<endl;

    Orbital_S_matrix_set(S_overlap,allOrbitals);
    core_hamiltonian_matrix_set(H_core,allOrbitals,allAtoms);
    two_electron_quad_tensor_set(quadTensor,allOrbitals);

    initial_guess(coeff,H_core,S_overlap);
    density_matrix_set(P_density,coeff,option.ELECTRON_NUMBER);

    // output initial guess
    if (option.SCF_INITIAL_PRINT){
        cout << endl<<"Overlap Integrals:"<<endl;
        gsl_matrix_print(S_overlap);

        cout <<endl<<"Kinetic Integrals:" <<endl;
        kinetic_matrix_set(debug_info,allOrbitals);
        gsl_matrix_print(debug_info);

        cout <<endl<<"Nuclear attraction Integrals:"<<endl;
        nuclear_attraction_energy_matrix_set(debug_info,allOrbitals,allAtoms);
        gsl_matrix_print(debug_info);

        cout << endl<<"Core Hamiltonian Matrix:"<<endl;
        gsl_matrix_print(H_core);

        cout <<endl<<"Coefficient matrix:"<<endl;
        gsl_matrix_print(coeff);

        cout<<endl<<"Initial density matrix:"<<endl;
        gsl_matrix_print(P_density);
    }

    Fock_matrix(Fock,quadTensor,P_density,H_core);
    fock_energy = HF_energy(quadTensor,P_density,H_core);

    // print initial fock matrix and HF energyLevel
    cout <<endl<<"Fock Matrix:"<<endl;
    gsl_matrix_print(Fock);
    cout <<endl<<"HF energyLevel:"<<fock_energy<<endl;

    // start SCF
    cout <<"=============== Start HF ==============="<<endl;
    bool IsConverged = false;
    int ConvergeCount = 0;
    for (int i = 0; i < option.SCF_MAX; ++i) {
        Fock_matrix(Fock,quadTensor,P_density,H_core);
        gsl_eigen_Lowdin_diag(Fock, S_overlap, energyLevel,coeff);
        density_matrix_set(P_density,coeff,option.ELECTRON_NUMBER);

        tmp_energy = HF_energy(quadTensor,P_density,H_core);

        if (fabs(tmp_energy-fock_energy)<option.MAX_ERR) ConvergeCount ++;
        else ConvergeCount=0;

        // reach converged
        if (ConvergeCount >= option.counter){
            IsConverged = true;
            break;
        }
        // not converged print message
        cout <<"iteration = "<<i+1<<", HF energy = "<<tmp_energy <<endl;
        if (option.SCF_FOCK_PRINT){
            cout <<endl<<"Fock Matrix:"<<endl;
            gsl_matrix_print(Fock);
        }
        if (option.SCF_COEF_PRINT){
            cout <<endl<<"Coefficient matrix:"<<endl;
            gsl_matrix_print(coeff);
        }

        // update energy
        fock_energy = tmp_energy;
    }

    // end SCF
    cout <<"=============== End HF ==============="<<endl;
    if (!IsConverged){
        cout <<"WARNING: SCF not converged."<<endl;
    }
    else{
        cout <<"INFO: SCF converged."<<endl;
    }

    *total_energy = fock_energy;
    gsl_matrix_free(H_core);
    gsl_matrix_free(S_overlap);
    gsl_matrix_free(P_density);
    gsl_matrix_free(Fock);
    gsl_matrix_free(debug_info);
    gsl_quad_tensor_free(quadTensor);
}

/// calculate S overlap matrix
/// \param S output matrix, S should be pre-allocated
/// \param orbitals all orbitals
void Orbital_S_matrix_set(gsl_matrix * S,std::vector<Orbital>& orbitals){
    int len = orbitals.size();
    // todo optimize symmetric
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            gsl_matrix_set(S,i,j,Orbital_SIntegral(orbitals[i],orbitals[j]));
        }
    }
}

