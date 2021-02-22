//
// Created by grasscube on 2021/2/21.
//

#include "../include/RHF.h"
#include "../include/gslextra.h"
#include "../include/basis.h"
#include "../include/Integral.h"
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <cmath>

/// single electron hamiltonian element
/// \param LOb left orbital
/// \param ROb right orbital
/// \param atoms list of atoms to calculate nuclear attraction integral
/// \return integral value
double single_electron_hamiltonian_element(Orbital &LOb,Orbital &ROb, const std::vector<Atom> &atoms){
    double result = 0;
    // calculate kinetic integral
    result += Orbital_Kinetic_Integral(LOb,ROb);
    for (auto &atom:atoms){
       result -= atom.n * Orbital_ZIntegral(LOb,ROb,atom.cartesian);
    }
    return result;
}

/// set core hamiltonian matrix to dest matrix
/// \param dest destiny matrix, should be pre-calloced
/// \param obs a list of all orbitals
/// \param atoms a list of all atoms
void core_hamiltonian_matrix_set(gsl_matrix * dest,std::vector<Orbital> &obs,const std::vector<Atom> &atoms){
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
double nuclei_repulsion(std::vector<Atom> &atoms){
    if (atoms.size() <=1){
        return 0;
    }
    double result = 0;
    double distance;
    for (auto &a:atoms){
        for(auto &b:atoms){
            if (&a !=&b){
                distance = sqrt(pow(a.cartesian[0]-b.cartesian[0],2)+pow(a.cartesian[1]-a.cartesian[1],2)+pow(a.cartesian[2],b.cartesian[2],2));
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