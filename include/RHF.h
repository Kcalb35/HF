_Pragma("once")

#include <gsl/gsl_matrix.h>
#include "basis.h"
#include "gslextra.h"

double single_electron_hamiltonian_element(Orbital &LOb,Orbital &ROb, const std::vector<Atom> &atoms);
void core_hamiltonian_matrix_set(gsl_matrix * dest,std::vector<Orbital> &obs,const std::vector<Atom> &atoms);

void initial_guess(gsl_matrix * coeff_matrix,gsl_matrix * core_hamiltonian_matrix,gsl_matrix * overlap_matrix);
void density_matrix_set(gsl_matrix * density_matrix,gsl_matrix * coeff_matrix,int electron_number);

double Fock_matrix_element(gsl_quad_tensor * tensor,gsl_matrix * density_matrix,gsl_matrix * h_core_matrix,int i,int j);
void Fock_matrix(gsl_matrix * fock,gsl_quad_tensor * v,gsl_matrix * density,gsl_matrix * h_core);

double nuclei_repulsion(std::vector<Atom> &atoms);
double HF_energy(gsl_quad_tensor *v,gsl_matrix * density,gsl_matrix * h_core);
