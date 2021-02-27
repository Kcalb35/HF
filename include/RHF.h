_Pragma("once")

#include <gsl/gsl_matrix.h>
#include "basis.h"
#include "gslextra.h"
#include "input.h"
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#define THREADS_NUMBER 6

static boost::asio::thread_pool pool(THREADS_NUMBER);

double single_electron_hamiltonian_element(Orbital &LOb, Orbital &ROb, std::vector<Atom> &atoms);
void core_hamiltonian_matrix_set(gsl_matrix * dest, std::vector<Orbital> &obs, std::vector<Atom> &atoms);

void initial_guess(gsl_matrix * coeff_matrix,gsl_matrix * core_hamiltonian_matrix,gsl_matrix * overlap_matrix);
void density_matrix_set(gsl_matrix * density_matrix,gsl_matrix * coeff_matrix,int electron_number);

double Fock_matrix_element(gsl_quad_tensor * tensor,gsl_matrix * density_matrix,gsl_matrix * h_core_matrix,int i,int j);
void Fock_matrix(gsl_matrix * fock,gsl_quad_tensor * v,gsl_matrix * density,gsl_matrix * h_core);

double nuclei_repulsion(const std::vector<Atom> &atoms);
double HF_energy(gsl_quad_tensor *v,gsl_matrix * density,gsl_matrix * h_core);
void two_electron_quad_tensor_set(gsl_quad_tensor * dest,std::vector<Orbital> &obs);

void kinetic_matrix_set(gsl_matrix *m,std::vector<Orbital> &obs);
double nuclear_attraction_energy(Orbital &LOb, Orbital &ROb, std::vector<Atom> &atoms);
void nuclear_attraction_energy_matrix_set(gsl_matrix *m,std::vector<Orbital> obs,std::vector<Atom> &atoms);
void Orbital_S_matrix_set(gsl_matrix * S,std::vector<Orbital>& orbitals);
void RHF_SCF(double * total_energy, gsl_vector * energyLevel, gsl_matrix * coeff, std::vector<Atom> &allAtoms, std::vector<Orbital> &allOrbitals, const HFoption &option);

