#include "types.h"
extern "C"
{
void compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, double *hamiltonian, bool compute_hamiltonian_only);
void compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, double *hamiltonian, bool compute_hamiltonian_only);
}
