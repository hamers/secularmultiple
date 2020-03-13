/*
*/

#include "types.h"
#include "postnewtonian.h"
#include <stdio.h>
extern "C"
{
void compute_EOM_pairwise_1PN(ParticlesMap *particlesMap, int binary_index, double *hamiltonian, bool compute_hamiltonian_only)
{
    Particle *binary = (*particlesMap)[binary_index];
    double e = binary->e;
    double a = binary->a;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;
    double mt = m1+m2;

    double hamiltonian_1PN = -3.0*CONST_G_P2*m1*m2*mt/(a*a*CONST_C_LIGHT_P2*j);
    *hamiltonian += hamiltonian_1PN;
        
    if (compute_hamiltonian_only == true)
    {
        return;
    }
    
    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);
    
    double GMdiva = CONST_G*mt/a;
    double Z_1PN = 3.0*sqrt(GMdiva)*GMdiva/(a*CONST_C_LIGHT_P2*j_p2);
    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += e*Z_1PN*q_vec_unit[i];
    }

    #ifdef DEBUG
    printf("postnewtonian.cpp -- compute_EOM_pairwise_1PN -- hamiltonian_1PN %g de_vec_dt += %g %g %g\n",hamiltonian_1PN,e*Z_1PN*q_vec_unit[0],e*Z_1PN*q_vec_unit[1],e*Z_1PN*q_vec_unit[2]);
    #endif
}

void compute_EOM_pairwise_25PN(ParticlesMap *particlesMap, int binary_index, double *hamiltonian, bool compute_hamiltonian_only)
{
    Particle *binary = (*particlesMap)[binary_index];
    double e = binary->e;
    double e_p2 = e*e;
    double a = binary->a;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p4 = binary->j_p4;
    double m1 = binary->child1_mass;
    double m2 = binary->child2_mass;
    double mt = m1+m2;

    double a_p3 = a*a*a;
    double GMdiva = CONST_G*mt/a;
    double c_common = CONST_G_P3*m1*m2/(CONST_C_LIGHT_P5*a_p3*j_p4);
    double f_e = 1.0 + c_121div304*e_p2;
    double f_h = 1.0 + c_7div8*e_p2;

    double de_dt = -c_304div15*c_common*mt*e*f_e/(a*j);
    double dh_dt = -c_32div5*c_common*m1*m2*sqrt(GMdiva)*f_h;

    for (int i=0; i<3; i++)
    {
        binary->de_vec_dt[i] += de_dt*e_vec_unit[i];
        binary->dh_vec_dt[i] += dh_dt*h_vec_unit[i];
    }

    #ifdef DEBUG
    printf("postnewtonian.cpp -- compute_EOM_pairwise_25PN -- de_vec_dt += %g %g %g dh_vec_dt += %g %g %g\n",de_dt*e_vec_unit[0],de_dt*e_vec_unit[1],de_dt*e_vec_unit[2],dh_dt*h_vec_unit[0],dh_dt*h_vec_unit[1],dh_dt*h_vec_unit[2]);
    #endif
    
}
}
