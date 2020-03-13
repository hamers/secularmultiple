#include "types.h"

extern "C"
{
    
int highest_particle_index = 0;
ParticlesMap particlesMap;

// Default constants //
double CONST_G = 4.0*M_PI*M_PI; 
double CONST_G_P2 = CONST_G*CONST_G;
double CONST_G_P3 = CONST_G_P2*CONST_G;
double CONST_C_LIGHT = 63239.72638679138;
double CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
double CONST_C_LIGHT_P4 = CONST_C_LIGHT_P2*CONST_C_LIGHT_P2;
double CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;
double CONST_MSUN = 1.0;
double CONST_R_SUN = 0.004649130343817401;
double CONST_L_SUN = 0.0002710404109745588;
double CONST_KM_PER_S = 0.210862;
double CONST_PER_PC3 = 1.14059e-16;

// Default parameters //
double relative_tolerance = 1.0e-12;
double absolute_tolerance_eccentricity_vectors = 1.0e-10;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = true;
bool include_hexadecupole_order_binary_pair_terms = true;
bool include_dotriacontupole_order_binary_pair_terms = true;
bool include_double_averaging_corrections = false;
int orbital_phases_random_seed = 0;

double epsilon = 1.0e-15; /* used for tiny numbers close to machine precision */

}
