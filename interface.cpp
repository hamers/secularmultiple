#include <stdio.h>
#include <math.h>
#include "interface.h"

#include "src/types.h"
//#include "interface.h"
#include "src/evolve.h"


extern "C"
{
    
int highest_particle_index = 0;
ParticlesMap particlesMap;

double relative_tolerance = 1.0e-16;
double absolute_tolerance_eccentricity_vectors = 1.0e-14;
bool include_quadrupole_order_terms = true;
bool include_octupole_order_binary_pair_terms = true;
bool include_octupole_order_binary_triplet_terms = false;
bool include_hexadecupole_order_binary_pair_terms = false;
bool include_dotriacontupole_order_binary_pair_terms = false;
int orbital_phases_random_seed = 0;


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

/*******************
/* basic interface *
 ******************/
 
 
int add_particle(int *index, int is_binary, int is_external)
{
    *index = highest_particle_index;
    Particle *p = new Particle(highest_particle_index, is_binary);
    particlesMap[highest_particle_index] = p;

    p->is_external = is_external;

    highest_particle_index++;
    return 0;
}


int delete_particle(int index)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    particlesMap.erase(index);

    return 0;
}

int set_children(int index, int child1, int child2)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    p->child1 = child1;
    p->child2 = child2;
    //printf("c1 %d c2 %d\n",child1,child2);
    return 0;
}
int get_children(int index, int *child1, int *child2)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *child1 = p->child1;
    *child2 = p->child2;
    
    return 0;
}

int set_mass(int index, double mass)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    p->mass = mass;

    return 0;
}
int get_mass(int index, double *mass)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
    Particle *p = particlesMap[index];
    *mass = p->mass;
    return 0;
}
int set_mass_dot(int index, double mass_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    p->mass_dot = mass_dot;

    return 0;
}
int get_mass_dot(int index, double *mass_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle *p = particlesMap[index];
    *mass_dot = p->mass_dot;

    return 0;
}

int set_radius(int index, double radius, double radius_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->radius = radius;
    p->radius_dot = radius_dot;
    
    return 0;
}
int get_radius(int index, double *radius, double *radius_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *radius = p->radius;
    *radius_dot = p->radius_dot;
    
    return 0;
}


/****************************
/* orbital vectors/elements *
 ****************************/

int set_orbital_vectors(int index, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];
    p->e_vec_x = e_vec_x;
    p->e_vec_y = e_vec_y;
    p->e_vec_z = e_vec_z;
    p->h_vec_x = h_vec_x;
    p->h_vec_y = h_vec_y;
    p->h_vec_z = h_vec_z;
    //printf("ok %g %g %g %g %g %g\n",e_vec_x,e_vec_y,e_vec_z,h_vec_x,h_vec_y,h_vec_z);
    return 0;
}
int get_orbital_vectors(int index, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle *p = particlesMap[index];
    *e_vec_x = p->e_vec_x;
    *e_vec_y = p->e_vec_y;
    *e_vec_z = p->e_vec_z;
    *h_vec_x = p->h_vec_x;
    *h_vec_y = p->h_vec_y;
    *h_vec_z = p->h_vec_z;

    return 0;
}

int set_orbital_elements(int index, double semimajor_axis, double eccentricity, double true_anomaly, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node, int sample_orbital_phase_randomly)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index];    
    
    if (p->is_binary == false)
    {
        return 0;
    }

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_vectors_from_orbital_elements(p->child1_mass, p->child2_mass, semimajor_axis, eccentricity, \
        inclination, argument_of_pericenter, longitude_of_ascending_node, \
        &(p->e_vec_x), &(p->e_vec_y), &(p->e_vec_z), &(p->h_vec_x), &(p->h_vec_y), &(p->h_vec_z) );
    
    p->true_anomaly = true_anomaly;
    p->sample_orbital_phase_randomly = sample_orbital_phase_randomly;
//    printf("soe a %g e %g TA %g I %g AP %g LAN %g SOPR %d\n",semimajor_axis,eccentricity,true_anomaly,inclination,argument_of_pericenter,longitude_of_ascending_node,sample_orbital_phase_randomly);
//    printf("set_orbital_elements %g %g %g %g %g %g\n",p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z);
    
    return 0;
}
int get_orbital_elements(int index, double *semimajor_axis, double *eccentricity, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];    
    
    if (p->is_binary == false)
    {
        return 0;
    }

    double h_tot_vec[3];
    compute_h_tot_vector(&particlesMap,h_tot_vec);

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
        p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z,
        semimajor_axis, eccentricity, inclination, argument_of_pericenter, longitude_of_ascending_node);
    return 0;
}


int get_level(int index, int *value)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *value = p->level;
    
    return 0;
}

int set_stellar_type(int index, int value)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->stellar_type = value;
    
    return 0;
}
int get_stellar_type(int index, int *value)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *value = p->stellar_type;
    
    return 0;
}



/*******************************
 * instantaneous perturbations *
 * ****************************/

int set_instantaneous_perturbation_properties(int index, double delta_mass, double delta_x, double delta_y, double delta_z, double delta_vx, double delta_vy, double delta_vz)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
    
    Particle *p = particlesMap[index];

    p->instantaneous_perturbation_delta_mass = delta_mass;
    p->instantaneous_perturbation_delta_x = delta_x;
    p->instantaneous_perturbation_delta_y = delta_y;
    p->instantaneous_perturbation_delta_z = delta_z;
    p->instantaneous_perturbation_delta_vx = delta_vx;
    p->instantaneous_perturbation_delta_vy = delta_vy;
    p->instantaneous_perturbation_delta_vz = delta_vz;
    
    return 0;
}


/************
 * external *
 * *********/

int set_external_particle_properties(int index, double external_t_ref, double e, double external_r_p, double INCL, double AP, double LAN)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];    
    
    /* determine masses in all binaries */
//    int N_bodies, N_binaries, N_root_finding;
//    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
//    set_binary_masses_from_body_masses(&particlesMap);
    
    p->external_t_ref = external_t_ref;
    p->external_e = e;
    p->external_r_p = external_r_p;
    /* e & h vectors for external particles are understood to be unit vectors */
    compute_orbital_vectors_from_orbital_elements_unit(INCL,AP,LAN,&(p->e_vec_x), &(p->e_vec_y), &(p->e_vec_z), &(p->h_vec_x), &(p->h_vec_y), &(p->h_vec_z) ); 
    
    //printf("set_external_particle_properties inputs %g %g %g %g %g\n",external_t_ref, e, external_r_p, INCL, AP, LAN);
    //printf("set_external_particle_properties OE %g %g %g %g %g %g\n",p->e_vec_x,p->e_vec_y,p->e_vec_z,p->h_vec_x,p->h_vec_y,p->h_vec_z);
    
    return 0;
}




/****************
/* spin vectors *
 ****************/

int set_spin_vector(int index, double spin_vec_x, double spin_vec_y, double spin_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->spin_vec_x = spin_vec_x;
    p->spin_vec_y = spin_vec_y;
    p->spin_vec_z = spin_vec_z;
    //printf("set spin %g %g %g\n",spin_vec_x,spin_vec_y,spin_vec_z);
    return 0;
}
int get_spin_vector(int index, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x = p->spin_vec_x;
    *spin_vec_y = p->spin_vec_y;
    *spin_vec_z = p->spin_vec_z;
    
    return 0;
}

int set_spin_vector_dot(int index, double spin_vec_x_dot, double spin_vec_y_dot, double spin_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->spin_vec_x_dot = spin_vec_x_dot;
    p->spin_vec_y_dot = spin_vec_y_dot;
    p->spin_vec_z_dot = spin_vec_z_dot;
    
    return 0;
}
int get_spin_vector_dot(int index, double *spin_vec_x_dot, double *spin_vec_y_dot, double *spin_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x_dot = p->spin_vec_x_dot;
    *spin_vec_y_dot = p->spin_vec_y_dot;
    *spin_vec_z_dot = p->spin_vec_z_dot;
    
    return 0;
}


int set_orbital_vectors_dot(int index, double e_vec_x_dot, double e_vec_y_dot, double e_vec_z_dot, \
    double h_vec_x_dot, double h_vec_y_dot, double h_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->e_vec_x_dot = e_vec_x_dot;
    p->e_vec_y_dot = e_vec_y_dot;
    p->e_vec_z_dot = e_vec_z_dot;
    p->h_vec_x_dot = h_vec_x_dot;
    p->h_vec_y_dot = h_vec_y_dot;
    p->h_vec_z_dot = h_vec_z_dot;
    
    return 0;
}
int get_orbital_vectors_dot(int index, double *e_vec_x_dot, double *e_vec_y_dot, double *e_vec_z_dot, \
    double *h_vec_x_dot, double *h_vec_y_dot, double *h_vec_z_dot)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *e_vec_x_dot = p->e_vec_x_dot;
    *e_vec_y_dot = p->e_vec_y_dot;
    *e_vec_z_dot = p->e_vec_z_dot;
    *h_vec_x_dot = p->h_vec_x_dot;
    *h_vec_y_dot = p->h_vec_y_dot;
    *h_vec_z_dot = p->h_vec_z_dot;
    
    return 0;
}


int set_position_vector(int index, double x, double y, double z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->x = x;
    p->x = y;
    p->x = z;
   
    return 0;
}
int get_position_vector(int index, double *x, double *y, double *z)
{
    //printf("get_position_vector\n");
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    set_positions_and_velocities(&particlesMap);
    
    Particle * p = particlesMap[index];
    *x = p->x;
    *y = p->x;
    *z = p->x;
    
    return 0;
}

int set_velocity_vector(int index, double x, double y, double z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    p->vx = x;
    p->vy = y;
    p->vz = z;
   
    return 0;
}
int get_velocity_vector(int index, double *x, double *y, double *z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    set_positions_and_velocities(&particlesMap);
    
    Particle * p = particlesMap[index];
    *x = p->vx;
    *y = p->vy;
    *z = p->vz;
    
    return 0;
}

/************
/* PN terms *
 ************/

int set_PN_terms(int index, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle * p = particlesMap[index];
    p->include_pairwise_1PN_terms = include_pairwise_1PN_terms;
    p->include_pairwise_25PN_terms = include_pairwise_25PN_terms;
        
    return 0;
}
int get_PN_terms(int index, int *include_pairwise_1PN_terms, int *include_pairwise_25PN_terms)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];
    
    *include_pairwise_1PN_terms = p->include_pairwise_1PN_terms;
    *include_pairwise_25PN_terms = p->include_pairwise_25PN_terms;
        
    return 0;
}


/*********
/* tides *
 *********/
int set_tides_terms(int index, int include_tidal_friction_terms, int tides_method, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, 
    double tides_apsidal_motion_constant, double tides_gyration_radius, double tides_viscous_time_scale, int tides_viscous_time_scale_prescription, double convective_envelope_mass, double convective_envelope_radius, double luminosity)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    
    p->include_tidal_friction_terms = include_tidal_friction_terms;
    p->tides_method = tides_method;
    p->include_tidal_bulges_precession_terms = include_tidal_bulges_precession_terms;
    p->include_rotation_precession_terms = include_rotation_precession_terms;
    p->minimum_eccentricity_for_tidal_precession = minimum_eccentricity_for_tidal_precession;
    p->tides_apsidal_motion_constant = tides_apsidal_motion_constant;
    p->tides_gyration_radius = tides_gyration_radius;
    p->tides_viscous_time_scale = tides_viscous_time_scale;
    p->tides_viscous_time_scale_prescription = tides_viscous_time_scale_prescription;
    p->convective_envelope_mass = convective_envelope_mass;
    p->convective_envelope_radius = convective_envelope_radius;
    p->luminosity = luminosity;
    //printf("set tides1 %d %d %d %g\n",include_tidal_friction_terms,include_tidal_bulges_precession_terms,include_rotation_precession_terms,minimum_eccentricity_for_tidal_precession);
    //printf("set tides2 %g %g %g %d %g %g %g\n",tides_apsidal_motion_constant,tides_gyration_radius,tides_viscous_time_scale,tides_viscous_time_scale_prescription,convective_envelope_mass,convective_envelope_radius,luminosity);
    return 0;
}
int get_tides_terms(int index, int *include_tidal_friction_terms, int *tides_method, int *include_tidal_bulges_precession_terms, int *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession, 
    double *tides_apsidal_motion_constant, double *tides_gyration_radius, double *tides_viscous_time_scale, int *tides_viscous_time_scale_prescription, double *convective_envelope_mass, double *convective_envelope_radius, double *luminosity)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle *p = particlesMap[index];
    
    *include_tidal_friction_terms = p->include_tidal_friction_terms;
    *tides_method = p->tides_method;
    *include_tidal_bulges_precession_terms = p->include_tidal_bulges_precession_terms;
    *include_rotation_precession_terms = p->include_rotation_precession_terms;
    *minimum_eccentricity_for_tidal_precession = p->minimum_eccentricity_for_tidal_precession;
    *tides_apsidal_motion_constant = p->tides_apsidal_motion_constant;
    *tides_gyration_radius = p->tides_gyration_radius;
    *tides_viscous_time_scale = p->tides_viscous_time_scale;
    *tides_viscous_time_scale_prescription = p->tides_viscous_time_scale_prescription;
    *convective_envelope_mass = p->convective_envelope_mass;
    *convective_envelope_radius = p->convective_envelope_radius;
    *luminosity = p->luminosity;
    return 0;
}

/****************
 * VRR          *
 ****************/

int set_VRR_properties(int index, int VRR_model, int VRR_include_mass_precession, double VRR_mass_precession_rate, 
    double VRR_Omega_vec_x, double VRR_Omega_vec_y, double VRR_Omega_vec_z, 
    double VRR_eta_20_init, double VRR_eta_a_22_init, double VRR_eta_b_22_init, double VRR_eta_a_21_init, double VRR_eta_b_21_init,
    double VRR_eta_20_final, double VRR_eta_a_22_final, double VRR_eta_b_22_final, double VRR_eta_a_21_final, double VRR_eta_b_21_final,
	double VRR_initial_time, double VRR_final_time)
{
	
	if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    
    p->VRR_model = VRR_model;
    p->VRR_include_mass_precession = VRR_include_mass_precession;
	p->VRR_mass_precession_rate = VRR_mass_precession_rate;
	p->VRR_Omega_vec_x = VRR_Omega_vec_x;
	p->VRR_Omega_vec_y = VRR_Omega_vec_y;
	p->VRR_Omega_vec_z = VRR_Omega_vec_z;
	p->VRR_eta_20_init = VRR_eta_20_init;
	p->VRR_eta_a_22_init = VRR_eta_a_22_init;
	p->VRR_eta_b_22_init = VRR_eta_b_22_init;
	p->VRR_eta_a_21_init = VRR_eta_a_21_init;
	p->VRR_eta_b_21_init = VRR_eta_b_21_init;
	p->VRR_eta_20_final = VRR_eta_20_final;
	p->VRR_eta_a_22_final = VRR_eta_a_22_final;
	p->VRR_eta_b_22_final = VRR_eta_b_22_final;
	p->VRR_eta_a_21_final = VRR_eta_a_21_final;
	p->VRR_eta_b_21_final = VRR_eta_b_21_final;
	p->VRR_initial_time = VRR_initial_time;
	p->VRR_final_time = VRR_final_time;

	return 0;
}

/****************
/* root finding *
 ****************/
int set_root_finding_terms(int index, int check_for_secular_breakdown, int check_for_dynamical_instability, int dynamical_instability_criterion, int dynamical_instability_central_particle, int dynamical_instability_K_parameter,
    int check_for_physical_collision_or_orbit_crossing, int check_for_minimum_periapse_distance, double check_for_minimum_periapse_distance_value, int check_for_RLOF_at_pericentre, int check_for_RLOF_at_pericentre_use_sepinsky_fit, int check_for_GW_condition)
{
    if (index > highest_particle_index)
    {
        return -1;
    }

    Particle *p = particlesMap[index];
    p->check_for_secular_breakdown = check_for_secular_breakdown;
    p->check_for_dynamical_instability = check_for_dynamical_instability;
    p->dynamical_instability_criterion = dynamical_instability_criterion;
    p->dynamical_instability_central_particle = dynamical_instability_central_particle;
    p->dynamical_instability_K_parameter = dynamical_instability_K_parameter;
    p->check_for_physical_collision_or_orbit_crossing = check_for_physical_collision_or_orbit_crossing;
    p->check_for_minimum_periapse_distance = check_for_minimum_periapse_distance;
    p->check_for_minimum_periapse_distance_value = check_for_minimum_periapse_distance_value;
    p->check_for_RLOF_at_pericentre = check_for_RLOF_at_pericentre;
    p->check_for_RLOF_at_pericentre_use_sepinsky_fit = check_for_RLOF_at_pericentre_use_sepinsky_fit;
    p->check_for_GW_condition = check_for_GW_condition;
    return 0;
}
int get_root_finding_terms(int index, int *check_for_secular_breakdown, int *check_for_dynamical_instability, int *dynamical_instability_criterion, int *dynamical_instability_central_particle, int *dynamical_instability_K_parameter,
    int *check_for_physical_collision_or_orbit_crossing, int *check_for_minimum_periapse_distance, double *check_for_minimum_periapse_distance_value, int *check_for_RLOF_at_pericentre, int *check_for_RLOF_at_pericentre_use_sepinsky_fit, int *check_for_GW_condition)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle *p = particlesMap[index];
    *check_for_secular_breakdown = p->check_for_secular_breakdown;
    *check_for_dynamical_instability = p->check_for_dynamical_instability;
    *dynamical_instability_criterion = p->dynamical_instability_criterion;
    *dynamical_instability_central_particle = p->dynamical_instability_central_particle;
    *dynamical_instability_K_parameter = p->dynamical_instability_K_parameter;
    *check_for_physical_collision_or_orbit_crossing = p->check_for_physical_collision_or_orbit_crossing;
    *check_for_minimum_periapse_distance = p->check_for_minimum_periapse_distance;
    *check_for_minimum_periapse_distance_value = p->check_for_minimum_periapse_distance_value;
    *check_for_RLOF_at_pericentre = p->check_for_RLOF_at_pericentre;
    *check_for_RLOF_at_pericentre_use_sepinsky_fit = p->check_for_RLOF_at_pericentre_use_sepinsky_fit;
    *check_for_GW_condition = p->check_for_GW_condition;
    return 0;
}

/* retrieve root finding state */
int set_root_finding_state(int index, int secular_breakdown_has_occurred, int dynamical_instability_has_occurred, int physical_collision_or_orbit_crossing_has_occurred, int minimum_periapse_distance_has_occurred, int RLOF_at_pericentre_has_occurred, int GW_condition_has_occurred)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];

    p->secular_breakdown_has_occurred = secular_breakdown_has_occurred;
    p->dynamical_instability_has_occurred = dynamical_instability_has_occurred;
    p->physical_collision_or_orbit_crossing_has_occurred = physical_collision_or_orbit_crossing_has_occurred;
    p->minimum_periapse_distance_has_occurred = minimum_periapse_distance_has_occurred;
    p->RLOF_at_pericentre_has_occurred = RLOF_at_pericentre_has_occurred;
    p->GW_condition_has_occurred = GW_condition_has_occurred;
    
    return 0;
}
int get_root_finding_state(int index, int *secular_breakdown_has_occurred, int *dynamical_instability_has_occurred, int *physical_collision_or_orbit_crossing_has_occurred, int* minimum_periapse_distance_has_occurred, int *RLOF_at_pericentre_has_occurred, int *GW_condition_has_occurred)
{
    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];

    *secular_breakdown_has_occurred = p->secular_breakdown_has_occurred;
    *dynamical_instability_has_occurred = p->dynamical_instability_has_occurred;
    *physical_collision_or_orbit_crossing_has_occurred = p->physical_collision_or_orbit_crossing_has_occurred;
    *minimum_periapse_distance_has_occurred = p->minimum_periapse_distance_has_occurred;
    *RLOF_at_pericentre_has_occurred = p->RLOF_at_pericentre_has_occurred;
    *GW_condition_has_occurred = p->GW_condition_has_occurred;
    
    return 0;
}


/********************
/* evolve interface *
 ********************/

int evolve_interface(double start_time, double time_step, double *output_time, double *hamiltonian, int *flag, int *error_code)
{
    int result = evolve(&particlesMap, start_time, time_step, output_time, hamiltonian, flag, error_code);
    
    return result;
}


/* set levels and masses */
int determine_binary_parents_levels_and_masses_interface()
{
    //printf("determine_binary_parents_levels_and_masses_interface\n");
    int N_bodies, N_binaries, N_root_finding;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding);
    set_binary_masses_from_body_masses(&particlesMap);
    
    return 0;
}

int apply_external_perturbation_assuming_integrated_orbits_interface()
{
    //printf("apply_external_perturbation_assuming_integrated_orbits_interface\n");
    apply_external_perturbation_assuming_integrated_orbits(&particlesMap);

    return 0;
}

int apply_user_specified_instantaneous_perturbation_interface()
{
    //printf("apply_user_specified_instantaneous_perturbation\n");
    apply_user_specified_instantaneous_perturbation(&particlesMap);
    
    return 0;
}


int clear_internal_particles()
{
    //printf("clear_internal_particles\n");
    particlesMap.clear();
	highest_particle_index = 0;
    return 0;
}


int set_positions_and_velocities_interface()
{
    set_positions_and_velocities(&particlesMap);
}

/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vec[3])
{
    for (int i=0; i<3; i++)
    {
        h_tot_vec[i] = 0.0;
    }
    
    ParticlesMapIterator it_p;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *p = (*it_p).second;
        if (p->is_binary == 1)
        {
            h_tot_vec[0] += p->h_vec_x;
            h_tot_vec[1] += p->h_vec_y;
            h_tot_vec[2] += p->h_vec_z;
        }
    }
    return 0;
//    printf("compute_h_tot_vector %g %g %g\n",h_tot_vec[0],h_tot_vec[1],h_tot_vec[2]);
}
    
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z)
{
    double cos_INCL = cos(inclination);
    double sin_INCL = sin(inclination);
    double cos_AP = cos(argument_of_pericenter);
    double sin_AP = sin(argument_of_pericenter);
    double cos_LAN = cos(longitude_of_ascending_node);
    double sin_LAN = sin(longitude_of_ascending_node);
           
    double h = (child1_mass*child2_mass*sqrt(CONST_G*semimajor_axis/(child1_mass+child2_mass)))*sqrt(1.0 - eccentricity*eccentricity);

    *e_vec_x = eccentricity*(cos_LAN*cos_AP - sin_LAN*sin_AP*cos_INCL);
    *e_vec_y = eccentricity*(sin_LAN*cos_AP + cos_LAN*sin_AP*cos_INCL);
    *e_vec_z = eccentricity*(sin_AP*sin_INCL);
    
    *h_vec_x = h*sin_LAN*sin_INCL;
    *h_vec_y = -h*cos_LAN*sin_INCL;
    *h_vec_z = h*cos_INCL;

    return 0;
}

int compute_orbital_vectors_from_orbital_elements_unit(double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_hat_vec_x, double *e_hat_vec_y, double *e_hat_vec_z, double *h_hat_vec_x, double *h_hat_vec_y, double *h_hat_vec_z)
{
    double cos_INCL = cos(inclination);
    double sin_INCL = sin(inclination);
    double cos_AP = cos(argument_of_pericenter);
    double sin_AP = sin(argument_of_pericenter);
    double cos_LAN = cos(longitude_of_ascending_node);
    double sin_LAN = sin(longitude_of_ascending_node);
           
    *e_hat_vec_x = (cos_LAN*cos_AP - sin_LAN*sin_AP*cos_INCL);
    *e_hat_vec_y = (sin_LAN*cos_AP + cos_LAN*sin_AP*cos_INCL);
    *e_hat_vec_z = (sin_AP*sin_INCL);
    
    *h_hat_vec_x = sin_LAN*sin_INCL;
    *h_hat_vec_y = -cos_LAN*sin_INCL;
    *h_hat_vec_z = cos_INCL;

    return 0;
}

int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node)
{
    double e_vec[3] = {e_vec_x,e_vec_y,e_vec_z};
    double h_vec[3] = {h_vec_x,h_vec_y,h_vec_z};
    double eccentricity_squared = norm3_squared(e_vec);
    *eccentricity = sqrt(eccentricity_squared);
    double h_squared = norm3_squared(h_vec);
    *semimajor_axis = h_squared*(child1_mass+child2_mass)/( CONST_G*child1_mass*child1_mass*child2_mass*child2_mass*(1.0 - eccentricity_squared) );
    double h = sqrt(h_squared);
    
//    double x_vec[3] = {1.0,0.0,0.0};
//    double y_vec[3] = {0.0,1.0,0.0};
//    double z_vec[3] = {0.0,0.0,1.0};

    double h_tot = norm3(h_tot_vec);
//    printf("h_tot %g x %g y %g z %g\n",h_tot,h_tot_vec[0],h_tot_vec[1],h_tot_vec[2]);
    double x_vec[3], y_vec[3], z_vec[3];
    for (int i=0; i<3; i++)
    {
        z_vec[i] = h_tot_vec[i]/h_tot;
    }

//    printf("test %g %g %g\n",z_vec[0],z_vec[1],z_vec[2]);
    z_vec[0] = 0.0;
    z_vec[1] = 0.0;
    z_vec[2] = 1.0;

    /* the above assumes that the total angular momentum vector does not change (i.e. no SNe effects etc.) */
    
    double f = 1.0/sqrt( z_vec[0]*z_vec[0] + z_vec[2]*z_vec[2] );
    x_vec[0] = z_vec[2]*f;
    x_vec[1] = 0.0;
    x_vec[2] = -z_vec[0]*f;
    cross3(z_vec,x_vec,y_vec);

    double cos_INCL = dot3(h_vec,z_vec)/h;

    double LAN_vec[3],LAN_vec_unit[3];
    cross3(z_vec,h_vec,LAN_vec);
    double LAN_vec_norm = norm3(LAN_vec);

    double e_vec_unit[3],h_vec_unit[3];

    for (int i=0; i<3; i++)
    {
        LAN_vec_unit[i] = LAN_vec[i]/LAN_vec_norm;
        e_vec_unit[i] = e_vec[i]/(*eccentricity);
        h_vec_unit[i] = h_vec[i]/h;
    }

    double sin_LAN = dot3(LAN_vec_unit,y_vec);
    double cos_LAN = dot3(LAN_vec_unit,x_vec);

    double e_vec_unit_cross_h_vec_unit[3];
    cross3(e_vec_unit,h_vec_unit,e_vec_unit_cross_h_vec_unit);
    double sin_AP = dot3(LAN_vec_unit,e_vec_unit_cross_h_vec_unit);
    double cos_AP = dot3(LAN_vec_unit,e_vec_unit);

    *inclination = acos(cos_INCL);
    *argument_of_pericenter = atan2(sin_AP,cos_AP);
    *longitude_of_ascending_node = atan2(sin_LAN,cos_LAN);
    
    return 0;
}

int get_inclination_relative_to_parent(int index, double *inclination_relative_to_parent)
{

    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];
    if (p->is_binary == 0)
    {
        *inclination_relative_to_parent = 0.0;
        return 0;
    }
    if (p->parent == -1)
    {
        *inclination_relative_to_parent = 0.0;
        return 0;
    }

    Particle *parent = particlesMap[p->parent];
    

    double h1_vec[3] = {p->h_vec_x,p->h_vec_y,p->h_vec_z};
    double h2_vec[3] = {parent->h_vec_x,parent->h_vec_y,parent->h_vec_z};

    double h1 = norm3(h1_vec);
    double h2 = norm3(h2_vec);
    
    *inclination_relative_to_parent = acos( dot3(h1_vec,h2_vec)/(h1*h2) );
    
    return 0;
}


void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly)
{
    double eccentric_anomaly;
    double eccentric_anomaly_next = mean_anomaly;
    double epsilon = 1e-10;
    double error = 2.0*epsilon;
    int j = 0;
    while (error > epsilon || j < 15)
    {
        j += 1;
        eccentric_anomaly = eccentric_anomaly_next;
        eccentric_anomaly_next = eccentric_anomaly - (eccentric_anomaly - eccentricity*sin(eccentric_anomaly) - mean_anomaly)/(1.0 - eccentricity*cos(eccentric_anomaly));
        error = fabs(eccentric_anomaly_next - eccentric_anomaly);
    }
    *cos_eccentric_anomaly = cos(eccentric_anomaly);
    *sin_eccentric_anomaly = sin(eccentric_anomaly);
}

void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly)
{
    *cos_true_anomaly = (cos_eccentric_anomaly - eccentricity)/(1.0 - eccentricity*cos_eccentric_anomaly);
    *sin_true_anomaly = sqrt(1.0 - eccentricity*eccentricity)*sin_eccentric_anomaly/(1.0 - eccentricity*cos_eccentric_anomaly);
}

double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity)
{
    double cos_eccentric_anomaly,sin_eccentric_anomaly;
    double cos_true_anomaly,sin_true_anomaly;
    
    compute_eccentric_anomaly_from_mean_anomaly(mean_anomaly,eccentricity,&cos_eccentric_anomaly,&sin_eccentric_anomaly);
    compute_true_anomaly_from_eccentric_anomaly(cos_eccentric_anomaly,sin_eccentric_anomaly,eccentricity,&cos_true_anomaly,&sin_true_anomaly);
    double true_anomaly = atan2(sin_true_anomaly,cos_true_anomaly);

    return true_anomaly;
}

double sample_random_true_anomaly(double eccentricity,int seed)
{
    srand(seed);
    double x = ((double) rand() / (RAND_MAX));
    double mean_anomaly = (2.0*x - 1.0)*M_PI;
    double true_anomaly = compute_true_anomaly_from_mean_anomaly(mean_anomaly,eccentricity);

    return true_anomaly;
}

void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3])
{
    double total_mass = child1_mass + child2_mass;
    
    double e = norm3(e_vec);
    double h = norm3(h_vec);

    double e_vec_unit[3],q_vec_unit[3],q_vec[3];
    cross3(h_vec,e_vec,q_vec);
    double q = norm3(q_vec);

    int i;
    for (i=0; i<3; i++)
    {        
        e_vec_unit[i] = e_vec[i]/e;
        q_vec_unit[i] = q_vec[i]/q;        
    }
    
    double e_p2 = e*e;
    double j_p2 = 1.0 - e_p2;
   
    double a = h*h*total_mass/( CONST_G*child1_mass*child2_mass*child1_mass*child2_mass*j_p2 );

    double cos_f = cos(true_anomaly);
    double sin_f = sin(true_anomaly);
    
    double r_norm = a*j_p2/(1.0 + e*cos_f);
    double v_norm = sqrt( CONST_G*total_mass/(a*j_p2) );
    
    for (i=0; i<3; i++)
    {
        r[i] = r_norm*( cos_f*e_vec_unit[i] + sin_f*q_vec_unit[i]);
        v[i] = v_norm*( -sin_f*e_vec_unit[i] + (e + cos_f)*q_vec_unit[i] );
    }
}

void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3])
{
    double total_mass = child1_mass + child2_mass;
       
    double v_dot_v = dot3(v,v);
    double r_dot_v = dot3(r,v);
    double r_norm = norm3(r);
    for (int i=0; i<3; i++)
    {
        e_vec[i] = (r[i]*v_dot_v - v[i]*r_dot_v)/(CONST_G*total_mass) - r[i]/r_norm;
    }

    double mu = child1_mass*child2_mass/total_mass;
    cross3(r,v,h_vec);
    for (int i=0; i<3; i++)
    {
        h_vec[i] *= mu;
    }
}


int get_de_dt(int index, double *de_dt)
{

    if (index > highest_particle_index)
    {
        return -1;
    }
  
    Particle * p = particlesMap[index];
    if (p->is_binary == 0)
    {
        *de_dt = 0.0;
        return 0;
    }

    *de_dt = dot3(p->e_vec_unit,p->de_vec_dt);

    return 0;
}



void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3])
{
    r[0] = p->x;
    r[1] = p->y;
    r[2] = p->z;
    v[0] = p->vx;
    v[1] = p->vy;
    v[2] = p->vz;
}
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3])
{
    p->x = r[0];
    p->y = r[1];
    p->z = r[2];
    p->vx = v[0];
    p->vy = v[1];
    p->vz = v[2];
}
void get_e_and_h_vectors_from_particle(Particle *p, double e_vec[3], double h_vec[3])
{
    e_vec[0] = p->e_vec_x;
    e_vec[1] = p->e_vec_y;
    e_vec[2] = p->e_vec_z;
    h_vec[0] = p->h_vec_x;    
    h_vec[1] = p->h_vec_y;    
    h_vec[2] = p->h_vec_z;    
}
void set_e_and_h_vectors_in_particle(Particle *p, double e_vec[3], double h_vec[3])
{
    p->e_vec_x = e_vec[0];
    p->e_vec_y = e_vec[1];
    p->e_vec_z = e_vec[2];
    p->h_vec_x = h_vec[0];
    p->h_vec_y = h_vec[1];
    p->h_vec_z = h_vec[2];
}


/************************
/* interface parameters *
 ************************/
 
 
int set_constants(double CONST_G_, double CONST_C_, double CONST_MSUN_, double CONST_R_SUN_, double CONST_L_SUN_)
{
    CONST_G = CONST_G_;
    CONST_G_P2 = CONST_G*CONST_G;
    CONST_G_P3 = CONST_G_P2*CONST_G;
    
    CONST_C_LIGHT = CONST_C_;
    CONST_C_LIGHT_P2 = CONST_C_LIGHT*CONST_C_LIGHT;
    CONST_C_LIGHT_P4 = CONST_C_LIGHT_P2*CONST_C_LIGHT_P2;
    CONST_C_LIGHT_P5 = CONST_C_LIGHT_P4*CONST_C_LIGHT;

    CONST_MSUN = CONST_MSUN_;
    CONST_R_SUN = CONST_R_SUN_;
    CONST_L_SUN = CONST_L_SUN_;
    
    //printf("CONSTS %g %g %g\n",CONST_G,CONST_C_LIGHT,CONST_MSUN);
    
    return 0;
}

int get_relative_tolerance(double *value)
{
    *value = relative_tolerance;
    return 0;
}
int set_relative_tolerance(double value)
{
    relative_tolerance = value;
    return 0;
}
int get_absolute_tolerance_eccentricity_vectors(double *value)
{
    *value = absolute_tolerance_eccentricity_vectors;
    return 0;
}
int set_absolute_tolerance_eccentricity_vectors(double value)
{
    absolute_tolerance_eccentricity_vectors = value;
    return 0;
}

int get_include_quadrupole_order_terms(int *value){
    *value = include_quadrupole_order_terms ? 1 : 0;
    return 0;
}
int set_include_quadrupole_order_terms(int value){
    include_quadrupole_order_terms = value == 1;
    return 0;
}

int get_include_octupole_order_binary_pair_terms(int *value){
    *value = include_octupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_order_binary_pair_terms(int value){
    include_octupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_include_octupole_order_binary_triplet_terms(int *value){
    *value = include_octupole_order_binary_triplet_terms ? 1 : 0;
    return 0;
}
int set_include_octupole_order_binary_triplet_terms(int value){
    include_octupole_order_binary_triplet_terms = value == 1;
    return 0;
}

int get_include_hexadecupole_order_binary_pair_terms(int *value){
    *value = include_hexadecupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_hexadecupole_order_binary_pair_terms(int value){
    include_hexadecupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_include_dotriacontupole_order_binary_pair_terms(int *value){
    *value = include_dotriacontupole_order_binary_pair_terms ? 1 : 0;
    return 0;
}
int set_include_dotriacontupole_order_binary_pair_terms(int value){
    include_dotriacontupole_order_binary_pair_terms = value == 1;
    return 0;
}

int get_orbital_phases_random_seed(int *value)
{
    *value = orbital_phases_random_seed;
    return 0;
}
int set_orbital_phases_random_seed(int value)
{
    orbital_phases_random_seed = value;
    return 0;
}

}


