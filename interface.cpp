#include <stdio.h>
#include <math.h>
#include "interface.h"

#include "src/types.h"
//#include "interface.h"
#include "src/evolve.h"


extern "C"
{
    
/*******************
/* basic interface *
 ******************/
 
 
int add_particle(int *index, bool is_binary, bool is_external)
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

int set_integration_method(int index, int integration_method, bool KS_use_perturbing_potential)
{
    if (index > highest_particle_index)
    {
      return -1;
    }

    Particle * p = particlesMap[index];
    p->integration_method = integration_method;
    p->KS_use_perturbing_potential = KS_use_perturbing_potential;
    
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

    p->e_vec[0] = e_vec_x;
    p->e_vec[1] = e_vec_y;
    p->e_vec[2] = e_vec_z;
    p->h_vec[0] = h_vec_x;
    p->h_vec[1] = h_vec_y;
    p->h_vec[2] = h_vec_z;
    
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

    *e_vec_x = p->e_vec[0];
    *e_vec_y = p->e_vec[1];
    *e_vec_z = p->e_vec[2];
    *h_vec_x = p->h_vec[0];
    *h_vec_y = p->h_vec[1];
    *h_vec_z = p->h_vec[2];

    return 0;
}

int set_orbital_elements(int index, double semimajor_axis, double eccentricity, double true_anomaly, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node, bool sample_orbital_phase_randomly)
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
    int N_bodies, N_binaries, N_root_finding,N_ODE_equations;

    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);

    set_binary_masses_from_body_masses(&particlesMap);

    compute_orbital_vectors_from_orbital_elements(p->child1_mass, p->child2_mass, semimajor_axis, eccentricity, \
        inclination, argument_of_pericenter, longitude_of_ascending_node, \
        &(p->e_vec[0]), &(p->e_vec[1]), &(p->e_vec[2]), &(p->h_vec[0]), &(p->h_vec[1]), &(p->h_vec[2]) );
    
    p->true_anomaly = true_anomaly;
    p->sample_orbital_phase_randomly = sample_orbital_phase_randomly;
    //printf("soe a %g e %g TA %g I %g AP %g LAN %g SOPR %d\n",semimajor_axis,eccentricity,true_anomaly,inclination,argument_of_pericenter,longitude_of_ascending_node,sample_orbital_phase_randomly);
    //printf("set_orbital_elements %g %g %g %g %g %g\n",p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2]);
    
    return 0;
}
int get_orbital_elements(int index, double *semimajor_axis, double *eccentricity, double *true_anomaly, \
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
    
    *true_anomaly = p->true_anomaly;
    
    double h_tot_vec[3];
    compute_h_tot_vector(&particlesMap,h_tot_vec);

    /* determine masses in all binaries */
    int N_bodies, N_binaries, N_root_finding, N_ODE_equations;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);
    set_binary_masses_from_body_masses(&particlesMap);
    
    compute_orbital_elements_from_orbital_vectors(p->child1_mass, p->child2_mass, h_tot_vec, \
        p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2],
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

int set_instantaneous_perturbation_properties(int index, double delta_mass, double delta_X, double delta_Y, double delta_Z, double delta_VX, double delta_VY, double delta_VZ)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
    
    Particle *p = particlesMap[index];

    p->instantaneous_perturbation_delta_mass = delta_mass;
    p->instantaneous_perturbation_delta_X = delta_X;
    p->instantaneous_perturbation_delta_Y = delta_Y;
    p->instantaneous_perturbation_delta_Z = delta_Z;
    p->instantaneous_perturbation_delta_VX = delta_VX;
    p->instantaneous_perturbation_delta_VY = delta_VY;
    p->instantaneous_perturbation_delta_VZ = delta_VZ;
    
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
    
    p->external_t_ref = external_t_ref;
    p->external_e = e;
    p->external_r_p = external_r_p;
    
    /* e & h vectors for external particles are understood to be unit vectors */
    compute_orbital_vectors_from_orbital_elements_unit(INCL,AP,LAN,&(p->e_vec[0]), &(p->e_vec[1]), &(p->e_vec[2]), &(p->h_vec[0]), &(p->h_vec[1]), &(p->h_vec[2]) ); 
    
    //printf("set_external_particle_properties inputs %g %g %g %g %g\n",external_t_ref, e, external_r_p, INCL, AP, LAN);
    //printf("set_external_particle_properties OE %g %g %g %g %g %g\n",p->e_vec[0],p->e_vec[1],p->e_vec[2],p->h_vec[0],p->h_vec[1],p->h_vec[2]);
    
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
    p->spin_vec[0] = spin_vec_x;
    p->spin_vec[1] = spin_vec_y;
    p->spin_vec[2] = spin_vec_z;

    return 0;
}
int get_spin_vector(int index, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *spin_vec_x = p->spin_vec[0];
    *spin_vec_y = p->spin_vec[1];
    *spin_vec_z = p->spin_vec[2];
    
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

int get_relative_position_and_velocity(int index, double *x, double *y, double *z, double *vx, double *vy, double *vz)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];
    *x = p->r_vec[0];
    *y = p->r_vec[1];
    *z = p->r_vec[2];
    *vx = p->v_vec[0];
    *vy = p->v_vec[1];
    *vz = p->v_vec[2];
   
    return 0;
}
int get_absolute_position_and_velocity(int index, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ)
{
    if (index > highest_particle_index)
    {
      return -1;
    }
  
    Particle * p = particlesMap[index];

    *X = p->R_vec[0];
    *Y = p->R_vec[1];
    *Z = p->R_vec[2];
    *VX = p->V_vec[0];
    *VY = p->V_vec[1];
    *VZ = p->V_vec[2];
   
    return 0;
}


/************
/* PN terms *
 ************/

int set_PN_terms(int index, bool include_pairwise_1PN_terms, bool include_pairwise_25PN_terms)
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
int get_PN_terms(int index, bool *include_pairwise_1PN_terms, bool *include_pairwise_25PN_terms)
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
int set_tides_terms(int index, bool include_tidal_friction_terms, int tides_method, bool include_tidal_bulges_precession_terms, bool include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, 
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
int get_tides_terms(int index, bool *include_tidal_friction_terms, int *tides_method, bool *include_tidal_bulges_precession_terms, bool *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession, 
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
int set_root_finding_terms(int index, bool check_for_secular_breakdown, bool check_for_dynamical_instability, int dynamical_instability_criterion, int dynamical_instability_central_particle, double dynamical_instability_K_parameter,
    bool check_for_physical_collision_or_orbit_crossing, bool check_for_minimum_periapse_distance, double check_for_minimum_periapse_distance_value, bool check_for_RLOF_at_pericentre, bool check_for_RLOF_at_pericentre_use_sepinsky_fit, bool check_for_GW_condition)
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
int get_root_finding_terms(int index, bool *check_for_secular_breakdown, bool *check_for_dynamical_instability, int *dynamical_instability_criterion, int *dynamical_instability_central_particle, double *dynamical_instability_K_parameter,
    bool *check_for_physical_collision_or_orbit_crossing, bool *check_for_minimum_periapse_distance, double *check_for_minimum_periapse_distance_value, bool *check_for_RLOF_at_pericentre, bool *check_for_RLOF_at_pericentre_use_sepinsky_fit, bool *check_for_GW_condition)
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
int set_root_finding_state(int index, bool secular_breakdown_has_occurred, bool dynamical_instability_has_occurred, bool physical_collision_or_orbit_crossing_has_occurred, bool minimum_periapse_distance_has_occurred, bool RLOF_at_pericentre_has_occurred, bool GW_condition_has_occurred)
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
int get_root_finding_state(int index, bool *secular_breakdown_has_occurred, bool *dynamical_instability_has_occurred, bool *physical_collision_or_orbit_crossing_has_occurred, bool *minimum_periapse_distance_has_occurred, bool *RLOF_at_pericentre_has_occurred, bool *GW_condition_has_occurred)
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
    int N_bodies, N_binaries, N_root_finding, N_ODE_equations;
    determine_binary_parents_and_levels(&particlesMap, &N_bodies, &N_binaries, &N_root_finding,&N_ODE_equations);
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
    return 0;
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

int set_parameters(double relative_tolerance_, double absolute_tolerance_eccentricity_vectors_, 
    bool include_quadrupole_order_terms_, bool include_octupole_order_binary_pair_terms_, bool include_octupole_order_binary_triplet_terms_,
    bool include_hexadecupole_order_binary_pair_terms_, bool include_dotriacontupole_order_binary_pair_terms_, bool include_double_averaging_corrections_)
{
    relative_tolerance = relative_tolerance_;
    absolute_tolerance_eccentricity_vectors = absolute_tolerance_eccentricity_vectors_;
    include_quadrupole_order_terms = include_quadrupole_order_terms_;
    include_octupole_order_binary_pair_terms = include_octupole_order_binary_pair_terms_;
    include_octupole_order_binary_triplet_terms = include_octupole_order_binary_triplet_terms_;
    include_hexadecupole_order_binary_pair_terms = include_hexadecupole_order_binary_pair_terms_;
    include_dotriacontupole_order_binary_pair_terms = include_dotriacontupole_order_binary_pair_terms_;
    include_double_averaging_corrections = include_double_averaging_corrections_;
    //printf("PARAMS %g %g %d %d %d %d %d\n",relative_tolerance,absolute_tolerance_eccentricity_vectors,include_quadrupole_order_terms,include_octupole_order_binary_pair_terms,include_octupole_order_binary_triplet_terms,include_hexadecupole_order_binary_pair_terms,include_dotriacontupole_order_binary_pair_terms);
    
    return 0;
}

}


