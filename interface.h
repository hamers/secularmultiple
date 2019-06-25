#include "src/types.h"
extern "C"
{

/*******************
/* basic interface *
 ******************/
int add_particle(int *index, int is_binary, int is_external);
int delete_particle(int index);

int set_children(int index, int child1, int child2);
int get_children(int index, int *child1, int *child2);

int set_mass(int index, double mass);
int get_mass(int index, double *mass, double *mass_dot);

int set_mass_dot(int index, double mass_dot);
int get_mass_dot(int index, double *mass_dot);

int set_radius(int index, double radius, double radius_dot);
int get_radius(int index, double *radius, double *radius_dot);

int get_level(int index, int *level);

int set_stellar_type(int index, int value);
int get_stellar_type(int index, int *stellar_type);

int set_true_anomaly(int index, double value);
int get_true_anomaly(int index, double *value);

int set_sample_orbital_phases_randomly(int index, int value);
int get_sample_orbital_phases_randomly(int index, int *value);


/*******************************
 * instantaneous perturbations *
 * ****************************/

int set_instantaneous_perturbation_properties(int index, double delta_mass, double delta_x, double delta_y, double delta_z, double delta_vx, double delta_vy, double delta_vz);


/************
 * external *
 * *********/
 
int set_external_particle_properties(int index, double t_ref, double e, double r_p, double INCL, double AP, double LAN);


/****************
/* spin vectors *
 ****************/
int set_spin_vector(int index, double spin_vec_x, double spin_vec_y, double spin_vec_z);
int get_spin_vector(int index, double *spin_vec_x, double *spin_vec_y, double *spin_vec_z);

int set_spin_vector_dot(int index, double spin_vec_x_dot, double spin_vec_y_dot, double spin_vec_z_dot);
int get_spin_vector_dot(int index, double *spin_vec_x_dot, double *spin_vec_y_dot, double *spin_vec_z_dot);

/****************************
/* orbital vectors/elements *
 ****************************/

int set_orbital_vectors(int index, double e_vec_x, double e_vec_y, double e_vec_z, \
    double h_vec_x, double h_vec_y, double h_vec_z);
int get_orbital_vectors(int index, double *e_vec_x, double *e_vec_y, double *e_vec_z, \
    double *h_vec_x, double *h_vec_y, double *h_vec_z);
    
int set_orbital_vectors_dot(int index, double e_vec_x_dot, double e_vec_y_dot, double e_vec_z_dot, \
    double h_vec_x_dot, double h_vec_y_dot, double h_vec_z_dot);
int get_orbital_vectors_dot(int index, double *e_vec_x_dot, double *e_vec_y_dot, double *e_vec_z_dot, \
    double *h_vec_x_dot, double *h_vec_y_dot, double *h_vec_z_dot);

int set_orbital_elements(int index, double semimajor_axis, double eccentricity, double true_anomaly, \
    double inclination, double argument_of_pericenter, double longitude_of_ascending_node, int sample_orbital_phase_randomly);
int get_orbital_elements(int index, double *semimajor_axis, double *eccentricity, \
    double *inclination, double *argument_of_pericenter, double *longitude_of_ascending_node);


int set_position_vector(int index, double x, double y, double z);
int get_position_vector(int index, double *x, double *y, double *z);

int set_velocity_vector(int index, double x, double y, double z);
int get_velocity_vector(int index, double *x, double *y, double *z);

    
/************
/* PN terms *
 ************/
int set_PN_terms(int index, int include_pairwise_1PN_terms, int include_pairwise_25PN_terms);
int get_PN_terms(int index, int *include_pairwise_1PN_terms, int *include_pairwise_25PN_terms);


/*********
/* tides *
 *********/
int set_tides_terms(int index, int include_tidal_friction_terms, int tides_method, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, 
    double tides_apsidal_motion_constant, double tides_gyration_radius, double tides_viscous_time_scale, int tides_viscous_time_scale_prescription, double convective_envelope_mass, double convective_envelope_radius, double luminosity);
int get_tides_terms(int index, int *include_tidal_friction_terms, int *tides_method, int *include_tidal_bulges_precession_terms, int *include_rotation_precession_terms, double *minimum_eccentricity_for_tidal_precession, 
    double *tides_apsidal_motion_constant, double *tides_gyration_radius, double *tides_viscous_time_scale, int *tides_viscous_time_scale_prescription, double *convective_envelope_mass, double *convective_envelope_radius, double *luminosity);

/*******
 * VRR *
 * *****/
 
int set_VRR_properties(int index, int VRR_model, int VRR_include_mass_precession, double VRR_mass_precession_rate, 
    double VRR_Omega_vec_x, double VRR_Omega_vec_y, double VRR_Omega_vec_z, 
    double VRR_eta_20_init, double VRR_eta_a_22_init, double VRR_eta_b_22_init, double VRR_eta_a_21_init, double VRR_eta_b_21_init,
    double VRR_eta_20_final, double VRR_eta_a_22_final, double VRR_eta_b_22_final, double VRR_eta_a_21_final, double VRR_eta_b_21_final,
	double VRR_initial_time, double VRR_final_time);

/****************
/* root finding *
 ****************/

int set_root_finding_terms(int index, int check_for_secular_breakdown, int check_for_dynamical_instability, int dynamical_instability_criterion, int dynamical_instability_central_particle, int dynamical_instability_K_parameter,
    int check_for_physical_collision_or_orbit_crossing, int check_for_minimum_periapse_distance, double check_for_minimum_periapse_distance_value, int check_for_RLOF_at_pericentre, int check_for_RLOF_at_pericentre_use_sepinsky_fit);
int get_root_finding_terms(int index, int *check_for_secular_breakdown, int *check_for_dynamical_instability, int *dynamical_instability_criterion, int *dynamical_instability_central_particle, int *dynamical_instability_K_parameter,
    int *check_for_physical_collision_or_orbit_crossing, int *check_for_minimum_periapse_distance, double *check_for_minimum_periapse_distance_value, int *check_for_RLOF_at_pericentre, int *check_for_RLOF_at_pericentre_use_sepinsky_fit);
int set_root_finding_state(int index, int secular_breakdown_has_occurred, int dynamical_instability_has_occurred, int physical_collision_or_orbit_crossing_has_occurred, int minimum_periapse_distance_has_occurred, int RLOF_at_pericentre_has_occurred);
int get_root_finding_state(int index, int *secular_breakdown_has_occurred, int *dynamical_instability_has_occurred, int *physical_collision_or_orbit_crossing_has_occurred, int* minimum_periapse_distance_has_occurred, int *RLOF_at_pericentre_has_occurred);


/***********************
/* interface functions *
 ***********************/
int evolve_interface(double start_time, double time_step, double *output_time, double *hamiltonian, int *flag, int *error_code);
int determine_binary_parents_levels_and_masses_interface();
int apply_external_perturbation_assuming_integrated_orbits_interface();
int apply_user_specified_instantaneous_perturbation_interface();
int set_positions_and_velocities_interface();

/**********************************************
/* orbital element/vector conversion routines *
 **********************************************/
int compute_h_tot_vector(ParticlesMap* particlesMap, double h_tot_vector[3]);
int compute_orbital_vectors_from_orbital_elements(double child1_mass, double child2_mass, double semimajor_axis, double eccentricity, double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_vec_x, double *e_vec_y, double *e_vec_z, double *h_vec_x, double *h_vec_y, double *h_vec_z);
int compute_orbital_vectors_from_orbital_elements_unit(double inclination, double argument_of_pericenter,double longitude_of_ascending_node, double *e_hat_vec_x, double *e_hat_vec_y, double *e_hat_vec_z, double *h_hat_vec_x, double *h_hat_vec_y, double *h_hat_vec_z);
int compute_orbital_elements_from_orbital_vectors(double child1_mass, double child2_mass, double h_tot_vec[3], double e_vec_x, double e_vec_y, double e_vec_z, double h_vec_x, double h_vec_y, double h_vec_z, double *semimajor_axis, double *eccentricity, double *inclination, double *argument_of_pericenter,double *longitude_of_ascending_node);
int get_inclination_relative_to_parent(int index, double *inclination_relative_to_parent);

void compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity, double *cos_eccentric_anomaly, double *sin_eccentric_anomaly);
void compute_true_anomaly_from_eccentric_anomaly(double cos_eccentric_anomaly, double sin_eccentric_anomaly, double eccentricity, double *cos_true_anomaly, double *sin_true_anomaly);
double compute_true_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity);
double sample_random_true_anomaly(double eccentricity,int seed);

void from_orbital_vectors_to_cartesian(double child1_mass, double child2_mass, double e_vec[3], double h_vec[3], double true_anomaly, double r[3], double v[3]);
void from_cartesian_to_orbital_vectors(double child1_mass, double child2_mass, double r[3], double v[3], double e_vec[3], double h_vec[3]);

int get_de_dt(int index, double *de_dt);

void get_position_and_velocity_vectors_from_particle(Particle *p, double r[3], double v[3]);
void set_position_and_velocity_vectors_in_particle(Particle *p,  double r[3], double v[3]);
void get_e_and_h_vectors_from_particle(Particle *p, double e_vec[3], double h_vec[3]);
void set_e_and_h_vectors_in_particle(Particle *p, double e_vec[3], double h_vec[3]);

/************************
/* interface parameters *
 ************************/
extern double relative_tolerance;
extern double absolute_tolerance_eccentricity_vectors;
extern bool include_quadrupole_order_terms;
extern bool include_octupole_order_binary_pair_terms;
extern bool include_octupole_order_binary_triplet_terms;
extern bool include_hexadecupole_order_binary_pair_terms;
extern bool include_dotriacontupole_order_binary_pair_terms;
extern int orbital_phases_random_seed;

int set_constants(double CONST_G_, double CONST_C_, double CONST_MSUN_, double CONST_R_SUN_, double CONST_L_SUN_);

int get_relative_tolerance(double *value);
int set_relative_tolerance(double value);

int get_absolute_tolerance_eccentricity_vectors(double *value);
int set_absolute_tolerance_eccentricity_vectors(double value);

int get_include_quadrupole_order_terms(int *value);
int set_include_quadrupole_order_terms(int value);

int get_include_octupole_order_binary_pair_terms(int *value);
int set_include_octupole_order_binary_pair_terms(int value);

int get_include_octupole_order_binary_triplet_terms(int *value);
int set_include_octupole_order_binary_triplet_terms(int value);

int get_include_hexadecupole_order_binary_pair_terms(int *value);
int set_include_hexadecupole_order_binary_pair_terms(int value);

int get_include_dotriacontupole_order_binary_pair_terms(int *value);
int set_include_dotriacontupole_order_binary_pair_terms(int value);

int get_orbital_phases_random_seed(int *value);
int set_orbital_phases_random_seed(int value);

}
