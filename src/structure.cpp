/*
*/

#include "structure.h"
#include "evolve.h"
#include <stdio.h>

extern "C"
{

int determine_binary_parents_and_levels(ParticlesMap *particlesMap, int *N_bodies, int *N_binaries, int *N_root_finding, int *N_ODE_equations)
{
    *N_bodies = 0;
    *N_binaries = 0;
    *N_root_finding = 0;
    *N_ODE_equations = 0;
    
    /* determine parent for each particle */
    ParticlesMapIterator it_p,it_q;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        P_p->parent = -1;
    }

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        if (P_p->is_binary == true)
        {
            (*N_binaries)++;
            
            /* root finding */
            if (P_p->check_for_secular_breakdown == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_dynamical_instability == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_minimum_periapse_distance == true)
            {
                (*N_root_finding)++;
            }
            if (P_p->check_for_GW_condition == true)
            {
                (*N_root_finding)++;
            }

            /* parents and siblings */
            for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
            {
                Particle *P_q = (*it_q).second;
                if (P_q->index == P_p->child1)
                {
                    P_q->parent = P_p->index;
                    P_q->sibling = P_p->child2;
                }
                if (P_q->index == P_p->child2)
                {
                    
                    P_q->parent = P_p->index;
                    P_q->sibling = P_p->child1;
                }
            }
            
            /* count number of ODE eqs */
            if (P_p->integration_method==0)
            {
                (*N_ODE_equations) += 6;
            }
            else if (P_p->integration_method==1)
            {
                (*N_ODE_equations) += 10;
            }
            else if (P_p->integration_method==2)
            {
                (*N_ODE_equations) += 6;
            }

        }
        else
        {
            (*N_bodies)++;
            (*N_ODE_equations) += 5;
            
            /* root finding */
            if (P_p->check_for_RLOF_at_pericentre == true)
            {
                (*N_root_finding)++;
            }
            
        }
    }

    /* determine levels and set of parents for each particle */
    int highest_level = 0;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;

        P_p->connecting_child_in_parents.clear();
        P_p->parents.clear();
        P_p->level=0;
        
        int child = P_p->index;
        int parent = P_p->parent;

        if (parent != -1) /* if parent == -1, P_p is the `top' binary, for which level=0 */
        {
            while (parent != -1) /* search parents until reaching the top binary */
            {
                for (it_q = particlesMap->begin(); it_q != particlesMap->end(); it_q++)
                {
                    Particle *P_q = (*it_q).second;
                    if (P_q->index == parent)
                    {
                        if (child==P_q->child1)
                        {
                            P_p->connecting_child_in_parents.push_back(1);
                        }
                        else if (child==P_q->child2)
                        {
                            P_p->connecting_child_in_parents.push_back(2);
                        }
                        P_p->parents.push_back(parent);
                        P_p->level++;
                        
                        child = P_q->index;
                        parent = P_q->parent;
                        
                        #ifdef DEBUG
                        printf("structure.cpp -- determine_binary_parents_and_levels -- p %d q %d %d child %d\n",P_p->index,P_q->index,P_p->level,child);
                        #endif
                        break;
                    }
                }
            }
        }
        highest_level = max(P_p->level,highest_level);
        
    }

    /* write highest level to all particles -- needed for function set_binary_masses_from_body_masses */
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        P_p->highest_level = highest_level;
    }

    return 0;
}

void set_binary_masses_from_body_masses(ParticlesMap *particlesMap)
{

    /* set binary masses -- to ensure this happens correctly, do this from highest level to lowest level */
    ParticlesMapIterator it_p;

    int highest_level = particlesMap->begin()->second->highest_level;

    int level=highest_level;

    while (level > -1)
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *P_p = (*it_p).second;
            if ((P_p->is_binary == true) && (P_p->level == level))
            {
                Particle *P_p_child1 = (*particlesMap)[P_p->child1];
                Particle *P_p_child2 = (*particlesMap)[P_p->child2];
                
                /* these quantities are used in ODE_system.cpp */
                P_p->child1_mass = P_p_child1->mass;
                P_p->child2_mass = P_p_child2->mass;
                P_p->mass = P_p->child1_mass + P_p->child2_mass;

                P_p->child1_mass_dot = P_p_child1->mass_dot;
                P_p->child2_mass_dot = P_p_child2->mass_dot;
                P_p->mass_dot = P_p->child1_mass_dot + P_p->child2_mass_dot;

                P_p->child1_mass_plus_child2_mass = P_p->child1_mass + P_p->child2_mass;
                P_p->child1_mass_minus_child2_mass = P_p->child1_mass - P_p->child2_mass;
                P_p->child1_mass_times_child2_mass = P_p->child1_mass*P_p->child2_mass;
                
                P_p->mu = P_p->child1_mass * P_p->child2_mass / P_p->mass;
                
                #ifdef DEBUG
                printf("structure.cpp -- set_binary_masses_from_body_masses -- level %d m %g highest_level %d\n",level,P_p->mass,highest_level);
                #endif
            }
        }
        level--;
    }

    /* determine total system mass -- needed for hyperbolic external orbits */
    double total_system_mass;
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->level==0) /* lowest-level binary */
        {
            total_system_mass = P_p->child1_mass + P_p->child2_mass;
            break;
        }
    }

    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        
        P_p->total_system_mass = total_system_mass;
    }
}

void set_positions_and_velocities(ParticlesMap *particlesMap) /* TO DO: add to name of function: "_of_all_bodies" */
{
    /* Compute and set the positions and velocities of all bodies */
    /* By default, sample orbital phases randomly 
     * if particle.sample_orbital_phases_randomly == False: look for particle.true_anomaly */

    set_binary_masses_from_body_masses(particlesMap);

    double true_anomaly;

    int index = 0;
    int i;
    
    double r[3],v[3],r_parent[3],v_parent[3],r_child1[3],v_child1[3],r_child2[3],v_child2[3];
    double parent_mass,child1_mass,child2_mass;
    double e;
    
    /*/ Go from the top of the system (level=0) downwards */
    
    ParticlesMapIterator it;
    int highest_level = particlesMap->begin()->second->highest_level;
    int level = 0;
    while (level < highest_level)
    {
        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *parent = (*it).second;
            if ((parent->is_binary == true) && (parent->level == level))
            {
                Particle *child1 = (*particlesMap)[parent->child1];
                Particle *child2 = (*particlesMap)[parent->child2];
                
                child1_mass = child1->mass;
                child2_mass = child2->mass;
                parent_mass = parent->mass;
                
                if (parent->sample_orbital_phase_randomly == false)
                {
                    true_anomaly = parent->true_anomaly;
                }
                else
                {
                    e = norm3(parent->e_vec);
                    true_anomaly = sample_random_true_anomaly(e);
                    //printf("parent->sample_orbital_phases_randomly  %d %g \n",parent->sample_orbital_phases_randomly,true_anomaly);
                }
                
                from_orbital_vectors_to_cartesian(
                    child1_mass,child2_mass,
                    parent->e_vec,parent->h_vec,
                    true_anomaly,
                    r,v);
                //printf("fov %g %g %g %g %g %g\n",e_vec[0],e_vec[1],e_vec[2],r[0],r[1],r[2]);
                
                if (parent->level == 0)
                {
                     /* without loss of generality, set the initial CM of the system to the origin */
                    for (i=0; i<3; i++)
                    {
                        r_parent[i] = 0.0;
                        v_parent[i] = 0.0;
                    }
                }
                else
                {
                    get_position_and_velocity_vectors_from_particle(parent,r_parent,v_parent);
                }
                
                for (i=0; i<3; i++)
                {
                    r_child1[i] = r_parent[i] + (child2_mass/parent_mass)*r[i];
                    v_child1[i] = v_parent[i] + (child2_mass/parent_mass)*v[i];
                    
                    r_child2[i] = r_parent[i] - (child1_mass/parent_mass)*r[i];
                    v_child2[i] = v_parent[i] - (child1_mass/parent_mass)*v[i];
                }
                set_position_and_velocity_vectors_in_particle(child1,r_child1,v_child1);
                set_position_and_velocity_vectors_in_particle(child2,r_child2,v_child2);

                parent->true_anomaly = true_anomaly;
            }
        }
        level++;
    }
}

void update_masses_positions_and_velocities_of_all_binaries(ParticlesMap *particlesMap)
{
    set_binary_masses_from_body_masses(particlesMap);

    /* set binary positions and velocities -- to ensure this happens correctly, do this from highest level to lowest level */
 
    int i;
    double child1_mass,child2_mass;
    double r[3],v[3];
    double r_child1[3],v_child1[3];
    double r_child2[3],v_child2[3];
    
    ParticlesMapIterator it_p;
    int highest_level = particlesMap->begin()->second->highest_level;
    int level=highest_level;

    while (level > -1)
    {
        for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
        {
            Particle *p = (*it_p).second;
            if ((p->is_binary == true) && (p->level == level))
            {
                Particle *child1 = (*particlesMap)[p->child1];
                Particle *child2 = (*particlesMap)[p->child2];
                
                
                get_position_and_velocity_vectors_from_particle(child1,r_child1,v_child1);
                get_position_and_velocity_vectors_from_particle(child2,r_child2,v_child2);
                
                //printf("level %d %d %d\n",p->level,child1->is_binary,child2->is_binary);

                child1_mass = child1->mass;
                child2_mass = child2->mass;

                for (i=0; i<3; i++)
                {
                    r[i] = (r_child1[i]*child1_mass + r_child2[i]*child2_mass)/(child1_mass+child2_mass);
                    v[i] = (v_child1[i]*child1_mass + v_child2[i]*child2_mass)/(child1_mass+child2_mass);

                    //printf("r_child1[i] %g r[i] %g\n",r_child1[i],r[i]);
                    //printf("r_child2[i] %g r[i] %g\n",r_child2[i],r[i]);

                }
            
                set_position_and_velocity_vectors_in_particle(p,r,v);
//                printf("level %d m %g hl %d\n",level,P_p->mass,highest_level);
            }
        }
        level--;
    }
}


void update_orbital_vectors_in_binaries_from_positions_and_velocities(ParticlesMap *particlesMap)
{
    //printf("update_orbital_vectors_in_binaries_from_positions_and_velocities\n");

    int index = 0;
    int i;
    
    double r[3],v[3],r_child1[3],v_child1[3],r_child2[3],v_child2[3];
    double child1_mass,child2_mass;
    double e_vec[3],h_vec[3];
    
    /*/ Go from the top of the system (level=0) downwards */
    
    ParticlesMapIterator it;

    int highest_level = particlesMap->begin()->second->highest_level;
    int level = 0;
    while (level < highest_level)
    {
        for (it = particlesMap->begin(); it != particlesMap->end(); it++)
        {
            Particle *parent = (*it).second;
            if ((parent->is_binary == true) && (parent->level == level))
            {
                Particle *child1 = (*particlesMap)[parent->child1];
                Particle *child2 = (*particlesMap)[parent->child2];
                
                child1_mass = child1->mass;
                child2_mass = child2->mass;
                
                //get_position_and_velocity_vectors_from_particle(parent,r,v);
                get_position_and_velocity_vectors_from_particle(child1,r_child1,v_child1);
                get_position_and_velocity_vectors_from_particle(child2,r_child2,v_child2);
                //printf("level %d\n",parent->level);
                for (i=0; i<3; i++)
                {
                    r[i] = r_child1[i] - r_child2[i];
                    v[i] = v_child1[i] - v_child2[i];
                    //printf("r_child1[i] %g r[i] %g %g\n",r_child1[i],r[i],parent->x);
                    //printf("r_child2[i] %g r[i] %g\n",r_child2[i],r[i]);
                }
                
                from_cartesian_to_orbital_vectors(
                    child1_mass,child2_mass,
                    r,v,
                    parent->e_vec,parent->h_vec,&parent->true_anomaly);
            }
        }
        level++;
    }
}

}
