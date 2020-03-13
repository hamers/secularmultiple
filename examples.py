"""
Some examples illustrating the usage of SecularMultiple
Adrian Hamers, March 2020
"""

import numpy as np
import numpy.random as randomf

import argparse
import time

from secularmultiple import SecularMultiple,Particle,Tools

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def add_bool_arg(parser, name, default=False,help=None):
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--' + name, dest=name, action='store_true',help="Enable %s"%help)
    group.add_argument('--no-' + name, dest=name, action='store_false',help="Disable %s"%help)
    parser.set_defaults(**{name:default})

def parse_arguments():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--e",                           type=int,     dest="example",                        default=1,              help="Example number")
    
    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                         default=True,        help="Verbose terminal output")
    add_bool_arg(parser, 'plot',                            default=True,         help="Make plots")
    add_bool_arg(parser, 'show',                            default=True,         help="Show plots")
    
    args = parser.parse_args()

    return args
    
class examples():
    
    def example1(self,args):
        """
        Example of a three-planet system with tides in the innermost planet
        System parameters taken from http://adsabs.harvard.edu/abs/2011ApJ...735..109W
        Units used in SecularMultiple: 
        length -- AU
        time -- yr
        mass -- MSun
        
        """

        code = SecularMultiple() ### initialize the code
        CONST_G = code.CONST_G ### extract physical constants from the code
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        RJup = 0.1027922358015816*CONST_R_SUN
        MJup = 0.0009546386983890755
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        N=4
        m0 = 1.0
        m1 = 0.5*MJup
        m2 = MJup
        m3 = 1.5*MJup
        
        a1 = 1.0
        a2 = 6.0
        a3 = 16.0
        e1 = 0.066
        e2 = 0.188
        e3 = 0.334
        i1 = 4.5*np.pi/180.0
        i2 = 19.9*np.pi/180.0
        i3 = 7.9*np.pi/180.0
        AP1 = np.pi
        AP2 = 0.38*np.pi
        AP3 = np.pi
        LAN1 = 0.01
        LAN2 = np.pi
        LAN3 = 0.01

        R0 = 1.0*CONST_R_SUN
        R1 = 1.0*RJup
        R2 = 1.0*RJup
        R3 = 1.0*RJup

        masses = [m0,m1,m2,m3]
        radii = [R0,R1,R2,R3]
        semimajor_axes = [a1,a2,a3]
        eccentricities = [e1,e2,e3]
        inclinations = [i1,i2,i3]
        APs = [AP1,AP2,AP3]
        LANs = [LAN1,LAN2,LAN3]
        
        particles = Tools.create_nested_multiple(N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii=radii) 
        orbits = [x for x in particles if x.is_binary==True]
        N_orbits = len(orbits)
        
        particles[1].spin_vec_x = 0.0
        particles[1].spin_vec_y = 0.0
        particles[1].spin_vec_z = 4.0e-2/day

        k_L = 0.38
        k_AM = k_L/2.0
        rg = 0.25
        tau = 0.66*second
        T = R1**3/(CONST_G*m1*tau)
        t_V = 3.0*(1.0 + 1.0/k_L)*T

        particles[0].include_tidal_friction_terms = False
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = True
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-8

        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_viscous_time_scale = t_V
        particles[1].tides_gyration_radius = rg

        orbits[0].include_1PN_terms = False ### do not include 1PN terms here
        code.add_particles(particles)
        primary = code.particles[0]

        code.enable_tides = True
        code.enable_root_finding = True
        
        a_AU_print = [[] for x in range(N_orbits)]
        e_print = [[] for x in range(N_orbits)]
        INCL_print = [[] for x in range(N_orbits)]
        rel_INCL_print = [[] for x in range(N_orbits)]
        t_print = []
        
        t = 0.0
        Nsteps = 2000
        tend = 3.0e8
        dt = tend/float(Nsteps)
        import time
        start = time.time()
        while t<=tend:

            code.evolve_model(t)
            t+=dt
        
            if args.verbose==True:
                print('t/Myr',t*1e-6,'es',[o.e for o in orbits])
            for i in range(N_orbits):
                rel_INCL_print[i].append(orbits[i].INCL_parent)
                a_AU_print[i].append(orbits[i].a)
                e_print[i].append(orbits[i].e)
                INCL_print[i].append(orbits[i].INCL)
            t_print.append(t)
            
        if args.verbose==True:
            print('wall time',time.time()-start)
        
        t_print = np.array(t_print)
        for i in range(N_orbits):
            INCL_print[i] = np.array(INCL_print[i])
            rel_INCL_print[i] = np.array(rel_INCL_print[i])
            e_print[i] = np.array(e_print[i])
            a_AU_print[i] = np.array(a_AU_print[i])
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            fig=pyplot.figure(figsize=(8,8))
            plot1=fig.add_subplot(2,1,1,yscale="log")
            plot2=fig.add_subplot(2,1,2,yscale="linear")
            colors=['k','r','g']
            for i in range(N_orbits):
                color=colors[i]
                plot1.plot(1.0e-6*t_print,a_AU_print[i],color=color)
                plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0-e_print[i]),color=color)
                plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0+e_print[i]),color=color)
                plot2.plot(1.0e-6*t_print,INCL_print[i]*180.0/np.pi,color=color)
                
                plot1.set_xlabel("$t/\mathrm{Myr}$")
                plot2.set_xlabel("$t/\mathrm{Myr}$")
                plot1.set_ylabel("$r_i/\mathrm{AU}$")
                plot2.set_ylabel("$\mathrm{incl}_i/\mathrm{deg}$")
            fig.savefig("example1.pdf")
            if args.show==True:
                pyplot.show()

    def example2(self,args):
        """
        Lidov-Kozai problem of a planet around a star around a supermassive black hole.
        Includes perturbations from other stars in the form of vector resonant relaxation (VRR)
        
        """

        code = SecularMultiple() ### initialize the code
        CONST_G = code.CONST_G ### extract physical constants from the code
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        RJup = 0.1027922358015816*CONST_R_SUN
        MJup = 0.0009546386983890755
        day = 1.0/365.25
        second = day/(24.0*3600.0)
        meter = 1.0/1.496e+11

        ### Input parameters ###
        m1 = 1.0
        m2 = MJup
        m3 = 4.0e6
        
        a1 = 1.0e-1
        a2 = 1.0e4
        e1 = 0.01
        e2 = 0.1
        i1 = 4.5*np.pi/180.0
        i2 = 19.9*np.pi/180.0
        AP1 = np.pi
        AP2 = 0.38*np.pi
        LAN1 = 0.01
        LAN2 = np.pi

        R1 = 1.0*CONST_R_SUN
        R2 = 1.0*RJup
        R3 = CONST_G*m3/(CONST_C**2)

        m_star = 1.0
        gamma = 3.0/2.0
        VRR_model = 3

        ### Simulation parameters ###
        VRR_include_mass_precession = True
        include_inner_1PN_terms = True
        include_outer_1PN_terms = True
        
        ### Process parameters ###
        P1 = 2.0*np.pi*np.sqrt(a1**3/(CONST_G*(m1+m2)))
        P2 = 2.0*np.pi*np.sqrt(a2**3/(CONST_G*(m1+m2+m3)))

        masses = [m1,m2,m3]
        radii = [R1,R2,R3]
        semimajor_axes = [a1,a2]
        eccentricities = [e1,e2]
        inclinations = [i1,i2]
        APs = [AP1,AP2]
        LANs = [LAN1,LAN2]

        N = len(masses)
        particles = Tools.create_nested_multiple(N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii=radii) 
        orbits = [x for x in particles if x.is_binary==True]
        N_orbits = len(orbits)

        inner_orbit = orbits[0]
        outer_orbit = orbits[1]
        
        c1 = 4.8
        c2 = -2.9
        log10_sigma_h_km_s = (np.log10(m3) - c2)/c1
        
        sigma_h_km_s = pow(10.0,log10_sigma_h_km_s)
        sigma_h = 1.0e3*sigma_h_km_s*meter/second
        
        if args.verbose==True:
            print('sigma_h_km_s',sigma_h_km_s,'sigma_h',sigma_h)

#       K_12 = K_12_function(gamma)
#        K_32 = K_32_function(gamma)
#        C_NRR = ((3.0*numpy.pi)/(64.0))*1.0/( K_12 - (1.0/5.0)*K_32 + (5.0*numpy.pi/8.0)*(1.0/(2.0*gamma-1.0)) )
            
        r_h = CONST_G*m3*(1.0/(sigma_h**2*(1.0+gamma)))*(1.0 + (1.0 + gamma)/(gamma - 1.0))
        
        r_0 = r_h
        n_0 = (2.0*m3/m_star)*((3.0-gamma)/(4.0*np.pi*r_h**3))

        r = a2
        rho_star = compute_rho_star_r(r,gamma,n_0,r_0,m_star)
        n_star = compute_n_star_r(r,gamma,n_0,r_0,m_star)
        M_star = compute_M_star_r(r,gamma,n_0,r_0,m_star)
        N_star = compute_N_star_r(r,gamma,n_0,r_0,m_star)
        sigma_r = compute_sigma_r(r,gamma,n_0,r_0,m_star,m3,CONST_G)
        
        if args.verbose==True:
            print('n_star',n_star,'M_star',M_star,'N_star',N_star,'sigma_r',sigma_r)
        
        LK_timescale = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
        if args.verbose==True:
            print('LK_timescale',LK_timescale)
        VRR_mass_precession_timescale = (1.0/2.0)*pow(1.0-e2**2,-1.0/2.0)*(m3/M_star)*P2
        VRR_mass_precession_rate = 1.0/VRR_mass_precession_timescale
        VRR_timescale = (P2/2.0)*(m3/m_star)*1.0/np.sqrt(N_star)
        #VRR_timescale *= 0.1
                
        if args.verbose==True:
            print('VRR_mass_precession_timescale',VRR_mass_precession_timescale,'VRR_timescale',VRR_timescale)
       
        outer_orbit.VRR_include_mass_precession = VRR_include_mass_precession
        outer_orbit.VRR_mass_precession_rate = VRR_mass_precession_rate

        VRR_reorientation_timestep = np.sqrt(0.1)*VRR_timescale
        if args.verbose==True:
            print('VRR_reorientation_timestep',VRR_reorientation_timestep)

        outer_orbit.VRR_model = VRR_model
        reorientation_function(VRR_model,VRR_timescale,VRR_reorientation_timestep,outer_orbit)
        

        v_bin = np.sqrt(CONST_G*(m1+m2)/a1)
        q_sigma = (m1+m2)/m_star
        log_Lambda = np.log( 3.0*((1.0 + 1.0/q_sigma)/(1.0 + 2.0/q_sigma))*sigma_r**2/v_bin**2 )
        evaporation_timescale = np.sqrt( (1.0+q_sigma)/(2.0*np.pi*q_sigma) )*(m1+m2)*sigma_r/(8.0*np.sqrt(np.pi)*CONST_G*a1*m_star**2*n_star*log_Lambda)
        if args.verbose==True:
            print('evaporation_timescale',evaporation_timescale)
        

        inner_orbit.include_1PN_terms = include_inner_1PN_terms
        outer_orbit.include_1PN_terms = include_outer_1PN_terms
        code.add_particles(particles)
        primary = code.particles[0]

        code.enable_tides = False
        code.enable_root_finding = True
        code.enable_VRR = True
        
        
        a_AU_print = [[] for x in range(N_orbits)]
        e_print = [[] for x in range(N_orbits)]
        INCL_print = [[] for x in range(N_orbits)]
        rel_INCL_print = [[] for x in range(N_orbits)]
        t_print = []
        
        t = 0.0
        Nsteps = 1000
        tend = evaporation_timescale
        dt_fixed = tend/float(Nsteps)
        t_next_reorientation = VRR_reorientation_timestep
        
        import time
        start = time.time()
        while t<=tend:
            dt = dt_fixed
            if t+dt > t_next_reorientation:
                dt = t_next_reorientation - t
                t_next_reorientation += VRR_reorientation_timestep

                reorientation_function(VRR_model,VRR_timescale,t_next_reorientation,outer_orbit)

            t+=dt
            code.evolve_model(t)

            if args.verbose==True:
                print('t',t,'es',[o.e for o in orbits],'Omegas',[o.LAN for o in orbits])
            for i in range(N_orbits):
                rel_INCL_print[i].append(orbits[i].INCL_parent)
                a_AU_print[i].append(orbits[i].a)
                e_print[i].append(orbits[i].e)
                INCL_print[i].append(orbits[i].INCL)
            t_print.append(t)
            
        if args.verbose==True:
            print('wall time',time.time()-start)
        
        t_print = np.array(t_print)
        for i in range(N_orbits):
            INCL_print[i] = np.array(INCL_print[i])
            rel_INCL_print[i] = np.array(rel_INCL_print[i])
            e_print[i] = np.array(e_print[i])
            a_AU_print[i] = np.array(a_AU_print[i])
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            fig=pyplot.figure(figsize=(8,8))
            plot1=fig.add_subplot(2,1,1,yscale="log")
            plot2=fig.add_subplot(2,1,2,yscale="linear")
            colors=['k','r','g']
            for i in range(N_orbits):
                color=colors[i]
                plot1.plot(1.0e-6*t_print,a_AU_print[i],color=color)
                plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0-e_print[i]),color=color)
                plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0+e_print[i]),color=color)
                plot2.plot(1.0e-6*t_print,rel_INCL_print[i]*180.0/np.pi,color=color)
                
                plot1.set_xlabel("$t/\mathrm{Myr}$")
                plot2.set_xlabel("$t/\mathrm{Myr}$")
                plot1.set_ylabel("$r_i/\mathrm{AU}$")
                plot2.set_ylabel("$\mathrm{incl}_\mathrm{rel}/\mathrm{deg}$")
            fig.savefig("example2.pdf")
            pyplot.show()

    def example3(self,args):
        m1=1.0
        m2=1.0e-6
        m3=1.0
        e1=0
        e2=0.4
        a1=1.0
        a2=10.0

        i1=0.2
        i2=65.0*np.pi/180.0
        AP1=0
        AP2=0
        LAN1=0
        LAN2=0

        do_nbody=True
        particles = Tools.create_nested_multiple(3, [m1,m2,m3],[a1,a2],[e1,e2],[i1,i2],[AP1,AP2],[LAN1,LAN2])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]
        N_binaries = len(binaries)
        N_bodies = len(bodies)
        
        code = SecularMultiple()
        code.add_particles(particles)
        
        CONST_G = code.CONST_G
        P1=2.0*np.pi*np.sqrt(a1**3/(CONST_G*(m1+m2)))
        P2=2.0*np.pi*np.sqrt(a2**3/(CONST_G*(m1+m2+m3)))
        P_LK12 = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
        if args.verbose==True:
            print("Ps",P1*1e-6,P2*1e-6)
            print("P_LKs",P_LK12*1e-6)

        N = 5000
        tend = 1.e4

        integration_methods = [[0,0],[0,1]]
        #integration_methods = [0,0,0]
        KS_use_V = [[True,True],[True,True]]
        #KS_use_V = [True,True,False]
        terms = [[False,True,True,True,True,True],[False,True,True,True,True,True]]
       
        import time
        
        data_arrays = []
        
        for index_combination,integration_method in enumerate(integration_methods):
            if args.verbose==True:
                print("index_combination",index_combination)

            particles = Tools.create_nested_multiple(3, [m1,m2,m3],[a1,a2],[e1,e2],[i1,i2],[AP1,AP2],[LAN1,LAN2])
            bodies = [x for x in particles if x.is_binary==False]
            binaries = [x for x in particles if x.is_binary==True]

            binaries[0].integration_method = integration_methods[index_combination][0]
            binaries[1].integration_method = integration_methods[index_combination][1]

            binaries[0].KS_use_perturbing_potential = KS_use_V[index_combination][0]
            binaries[1].KS_use_perturbing_potential = KS_use_V[index_combination][1]

            code = SecularMultiple()
            code.add_particles(particles)

            code.enable_root_finding = True
            binaries[0].check_for_physical_collision_or_orbit_crossing=True
            bodies[0].radius=1.0e-5
            bodies[1].radius=1.0e-5
            
            code.include_double_averaging_corrections = terms[index_combination][0]
            code.include_quadrupole_order_terms = terms[index_combination][1]
            code.include_octupole_order_binary_pair_terms = terms[index_combination][2]
            code.include_octupole_order_binary_triplet_terms = terms[index_combination][3]
            code.include_hexadecupole_order_binary_pair_terms = terms[index_combination][4]
            code.include_dotriacontupole_order_binary_pair_terms = terms[index_combination][5]
            
            if args.verbose==True:
                print("Integration methods ",[x.integration_method for x in binaries],"KS_V",[x.KS_use_perturbing_potential for x in binaries],"terms",code.include_double_averaging_corrections,code.include_quadrupole_order_terms,code.include_octupole_order_binary_pair_terms,code.include_octupole_order_binary_triplet_terms,code.include_hexadecupole_order_binary_pair_terms,code.include_dotriacontupole_order_binary_pair_terms)
        
            a_print = [[] for x in range(N_binaries)]
            e_print = [[] for x in range(N_binaries)]
            i_print = [[] for x in range(N_binaries)]
            rel_INCL_print = [[] for x in range(N_binaries)]
            t_print = []

            start = time.time()
            t = 0.0
            dt = tend/float(N)
            while t<tend:
                t+=dt
                code.evolve_model(t)
                
                if args.verbose==True:
                    print('t',t,'es',[o.e for o in binaries])

                for i in range(N_binaries):
                    a_print[i].append([binaries[i].a])
                    e_print[i].append([binaries[i].e])
                    i_print[i].append(binaries[i].INCL)
                    rel_INCL_print[i].append(binaries[i].INCL_parent)
                t_print.append(t)
                
            wall_time = time.time()-start
            code.reset()
            
            t_print = np.array(t_print)
            for i in range(N_binaries):
                a_print[i] = np.array(a_print[i])
                e_print[i] = np.array(e_print[i])
                i_print[i] = np.array(i_print[i])

            data = {'t_print':t_print,'a_print':a_print,'e_print':e_print,'i_print':i_print,'wall_time':wall_time,'integration_methods':integration_methods[index_combination],'KS_use_perturbing_potential':KS_use_V[index_combination],'terms':terms[index_combination]}
            data_arrays.append(data)


        if HAS_MATPLOTLIB==True and args.plot==True:
            linestyles=['solid','dotted','dashed','-.']
            linewidth=2.0

            fig=pyplot.figure(figsize=(8,8))
            plot1=fig.add_subplot(2,1,1,yscale="linear")
            plot2=fig.add_subplot(2,1,2,yscale="log")

            linewidths=[1.5,2.5,1.5]
            colors =['k','tab:red','tab:orange']

            for index_combination,data in enumerate(data_arrays):
                linewidth=linewidths[index_combination]
                linestyle=linestyles[index_combination]

                N_binaries = len(data["a_print"])
                color=colors[index_combination]
                
                if index_combination==0:
                    label = "$\mathrm{Double\,averaged; \,WT=%s\,s}$"%round(data["wall_time"],1)
                elif index_combination==1:
                    label = "$\mathrm{Single\,averaged; \,WT=%s\,s}$"%round(data["wall_time"],1)
                
                for i in range(N_binaries):
                    if i!=0:
                        label=""
                        label_nb=""

                    plot1.plot(1.0e-6*data["t_print"],data["i_print"][i]*180.0/np.pi,color=color,linestyle=linestyle,linewidth=linewidth)
                    plot2.plot(1.0e-6*data["t_print"],1.0-data["e_print"][i],color=color,linestyle=linestyle,linewidth=linewidth,label=label)
            
                fontsize=18
                labelsize=18
                
            plot1.set_ylabel("$i_\mathrm{}\,(\mathrm{deg})$",fontsize=fontsize)
            plot2.set_ylabel("$1-e$",fontsize=fontsize)
            plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
           
            plots=[plot1,plot2]
            for plot in plots:
                plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

            plot2.set_ylim(5e-4,1.1e0)
            plot1.set_xticklabels([])
            
            ticks=plot1.get_yticks()
            plot1.set_yticks(ticks[2::])

            handles,labels = plot2.get_legend_handles_labels()
            plot2.legend(handles,labels,loc="lower left",fontsize=0.8*fontsize)

            fig.subplots_adjust(hspace=0.0,wspace=0.0)
            fig.savefig("example3.pdf")
            
            if args.show==True:
                pyplot.show()


### Functions below are used for example 2 ###
def reorientation_function(VRR_model,VRR_timescale,next_reorientation_time,orbit):

    print('='*50)
    print('reorientation_function')

    theta = 2.0*np.pi*randomf.random() - np.pi
    phi = 2.0*np.pi*randomf.random()
    r_hat_vec = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])

    if VRR_model == 1:

        Omega = 1.0/VRR_timescale
        
        orbit.VRR_Omega_vec_x = r_hat_vec[0]*Omega
        orbit.VRR_Omega_vec_y = r_hat_vec[1]*Omega
        orbit.VRR_Omega_vec_z = r_hat_vec[2]*Omega
    
    if VRR_model == 2:
        raise RuntimeError('VRR model 2 not supported')

    if VRR_model == 3:

        mu = 0.0
        sigma = (1.0/VRR_timescale)
        orbit.VRR_eta_20_init = orbit.VRR_eta_20_final
        orbit.VRR_eta_a_22_init = orbit.VRR_eta_a_22_final
        orbit.VRR_eta_b_22_init = orbit.VRR_eta_b_22_final
        orbit.VRR_eta_a_21_init = orbit.VRR_eta_a_21_final
        orbit.VRR_eta_b_21_init = orbit.VRR_eta_b_21_final

        orbit.VRR_eta_20_final = randomf.normal(mu,sigma)
        orbit.VRR_eta_a_22_final = randomf.normal(mu,sigma) 
        orbit.VRR_eta_b_22_final = randomf.normal(mu,sigma)
        orbit.VRR_eta_a_21_final = randomf.normal(mu,sigma)
        orbit.VRR_eta_b_21_final = randomf.normal(mu,sigma)
         
        orbit.VRR_initial_time = orbit.VRR_final_time
        orbit.VRR_final_time = next_reorientation_time

def compute_M_star_r(r,gamma,n_0,r_0,m_star):
    return 4.0*np.pi*(1.0/(3.0-gamma))*m_star*r_0**3*n_0*pow(r/r_0,3.0-gamma)

def compute_N_star_r(r,gamma,n_0,r_0,m_star):
    return compute_M_star_r(r,gamma,n_0,r_0,m_star)/m_star
    
def compute_n_star_r(r,gamma,n_0,r_0,m_star):
    return n_0*pow(r/r_0,-gamma)

def compute_rho_star_r(r,gamma,n_0,r_0,m_star):
    return m_star*compute_n_star_r(r,gamma,n_0,r_0,m_star)

def compute_sigma_r(r,gamma,n_0,r_0,m_star,M_MBH,CONST_G):
    return np.sqrt( CONST_G*M_MBH*(1.0/(r*(1.0+gamma)))*(1.0 + ((gamma+1.0)/(2.0*gamma-2.0))*compute_M_star_r(r,gamma,n_0,r_0,m_star)/M_MBH ) )
    

if __name__ == '__main__':
    args = parse_arguments()
    
    t=examples()
    print( 'Running example number',args.example,'; verbose =',args.verbose,'plot =',args.plot)
    function = getattr(t, 'example%s'%args.example)
    function(args)
