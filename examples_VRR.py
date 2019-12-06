"""
Some examples illustrating the usage of SecularMultiple
Adrian Hamers, June 2019
"""

import numpy as np
import numpy.random as randomf

import pickle

from secularmultiple import SecularMultiple,Particle,Tools

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class examples():
    def example1(self):
        """
        Lidov-Kozai problem of a planet around a star around a supermassive black hole.
        Includes perturbations from other stars in the form of vector resonant relaxation (VRR)
        
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
        meter = 1.0/1.496e+11

        ### Input parameters ###
        m1 = 1.0 ### stellar mass
        m2 = MJup ### planetary mass
        m3 = 4.0e6 ### GC MBH mass
        
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
        
        rho = C*r^{-gamma}
        
        
        c1 = 4.8
        c2 = -2.9
        log10_sigma_h_km_s = (np.log10(m3) - c2)/c1
        
        sigma_h_km_s = pow(10.0,log10_sigma_h_km_s)
        sigma_h = 1.0e3*sigma_h_km_s*meter/second
        print 'sigma_h_km_s',sigma_h_km_s,'sigma_h',sigma_h

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
        
        print 'n_star',n_star,'M_star',M_star,'N_star',N_star,'sigma_r',sigma_r
        
        LK_timescale = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
        print 'LK_timescale',LK_timescale
        VRR_mass_precession_timescale = (1.0/2.0)*pow(1.0-e2**2,-1.0/2.0)*(m3/M_star)*P2
        VRR_mass_precession_rate = 1.0/VRR_mass_precession_timescale
        VRR_timescale = (P2/2.0)*(m3/m_star)*1.0/np.sqrt(N_star)
        #VRR_timescale *= 0.1
                
        print 'VRR_mass_precession_timescale',VRR_mass_precession_timescale,'VRR_timescale',VRR_timescale
       
        outer_orbit.VRR_include_mass_precession = VRR_include_mass_precession
        outer_orbit.VRR_mass_precession_rate = VRR_mass_precession_rate

        VRR_reorientation_timestep = np.sqrt(0.1)*VRR_timescale
        print 'VRR_reorientation_timestep',VRR_reorientation_timestep

        outer_orbit.VRR_model = VRR_model
        reorientation_function(VRR_model,VRR_timescale,VRR_reorientation_timestep,outer_orbit)
        

        v_bin = np.sqrt(CONST_G*(m1+m2)/a1)
        q_sigma = (m1+m2)/m_star
        log_Lambda = np.log( 3.0*((1.0 + 1.0/q_sigma)/(1.0 + 2.0/q_sigma))*sigma_r**2/v_bin**2 )
        evaporation_timescale = np.sqrt( (1.0+q_sigma)/(2.0*np.pi*q_sigma) )*(m1+m2)*sigma_r/(8.0*np.sqrt(np.pi)*CONST_G*a1*m_star**2*n_star*log_Lambda)
        print 'evaporation_timescale',evaporation_timescale
        

        inner_orbit.include_pairwise_1PN_terms = include_inner_1PN_terms
        outer_orbit.include_pairwise_1PN_terms = include_outer_1PN_terms
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

            print 't',t,'es',[o.e for o in orbits],'Omegas',[o.LAN for o in orbits]
            for i in range(N_orbits):
                rel_INCL_print[i].append(orbits[i].INCL_parent)
                a_AU_print[i].append(orbits[i].a)
                e_print[i].append(orbits[i].e)
                INCL_print[i].append(orbits[i].INCL)
            t_print.append(t)
            
        print 'wall time',time.time()-start
        
        t_print = np.array(t_print)
        for i in range(N_orbits):
            INCL_print[i] = np.array(INCL_print[i])
            rel_INCL_print[i] = np.array(rel_INCL_print[i])
            e_print[i] = np.array(e_print[i])
            a_AU_print[i] = np.array(a_AU_print[i])
        
        from matplotlib import pyplot
        fig=pyplot.figure(figsize=(8,8))
        plot1=fig.add_subplot(2,1,1,yscale="log")
        plot2=fig.add_subplot(2,1,2,yscale="linear")
        colors=['k','r','g']
        for i in range(N_orbits):
            color=colors[i]
            plot1.plot(1.0e-6*t_print,a_AU_print[i],color=color)
            plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0-e_print[i]),color=color)
            plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0+e_print[i]),color=color)
            #plot2.plot(1.0e-6*t_print,INCL_print[i]*180.0/np.pi,color=color,linestyle='dotted')
            plot2.plot(1.0e-6*t_print,rel_INCL_print[i]*180.0/np.pi,color=color)
            
            plot1.set_xlabel("$t/\mathrm{Myr}$")
            plot2.set_xlabel("$t/\mathrm{Myr}$")
            plot1.set_ylabel("$r_i/\mathrm{AU}$")
            plot2.set_ylabel("$\mathrm{incl}_\mathrm{rel}/\mathrm{deg}$")
        fig.savefig("example1.pdf")
        pyplot.show()

    def example2(self):
        """
        Lidov-Kozai problem of a planet around an S-star-like star in the Galactic Center (GC)
        Includes perturbations from other stars in the form of vector resonant relaxation (VRR)
        Integrate for the duration of the expected lifetime of the star+planet binary due to incoherent encounters
        Check for physical collisions of the planet with the star
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
        m1 = 1.0 ### stellar mass
        m2 = MJup ### planetary mass
        m3 = 4.0e6 ### GC MBH mass
        
        a1_min = 1.0e-2
        a1_max = 1.0e0
        N_a1 = 5
        a2 = 1.0e4 ### typical for S-star 
        e1 = 0.01
        e2 = 0.1
        i1 = 0.01*np.pi/180.0
        i2 = 55.0*np.pi/180.0
        AP1 = 0.01*np.pi/180.0
        AP2 = 0.01*np.pi/180.0
        LAN1 = 0.01*np.pi/180.0
        LAN2 = 0.01*np.pi/180.0

        R1 = 1.0*CONST_R_SUN
        R2 = 1.0*RJup
        R3 = CONST_G*m3/(CONST_C**2)

        m_star = 1.0 ### average stellar mass in the background
        gamma = 3.0/2.0 ### slope of stellar background
        VRR_model = 3

        tmax = 1.0e8 ### maximum integration time

        ### Simulation parameters ###
        VRR_include_mass_precession = True
        include_inner_1PN_terms = False
        include_outer_1PN_terms = True
        
        a1_values = pow(10.0,np.linspace(np.log10(a1_min),np.log10(a1_max),N_a1))
        for index_a1,a1 in enumerate(a1_values):
            
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

            for o in orbits:
                o.check_for_physical_collision_or_orbit_crossing = True

            inner_orbit = orbits[0]
            outer_orbit = orbits[1]
            
            c1 = 4.8
            c2 = -2.9
            log10_sigma_h_km_s = (np.log10(m3) - c2)/c1
            
            sigma_h_km_s = pow(10.0,log10_sigma_h_km_s)
            sigma_h = 1.0e3*sigma_h_km_s*meter/second
            print 'sigma_h_km_s',sigma_h_km_s,'sigma_h',sigma_h

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
            
            print 'n_star',n_star,'M_star',M_star,'N_star',N_star,'sigma_r',sigma_r
        
            LK_timescale = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
            print 'LK_timescale',LK_timescale
            VRR_mass_precession_timescale = (1.0/2.0)*pow(1.0-e2**2,-1.0/2.0)*(m3/M_star)*P2
            VRR_mass_precession_rate = 1.0/VRR_mass_precession_timescale
            VRR_timescale = (P2/2.0)*(m3/m_star)*1.0/np.sqrt(N_star)
                
            print 'VRR_mass_precession_timescale',VRR_mass_precession_timescale,'VRR_timescale',VRR_timescale
           
            outer_orbit.VRR_include_mass_precession = VRR_include_mass_precession
            outer_orbit.VRR_mass_precession_rate = VRR_mass_precession_rate

            VRR_reorientation_timestep = np.sqrt(0.1)*VRR_timescale
            print 'VRR_reorientation_timestep',VRR_reorientation_timestep

            outer_orbit.VRR_model = VRR_model
            reorientation_function(VRR_model,VRR_timescale,VRR_reorientation_timestep,outer_orbit)

            v_bin = np.sqrt(CONST_G*(m1+m2)/a1)
            q_sigma = (m1+m2)/m_star
            log_Lambda = np.log( 3.0*((1.0 + 1.0/q_sigma)/(1.0 + 2.0/q_sigma))*sigma_r**2/v_bin**2 )
            evaporation_timescale = np.sqrt( (1.0+q_sigma)/(2.0*np.pi*q_sigma) )*(m1+m2)*sigma_r/(8.0*np.sqrt(np.pi)*CONST_G*a1*m_star**2*n_star*log_Lambda)
            print 'evaporation_timescale',evaporation_timescale
        
            inner_orbit.include_pairwise_1PN_terms = include_inner_1PN_terms
            outer_orbit.include_pairwise_1PN_terms = include_outer_1PN_terms

            tend = np.amin([tmax,evaporation_timescale])
            Nsteps = 1000
            data = run_simulation(tend,Nsteps,particles, \
                VRR_reorientation_timestep,VRR_model,VRR_timescale,outer_orbit, \
                enable_tides=False,enable_root_finding=True,enable_VRR=True)
            found_root,t_print,rel_INCL_print,e_print,a_AU_print,rp_AU_print = data
            
            print 'found_root',found_root
            
            from matplotlib import pyplot
            fig=pyplot.figure(figsize=(8,8))
            plot1=fig.add_subplot(2,1,1,yscale="log")
            plot2=fig.add_subplot(2,1,2,yscale="linear")
            colors=['k','r','g']
            for i in range(N_orbits):
                color=colors[i]
                plot1.plot(1.0e-6*t_print,a_AU_print[i],color=color)
                plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0-e_print[i]),color=color)
                plot1.plot(1.0e-6*t_print,a_AU_print[i]*(1.0+e_print[i]),color=color)
                #plot2.plot(1.0e-6*t_print,INCL_print[i]*180.0/np.pi,color=color,linestyle='dotted')
                plot2.plot(1.0e-6*t_print,rel_INCL_print[i]*180.0/np.pi,color=color)
                
                incl_rel = Tools.compute_mutual_inclination(i1,i2,LAN1,LAN2)
                emax_LK = np.sqrt(1.0 - (5.0/3.0)*np.cos(incl_rel)**2)
                plot1.axhline(y = a1*(1.0 - emax_LK),color='g',linestyle='dotted')
                plot1.axhline(y = R1+R2,color='r')
                plot1.set_xlabel("$t/\mathrm{Myr}$")
                plot2.set_xlabel("$t/\mathrm{Myr}$")
                plot1.set_ylabel("$r_i/\mathrm{AU}$")
                plot2.set_ylabel("$\mathrm{incl}_\mathrm{rel}/\mathrm{deg}$")
#            fig.savefig("example2.pdf")
            pyplot.show()


    def example3(self):
        """
        Similar to example 2, now recording maximum eccentricity e_max at different inner semimajor axes a1, and plotting
        e_max as a function of a1
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
        m1 = 1.0 ### stellar mass
        m2 = MJup ### planetary mass
        m3 = 4.0e6 ### GC MBH mass
        
        a1_min = 1.0e-2
        a1_max = 1.0e0
        N_a1 = 10
        a2 = 1.0e4 ### typical for S-star 
        e1 = 0.01
        e2 = 0.1
        i1 = 0.01*np.pi/180.0
        i2 = 55.0*np.pi/180.0
        AP1 = 0.01*np.pi/180.0
        AP2 = 0.01*np.pi/180.0
        LAN1 = 0.01*np.pi/180.0
        LAN2 = 0.01*np.pi/180.0

        R1 = 1.0*CONST_R_SUN
        R2 = 1.0*RJup
        R3 = CONST_G*m3/(CONST_C**2)

        m_star = 1.0 ### average stellar mass in the background
        gamma = 3.0/2.0 ### slope of stellar background
        VRR_model = 3

        tmax = 1.0e8 ### maximum integration time

        ### Simulation parameters ###
        VRR_include_mass_precession = True
        include_inner_1PN_terms = False
        include_outer_1PN_terms = True
        
        
        e_maxs = []
        a1_values = pow(10.0,np.linspace(np.log10(a1_min),np.log10(a1_max),N_a1))
        for index_a1,a1 in enumerate(a1_values):
            print '='*50
            print 'a1/AU',a1
            
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

            for o in orbits:
                o.check_for_physical_collision_or_orbit_crossing = True

            inner_orbit = orbits[0]
            outer_orbit = orbits[1]
            
            c1 = 4.8
            c2 = -2.9
            log10_sigma_h_km_s = (np.log10(m3) - c2)/c1
            
            sigma_h_km_s = pow(10.0,log10_sigma_h_km_s)
            sigma_h = 1.0e3*sigma_h_km_s*meter/second
            print 'sigma_h_km_s',sigma_h_km_s,'sigma_h',sigma_h

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
            
            print 'n_star',n_star,'M_star',M_star,'N_star',N_star,'sigma_r',sigma_r
        
            LK_timescale = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
            print 'LK_timescale/Myr',LK_timescale*1e-6
            VRR_mass_precession_timescale = (1.0/2.0)*pow(1.0-e2**2,-1.0/2.0)*(m3/M_star)*P2
            VRR_mass_precession_rate = 1.0/VRR_mass_precession_timescale
            VRR_timescale = (P2/2.0)*(m3/m_star)*1.0/np.sqrt(N_star)
                
            print 'VRR_mass_precession_timescale/Myr',VRR_mass_precession_timescale*1e-6,'VRR_timescale/Myr',VRR_timescale*1e-6
           
            outer_orbit.VRR_include_mass_precession = VRR_include_mass_precession
            outer_orbit.VRR_mass_precession_rate = VRR_mass_precession_rate

            VRR_reorientation_timestep = np.sqrt(0.1)*VRR_timescale
            print 'VRR_reorientation_timestep/Myr',VRR_reorientation_timestep*1e-6

            outer_orbit.VRR_model = VRR_model
            reorientation_function(VRR_model,VRR_timescale,VRR_reorientation_timestep,outer_orbit)

            v_bin = np.sqrt(CONST_G*(m1+m2)/a1)
            q_sigma = (m1+m2)/m_star
            log_Lambda = np.log( 3.0*((1.0 + 1.0/q_sigma)/(1.0 + 2.0/q_sigma))*sigma_r**2/v_bin**2 )
            evaporation_timescale = np.sqrt( (1.0+q_sigma)/(2.0*np.pi*q_sigma) )*(m1+m2)*sigma_r/(8.0*np.sqrt(np.pi)*CONST_G*a1*m_star**2*n_star*log_Lambda)
            print 'evaporation_timescale/Myr',evaporation_timescale*1e-6
        
            inner_orbit.include_pairwise_1PN_terms = include_inner_1PN_terms
            outer_orbit.include_pairwise_1PN_terms = include_outer_1PN_terms

            tend = np.amin([tmax,evaporation_timescale])
            Nsteps = 10000
            data = run_simulation(tend,Nsteps,particles, \
                VRR_reorientation_timestep,VRR_model,VRR_timescale,outer_orbit, \
                enable_tides=False,enable_root_finding=True,enable_VRR=True)
            found_root,t_print,rel_INCL_print,e_print,a_AU_print,rp_AU_print = data
            
            if found_root==True:
                e_maxs.append(1.0)
            else:
                e_maxs.append(np.amax(e_print))
                
            print 'found_root',found_root
            
        from matplotlib import pyplot
        fig=pyplot.figure(figsize=(8,8))
        plot1=fig.add_subplot(1,1,1,yscale="linear")
#        plot2=fig.add_subplot(2,1,2,yscale="linear")
        colors=['k','r','g']
        for i in range(N_orbits):
            color=colors[i]
            plot1.plot(a1_values,e_maxs,color='k')
            plot1.set_xlabel("$a_1/\mathrm{AU}$",fontsize=18)
            plot1.set_ylabel("$e_\mathrm{max}$",fontsize=18)
#            fig.savefig("example2.pdf")
        pyplot.show()


def run_simulation(tend,Nsteps,particles, \
    VRR_reorientation_timestep,VRR_model,VRR_timescale,outer_orbit, \
    enable_tides=False,enable_root_finding=True,enable_VRR=True):
    code = SecularMultiple() ### initialize the code
    
    code.add_particles(particles)
    primary = code.particles[0]

    code.enable_tides = enable_tides
    code.enable_root_finding = enable_root_finding
    code.enable_VRR = enable_VRR
            
    orbits = [x for x in particles if x.is_binary==True]
    N_orbits = len(orbits)

    a_AU_print = [[] for x in range(N_orbits)]
    e_print = [[] for x in range(N_orbits)]
    rp_AU_print = [[] for x in range(N_orbits)]
    INCL_print = [[] for x in range(N_orbits)]
    rel_INCL_print = [[] for x in range(N_orbits)]
    t_print = []
            
    t = 0.0
    
    dt_fixed = tend/float(Nsteps)
    t_next_reorientation = VRR_reorientation_timestep
        
    found_root = False
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

        #print 't',t,'es',[o.e for o in orbits],'Omegas',[o.LAN for o in orbits]
        
        
        if code.flag == 2:
            t = code.model_time
            print 'root found at t=',t
            found_root = True

        for i in range(N_orbits):
            rel_INCL_print[i].append(orbits[i].INCL_parent)
            a_AU_print[i].append(orbits[i].a)
            e_print[i].append(orbits[i].e)
            INCL_print[i].append(orbits[i].INCL)
            rp_AU_print[i].append( orbits[i].a*(1.0-orbits[i].e) )
        t_print.append(t)
        
        if found_root == True:
            break

    print 'wall time',time.time()-start
    
    t_print = np.array(t_print)
    for i in range(N_orbits):
        INCL_print[i] = np.array(INCL_print[i])
        rel_INCL_print[i] = np.array(rel_INCL_print[i])
        e_print[i] = np.array(e_print[i])
        a_AU_print[i] = np.array(a_AU_print[i])
        rp_AU_print[i] = np.array(rp_AU_print[i])

    code.reset()

    data = found_root,t_print,rel_INCL_print,e_print,a_AU_print,rp_AU_print

    return data
        

def reorientation_function(VRR_model,VRR_timescale,next_reorientation_time,orbit):

    #print('='*50)
    #print('reorientation_function')

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
    
def save(data,filename):
    with open(filename,'wb') as file:
        pickle.dump(data,file)

def load(filename):
    with open('output.pkl','rb') as file:
        data = pickle.load(file)
    return data



def run_and_save_simulation(filename):
    """
    Similar to example 2, now recording maximum eccentricity e_max at different inner semimajor axes a1, and plotting
    e_max as a function of a1
    """

    print 'run_simulation'
    
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
    m1 = 1.0 ### stellar mass
    m2 = MJup ### planetary mass
    m3 = 4.0e6 ### GC MBH mass
    
    a1_min = 1.0e-2
    a1_max = 1.0e0
    N_a1 = 2
    a2 = 1.0e4 ### typical for S-star 
    e1 = 0.01
    e2 = 0.1
    i1 = 0.01*np.pi/180.0
    i2 = 55.0*np.pi/180.0
    AP1 = 0.01*np.pi/180.0
    AP2 = 0.01*np.pi/180.0
    LAN1 = 0.01*np.pi/180.0
    LAN2 = 0.01*np.pi/180.0

    R1 = 1.0*CONST_R_SUN
    R2 = 1.0*RJup
    R3 = CONST_G*m3/(CONST_C**2)

    m_star = 1.0 ### average stellar mass in the background
    gamma = 3.0/2.0 ### slope of stellar background
    VRR_model = 3

    tmax = 1.0e2 ### maximum integration time

    ### Simulation parameters ###
    VRR_include_mass_precession = True
    include_inner_1PN_terms = False
    include_outer_1PN_terms = True
    
    data_container = []
    e_maxs = []
    a1_values = pow(10.0,np.linspace(np.log10(a1_min),np.log10(a1_max),N_a1))
    for index_a1,a1 in enumerate(a1_values):
        print '='*50
        print 'a1/AU',a1
        
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

        for o in orbits:
            o.check_for_physical_collision_or_orbit_crossing = True

        inner_orbit = orbits[0]
        outer_orbit = orbits[1]
        
        c1 = 4.8
        c2 = -2.9
        log10_sigma_h_km_s = (np.log10(m3) - c2)/c1
        
        sigma_h_km_s = pow(10.0,log10_sigma_h_km_s)
        sigma_h = 1.0e3*sigma_h_km_s*meter/second
        print 'sigma_h_km_s',sigma_h_km_s,'sigma_h',sigma_h

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
        
        print 'n_star',n_star,'M_star',M_star,'N_star',N_star,'sigma_r',sigma_r
    
        LK_timescale = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
        print 'LK_timescale/Myr',LK_timescale*1e-6
        VRR_mass_precession_timescale = (1.0/2.0)*pow(1.0-e2**2,-1.0/2.0)*(m3/M_star)*P2
        VRR_mass_precession_rate = 1.0/VRR_mass_precession_timescale
        VRR_timescale = (P2/2.0)*(m3/m_star)*1.0/np.sqrt(N_star)
            
        print 'VRR_mass_precession_timescale/Myr',VRR_mass_precession_timescale*1e-6,'VRR_timescale/Myr',VRR_timescale*1e-6
       
        outer_orbit.VRR_include_mass_precession = VRR_include_mass_precession
        outer_orbit.VRR_mass_precession_rate = VRR_mass_precession_rate

        VRR_reorientation_timestep = np.sqrt(0.1)*VRR_timescale
        print 'VRR_reorientation_timestep/Myr',VRR_reorientation_timestep*1e-6

        outer_orbit.VRR_model = VRR_model
        reorientation_function(VRR_model,VRR_timescale,VRR_reorientation_timestep,outer_orbit)

        v_bin = np.sqrt(CONST_G*(m1+m2)/a1)
        q_sigma = (m1+m2)/m_star
        log_Lambda = np.log( 3.0*((1.0 + 1.0/q_sigma)/(1.0 + 2.0/q_sigma))*sigma_r**2/v_bin**2 )
        evaporation_timescale = np.sqrt( (1.0+q_sigma)/(2.0*np.pi*q_sigma) )*(m1+m2)*sigma_r/(8.0*np.sqrt(np.pi)*CONST_G*a1*m_star**2*n_star*log_Lambda)
        print 'evaporation_timescale/Myr',evaporation_timescale*1e-6
    
        inner_orbit.include_pairwise_1PN_terms = include_inner_1PN_terms
        outer_orbit.include_pairwise_1PN_terms = include_outer_1PN_terms

        tend = np.amin([tmax,evaporation_timescale])
        Nsteps = 10000
        data = run_simulation(tend,Nsteps,particles, \
            VRR_reorientation_timestep,VRR_model,VRR_timescale,outer_orbit, \
            enable_tides=False,enable_root_finding=True,enable_VRR=True)
        found_root,t_print,rel_INCL_print,e_print,a_AU_print,rp_AU_print = data
        
        data_container.append(data)
        if found_root==True:
            e_maxs.append(1.0)
        else:
            e_maxs.append(np.amax(e_print))
            
        print 'found_root',found_root
        
    data = a1_values,e_maxs
    save(data,filename)
    


def plot_data(filename):
    data = load(filename)
    a1_values,e_maxs = data
    
    N_orbits = 2
    
    irel = 55.0*np.pi/180.0
    emax_LK = np.sqrt(1.0 - (5.0/3.0)*np.cos(irel)**2)

    code = SecularMultiple() ### initialize the code
    CONST_G = code.CONST_G ### extract physical constants from the code
    CONST_C = code.CONST_C
    CONST_R_SUN = code.CONST_R_SUN
    RJup = 0.1027922358015816*CONST_R_SUN
    MJup = 0.0009546386983890755
    day = 1.0/365.25
    second = day/(24.0*3600.0)
    meter = 1.0/1.496e+11
    
    R1 = 1.0*CONST_R_SUN
    R2 = 1.0*RJup

    e1_col = 1.0 - (R1+R2)/a1_values
    #e1_col = [1.0 - (R1+R2)/a1 for a1 in a1_values]

    from matplotlib import pyplot
    fig=pyplot.figure(figsize=(8,8))
    plot1=fig.add_subplot(1,1,1,yscale="linear")
#        plot2=fig.add_subplot(2,1,2,yscale="linear")
    colors=['k','r','g']
    for i in range(N_orbits):
        color=colors[i]
        plot1.plot(a1_values,e_maxs,color='k')
        plot1.set_xlabel("$a_1/\mathrm{AU}$",fontsize=18)
        plot1.set_ylabel("$e_\mathrm{max}$",fontsize=18)
    
    plot1.plot(a1_values,e1_col,color='r')
    plot1.axhline(y=emax_LK,color='k',linestyle='dotted')
    
    fig.savefig(filename + ".pdf")
    
    pyplot.show()


if __name__ == '__main__':
    import sys
    t=examples()
    if len(sys.argv)>1:
        i = int(sys.argv[1])
        if i<0:
            t.development_test()
        
        filename = "reference"
        if i==1:
            run_and_save_simulation(filename)
        if i==2:
            plot_data(filename)
            
        #else:
        #    print('Running example %s'%i)
        #    function = getattr(t, 'example%s'%i)
        #    function()
