"""
Some examples illustrating the usage of SecularMultiple
Adrian Hamers, June 2019
"""

import numpy as np
import numpy.random as randomf

from secularmultiple import SecularMultiple,Particle,Tools

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class examples():
    def example1(self):
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
        #I = rg*M*R**2
        #alpha = I/(mu*a0**2)
        T = R1**3/(CONST_G*m1*tau)
        t_V = 3.0*(1.0 + 1.0/k_L)*T
        #print 't_V',t_V,'M',M,'R',R

        particles[0].include_tidal_friction_terms = False
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = True
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-8

        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_viscous_time_scale = t_V
        particles[1].tides_gyration_radius = rg

        #binaries[0].include_1PN_terms = True
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
        
            print 't',t,'es',[o.e for o in orbits]
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
            plot2.plot(1.0e-6*t_print,INCL_print[i]*180.0/np.pi,color=color)
            
            plot1.set_xlabel("$t/\mathrm{Myr}$")
            plot2.set_xlabel("$t/\mathrm{Myr}$")
            plot1.set_ylabel("$r_i/\mathrm{AU}$")
            plot2.set_ylabel("$\mathrm{incl}_i/\mathrm{deg}$")
        fig.savefig("example1.pdf")
        pyplot.show()

    def example2(self):
        """
        Example of an (N-1)-planet system with constant spacing in mutual Hill sphere
        
        """

        code = SecularMultiple() ### initialize the code
        CONST_G = code.CONST_G ### extract physical constants from the code
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        RJup = 0.1027922358015816*CONST_R_SUN
        MJup = 0.0009546386983890755
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        N=5
        m0 = 1.0
        R0 = CONST_R_SUN
        a1 = 1.0
        mp = MJup
        Rp = RJup

        ei = 0.28
        ii = ei
        APi = 1.0e-10 ### arguments of periapsis
        LANi = 1.0e-10 ### longitudes of ascending node

        X = (1.0/2.0)*pow(2.0*mp/(3.0*m0),1.0/3.0)
        Delta = 10.0

        Delta_min = ei/X
        print 'Delta_min',Delta_min
        
        Nsteps = 100
        tend = 1.0e6
                
        masses = [m0]
        radii = [R0]
        semimajor_axes = []
        eccentricities = []
        inclinations = []
        APs = []
        LANs = []
        ai = a1        
        for i in range(N-1):
            masses.append(mp)
            radii.append(Rp)
            ai = ai*(1.0+Delta*X)/(1.0-Delta*X)
                
            semimajor_axes.append(ai)
            eccentricities.append(ei)
            inclinations.append(ii)
            APs.append(APi)
            LANs.append(LANi)
            
        print 'test',masses,semimajor_axes,eccentricities

        
        particles = Tools.create_nested_multiple(N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii=radii) 
        orbits = [x for x in particles if x.is_binary==True]
        for o in orbits:
            o.check_for_physical_collision_or_orbit_crossing = True
        N_orbits = len(orbits)
        

        #orbits[0].include_1PN_terms = True
        code.add_particles(particles)
        primary = code.particles[0]

        code.enable_tides = False
        code.enable_root_finding = True
        
        a_AU_print = [[] for x in range(N_orbits)]
        e_print = [[] for x in range(N_orbits)]
        INCL_print = [[] for x in range(N_orbits)]
        rel_INCL_print = [[] for x in range(N_orbits)]
        t_print = []
        
        t = 0.0
        dt = tend/float(Nsteps)
        import time
        start = time.time()
        while t<=tend:

            code.evolve_model(t)
            
            print 't',t,'es',[o.e for o in orbits]
            for i in range(N_orbits):
                rel_INCL_print[i].append(orbits[i].INCL_parent)
                a_AU_print[i].append(orbits[i].a)
                e_print[i].append(orbits[i].e)
                INCL_print[i].append(orbits[i].INCL)
            t_print.append(t)

            if code.flag == 2:
                t = code.model_time
                print 'root found at t=',t
                break

            t+=dt            
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
        colors=['k','r','g','y','b']
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
        fig.savefig("example2.pdf")

        pyplot.show()

    def example3(self):
        """
        Example of an (N-1)-planet system with constant spacing in mutual Hill sphere, Delta
        Calculate a series of systems with different Delta
        """

        code = SecularMultiple() ### initialize the code
        CONST_G = code.CONST_G ### extract physical constants from the code
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        RJup = 0.1027922358015816*CONST_R_SUN
        MJup = 0.0009546386983890755
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        N=5
        m0 = 1.0
        R0 = CONST_R_SUN
        a1 = 1.0
        mp = MJup
        Rp = RJup
        Porb1 = 2.0*np.pi*np.sqrt(a1**3/(CONST_G*(m0+mp)))

        ei = 0.35
        ii = ei
        APi = 1.0e-10
        LANi = 1.0e-10

        X = (1.0/2.0)*pow(2.0*mp/3.0,1.0/3.0)

        Delta_crit = ei/X
        print 'Delta_crit',Delta_crit
        Delta_min = np.amax([4.0,Delta_crit]) ### minimum allowed value of Delta to avoid orbit crossing at the onset
        Delta_max = 12.0 ### maximum Delta considered
        N_Delta = 10 ### number of points in Delta
        
        print 'Delta_min',Delta_min,'Delta_max',Delta_max
        
        Nsteps = 10
        tend = 1.0e7
                
        instability_times = []
        plot_Deltas = []
        Deltas = np.linspace(1.1*Delta_min,Delta_max,N_Delta)
        for index_Delta,Delta in enumerate(Deltas):
            print 'index_Delta',index_Delta,'Delta',Delta
            
            masses = [m0]
            radii = [R0]
            semimajor_axes = []
            eccentricities = []
            inclinations = []
            APs = []
            LANs = []
            ai = a1        
            for i in range(N-1):
                masses.append(mp)
                radii.append(Rp)

                semimajor_axes.append(ai)
                eccentricities.append(ei)
                inclinations.append(ii)
                APs.append(APi)
                LANs.append(LANi)

                ai = ai*(1.0+Delta*X)/(1.0-Delta*X)                
            
            instability_time = determine_stability_time(tend,Nsteps,N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii)
            instability_times.append(instability_time)
            plot_Deltas.append(Delta)
            
            print 'instability_time',instability_time
            
            if instability_time == tend:
                print 'instability time reached tend'
                break
                
        plot_Deltas = np.array(plot_Deltas)
        instability_times = np.array(instability_times)
        
        from matplotlib import pyplot
        fig=pyplot.figure(figsize=(8,6))
        plot1=fig.add_subplot(1,1,1,yscale="log")
        colors=['k','r','g','y','b']
        plot1.scatter(plot_Deltas,instability_times/Porb1,color='k')
        plot1.plot(plot_Deltas,instability_times/Porb1,color='k')
        plot1.set_xlabel("$\Delta$",fontsize=20)
        plot1.set_ylabel("$t_\mathrm{c}/t_0$",fontsize=20)
        plot1.set_title("$e_i=%s$"%ei,fontsize=20)
        
        fig.savefig("example3.pdf")
        pyplot.show()


    def example4(self):
        """
        Example of an (N-1)-planet system with constant spacing in mutual Hill sphere, Delta
        Calculate a series of systems with different Delta, and random arguments of periapsis and longitudes of the ascending node
        """

        code = SecularMultiple() ### initialize the code
        CONST_G = code.CONST_G ### extract physical constants from the code
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        RJup = 0.1027922358015816*CONST_R_SUN
        MJup = 0.0009546386983890755
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        N=5
        m0 = 1.0
        R0 = CONST_R_SUN
        a1 = 1.0
        mp = MJup
        Rp = RJup
        Porb1 = 2.0*np.pi*np.sqrt(a1**3/(CONST_G*(m0+mp)))

        ei = 0.35
        ii = ei
        #APi = 1.0e-10
        #LANi = 1.0e-10

        X = (1.0/2.0)*pow(2.0*mp/3.0,1.0/3.0)

        Delta_crit = ei/X
        print 'Delta_crit',Delta_crit
        Delta_min = np.amax([4.0,Delta_crit]) ### minimum allowed value of Delta to avoid orbit crossing at the onset
        Delta_max = 12.0 ### maximum Delta considered
        N_Delta = 4 ### number of points in Delta
        
        print 'Delta_min',Delta_min,'Delta_max',Delta_max
        
        Nsteps = 10
        tend = 1.0e7
        N_rand = 4
        
        instability_times_mean = []
        instability_times_std = []
        plot_Deltas = []
        Deltas = np.linspace(1.1*Delta_min,Delta_max,N_Delta)
        for index_Delta,Delta in enumerate(Deltas):
            print 'index_Delta',index_Delta,'Delta',Delta
            
            
            instability_times_temp = []
            for i_rand in range(N_rand):
                print '='*50
                print 'i_rand',i_rand
                
                randomf.seed(i_rand)
                
                masses = [m0]
                radii = [R0]
                semimajor_axes = []
                eccentricities = []
                inclinations = []
                APs = []
                LANs = []
                ai = a1        
                for i in range(N-1):
                    masses.append(mp)
                    radii.append(Rp)

                    semimajor_axes.append(ai)
                    eccentricities.append(ei)
                    inclinations.append(ii)
                    
                    APi = 2.0*np.pi*randomf.random()
                    APs.append(APi)
                    LANi = 2.0*np.pi*randomf.random()
                    LANs.append(LANi)

                    ai = ai*(1.0+Delta*X)/(1.0-Delta*X)                
                print 'masses',masses
                print 'radii',radii
                print 'semimajor_axes',semimajor_axes
                print 'eccentricities',eccentricities
                print 'inclinations',inclinations
                print 'APs',APs
                print 'LANs',LANs
            
                instability_time = determine_stability_time(tend,Nsteps,N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii)
                instability_times_temp.append(instability_time)
                
                print 'instability_time',instability_time
                
            if instability_time == tend:
                print 'instability time reached tend'
                break
                    
            instability_times_temp = np.array(instability_times_temp)
            
            instability_times_mean.append( np.mean( instability_times_temp ))
            instability_times_std.append( np.std( instability_times_temp ))

            plot_Deltas.append(Delta)

        instability_times_mean = np.array(instability_times_mean)
        instability_times_std = np.array(instability_times_std)
        plot_Deltas = np.array(plot_Deltas)
            
        #print 'instability_time',instability_time
        
        
        from matplotlib import pyplot
        fig=pyplot.figure(figsize=(8,6))
        plot1=fig.add_subplot(1,1,1,yscale="log")
        colors=['k','r','g','y','b']
        plot1.scatter(plot_Deltas,instability_times_mean/Porb1,color='k')
        plot1.errorbar(plot_Deltas,instability_times_mean/Porb1,yerr=instability_times_std/Porb1,color='k')
        
        plot1.set_xlabel("$\Delta$",fontsize=20)
        plot1.set_ylabel("$t_\mathrm{c}/t_0$",fontsize=20)
        plot1.set_title("$e_i=%s$"%ei,fontsize=20)
        
        fig.savefig("example4.pdf")
        pyplot.show()

    def example5(self):
        """
        Calculate Laplace-Lagrange (linear) solution for the Solar System
        
        """


        Mercury =   [0.38709927,      0.20563593,      7.00497902,      252.25032350,     77.45779628,     48.33076593]
        Venus =     [0.72333566,      0.00677672,      3.39467605,      181.97909950,    131.60246718,     76.67984255]
        Earth =     [1.00000261,      0.01671123,     -0.00001531,      100.46457166,    102.93768193,      0.0]
        Mars =      [1.52371034,      0.09339410,      1.84969142,       -4.55343205,    -23.94362959,     49.55953891]
        Jupiter =   [5.20288700,      0.04838624,      1.30439695,       34.39644051,     14.72847983,    100.47390909]
        Saturn =    [9.53667594,      0.05386179,      2.48599187,       49.95424423,     92.59887831,    113.66242448]
        Uranus =    [19.18916464,      0.04725744,      0.77263783,      313.23810451,    170.95427630,     74.01692503]
        Neptune =   [30.06992276,      0.00859048,      1.77004347,      -55.12002969,     44.96476227,    131.78422574]
        Pluto =     [39.48211675,      0.24882730,     17.14001206,      238.92903833,    224.06891629,    110.30393684]

        MJupiter_in_MSun = 0.0009546386983890755
        CONST_G = 4.0*np.pi**2
        
        names = ["$\mathrm{Mercury}$","$\mathrm{Venus}$","$\mathrm{Earth}$","$\mathrm{Mars}$","$\mathrm{Jupiter}$","$\mathrm{Saturn}$","$\mathrm{Uranus}$","$\mathrm{Neptune}$","$\mathrm{Pluto}$"]
        semimajor_axes = [Mercury[0], Venus[0], Earth[0], Mars[0], Jupiter[0], Saturn[0], Uranus[0], Neptune[0], Pluto[0]] ### AU
        eccentricities = [Mercury[1], Venus[1], Earth[1], Mars[1], Jupiter[1], Saturn[1], Uranus[1], Neptune[1], Pluto[1]]
        inclinations = np.array([Mercury[2], Venus[2], Earth[2], Mars[2], Jupiter[2], Saturn[2], Uranus[2], Neptune[2], Pluto[2]])*np.pi/180.0 ### radians
        longitudes_of_periapse = np.array([Mercury[4], Venus[4], Earth[4], Mars[4], Jupiter[4], Saturn[4], Uranus[4], Neptune[4], Pluto[4]])*np.pi/180.0 ### radians
        longitudes_of_ascending_node = np.array([Mercury[5], Venus[5], Earth[5], Mars[5], Jupiter[5], Saturn[5], Uranus[5], Neptune[5], Pluto[5]])*np.pi/180.0 ### radians
        arguments_of_pericenter = longitudes_of_periapse - longitudes_of_ascending_node
        
        
        planetary_masses = [0.00017, 	0.00256, 	0.00315, 	0.00034, 	1.0, 	0.299, 	0.046, 	0.054, 6.89626106e-6] ### MJupiter
        planetary_masses = [x*MJupiter_in_MSun for x in planetary_masses]
        M_star = 1.0

        N_p = len(planetary_masses)

        A_ij_array, B_ij_array, g_i, e_ij_bar, f_i, I_ij_bar = compute_Laplace_Lagrange_theory_quantities(CONST_G,M_star,planetary_masses,semimajor_axes)

        times_array_Myr = np.linspace(0.0,5.0,1000)
        times_array_yr = 1.0e6*times_array_Myr
        
        e_i_t, I_i_t, AP_i_t, LAN_i_t, h_i_t, k_i_t, p_i_t, q_i_t, e_max_i = compute_Laplace_Lagrange_theory_time_solution(CONST_G, \
            N_p, A_ij_array, B_ij_array, g_i, e_ij_bar, f_i, I_ij_bar, eccentricities, inclinations, arguments_of_pericenter, longitudes_of_ascending_node, times_array_yr)

        from matplotlib import pyplot

        #import distinct_colours
        
        pyplot.rc('text',usetex=True)
        pyplot.rc('legend',fancybox=True)          
        
        fig=pyplot.figure(figsize=(12,10))
        plot1=fig.add_subplot(2,1,1)
        plot2=fig.add_subplot(2,1,2)

        #colors = distinct_colours.get_distinct(12)
        colors = ['r','g','b','r','g','b','r','g','b']
        for index_planet in range(N_p):
            color=colors[index_planet]
            label = names[index_planet]
            plot1.plot(times_array_Myr,e_i_t[index_planet],label=label,color=color)
            plot2.plot(times_array_Myr,[x*180.0/np.pi for x in I_i_t[index_planet]],color=color)
        
        fontsize = 18
        labelsize = 18
        handles,labels = plot1.get_legend_handles_labels()
        plot1.legend(handles,labels,loc="upper right",fontsize=0.7*fontsize)

        plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
        plot1.set_ylabel("$e_i$",fontsize=fontsize)
        plot2.set_ylabel("$I_i/\mathrm{deg}$",fontsize=fontsize)

        plot1.tick_params(axis='both', which ='major', labelsize = labelsize)
        plot2.tick_params(axis='both', which ='major', labelsize = labelsize)
        
        ticks = plot1.get_yticks()
        plot1.set_yticks(ticks[1:])
        plot1.set_xticklabels([])
        
        fig.subplots_adjust(hspace=0.0,wspace=0.0)
        fig.savefig("example5.pdf")
        
        pyplot.show()

def compute_Laplace_Lagrange_theory_quantities(CONST_G,M,m_i_array,a_i_array):
    N_p = len(m_i_array)

    n_i_array = np.zeros((N_p))
    alpha_ij_array = np.zeros((N_p,N_p))
    alpha_bar_ij_array = np.zeros((N_p,N_p))
    epsilon_ij_array = np.zeros((N_p,N_p))
    for i,a_i in enumerate(a_i_array):
        m_i = m_i_array[i]
        n_i_array[i] = np.sqrt( CONST_G*(M + m_i)/(a_i**3) )
        
        for j,a_j in enumerate(a_i_array):
            m_j = m_i_array[j]

            epsilon_ij_array[i][j] = m_j/(M + m_i)

            if a_j > a_i:
                alpha_ij_array[i][j] = a_i/a_j
                alpha_bar_ij_array[i][j] = a_i/a_j
            elif a_j < a_i:
                alpha_ij_array[i][j] = a_j/a_i
                alpha_bar_ij_array[i][j] = 1.0

    A_ij_array = np.zeros((N_p,N_p))
    B_ij_array = np.zeros((N_p,N_p))
    for i in range(N_p):
        A_sum_j = 0.0
        n_i = n_i_array[i]
        
        for j in range(N_p):
            alpha_ij = alpha_ij_array[i][j]
            alpha_bar_ij = alpha_bar_ij_array[i][j]
            epsilon_ij = epsilon_ij_array[i][j]
            
            if j != i:
                A_sum_j += (n_i/4.0)*epsilon_ij*alpha_ij*alpha_bar_ij*laplace_coefficient(3.0/2.0,1.0,alpha_ij)
                A_ij_array[i][j] = -(n_i/4.0)*epsilon_ij*alpha_ij*alpha_bar_ij*laplace_coefficient(3.0/2.0,2.0,alpha_ij)
                B_ij_array[i][j] = (n_i/4.0)*epsilon_ij*alpha_ij*alpha_bar_ij*laplace_coefficient(3.0/2.0,1.0,alpha_ij)

        A_ij_array[i][i] = A_sum_j
        B_ij_array[i][i] = -A_ij_array[i][i]

    g_i,e_ij_bar = np.linalg.eig(A_ij_array)
    f_i,I_ij_bar = np.linalg.eig(B_ij_array)

    return A_ij_array, B_ij_array, g_i, e_ij_bar, f_i, I_ij_bar
    
def laplace_coefficient(s,j,alpha):
    from scipy import integrate
    args = s,j,alpha
    integral = (1.0/np.pi)*integrate.quad(laplace_coefficient_integrand, 0.0, 2.0*np.pi, args=(args,),epsrel=1.0e-10,epsabs=1.0e-10,limit=500)[0]
    return integral

def laplace_coefficient_integrand(x,args):
    s,j,alpha = args
    return np.cos(j*x)*pow(1.0 - 2.0*alpha*np.cos(x) + alpha**2, -s)

def compute_Laplace_Lagrange_theory_time_solution(CONST_G,N_p, A_ij_array, B_ij_array, g_i, e_ij_bar, f_i, I_ij_bar, e_0s, INCL_0s, AP_0s, LANS_0s, times_array):

    h_i_0 = np.zeros((N_p))
    k_i_0 = np.zeros((N_p))
    p_i_0 = np.zeros((N_p))
    q_i_0 = np.zeros((N_p))
    
    for i in range(N_p):
        e_0 = e_0s[i]
        INCL_0 = INCL_0s[i]
        AP_0 = AP_0s[i]
        LAN_0 = LANS_0s[i]
        
        h_i_0[i] = e_0*np.sin(AP_0 + LAN_0)
        k_i_0[i] = e_0*np.cos(AP_0 + LAN_0)

        p_i_0[i] = INCL_0*np.sin(LAN_0)
        q_i_0[i] = INCL_0*np.cos(LAN_0)

    inv_e_ij_bar = np.linalg.inv(e_ij_bar)
    S_i_sin_beta_i = np.dot(inv_e_ij_bar,h_i_0)
    S_i_cos_beta_i = np.dot(inv_e_ij_bar,k_i_0)
    
    tan_beta_i = S_i_sin_beta_i/S_i_cos_beta_i
    beta_i = np.arctan(tan_beta_i)
    S_i = S_i_cos_beta_i/np.cos(beta_i)
    e_ij = S_i*e_ij_bar

    inv_I_ij_bar = np.linalg.inv(I_ij_bar)
    T_i_sin_gamma_i = np.dot(inv_I_ij_bar,p_i_0)
    T_i_cos_gamma_i = np.dot(inv_I_ij_bar,q_i_0)

    tan_gamma_i = T_i_sin_gamma_i/T_i_cos_gamma_i
    gamma_i = np.arctan(tan_gamma_i)
    T_i = T_i_cos_gamma_i/np.cos(gamma_i)
    I_ij = T_i*I_ij_bar

    N_t = len(times_array)

    h_i_t = np.zeros((N_p,N_t))
    k_i_t = np.zeros((N_p,N_t))    
    p_i_t = np.zeros((N_p,N_t))
    q_i_t = np.zeros((N_p,N_t))    

    e_max_i = np.zeros((N_p))

    for i in range(N_p):
       
        for index_t, t in enumerate(times_array):
            h_sum_l = 0.0
            k_sum_l = 0.0

            p_sum_l = 0.0
            q_sum_l = 0.0

            e_max_i_partial = 0.0
            for l in range(N_p):
                e_il = e_ij[i][l]
                g_l = g_i[l]
                beta_l = beta_i[l]

                I_il = I_ij[i][l]
                f_l = f_i[l]
                gamma_l = gamma_i[l]

                h_sum_l += e_il*np.sin(g_l*t + beta_l)
                k_sum_l += e_il*np.cos(g_l*t + beta_l)

                p_sum_l += I_il*np.sin(f_l*t + gamma_l)
                q_sum_l += I_il*np.cos(f_l*t + gamma_l)

                e_max_i_partial += np.fabs(e_il)

            e_max_i[i] = e_max_i_partial
        
            h_i_t[i][index_t] = h_sum_l
            k_i_t[i][index_t] = k_sum_l

            p_i_t[i][index_t] = p_sum_l
            q_i_t[i][index_t] = q_sum_l
   
    e_i_t = np.sqrt(h_i_t**2 + k_i_t**2)
    I_i_t = np.sqrt(p_i_t**2 + q_i_t**2)

    sin_LAN_i_t = p_i_t/I_i_t
    cos_LAN_i_t = q_i_t/I_i_t
    LAN_i_t = np.arctan2( sin_LAN_i_t, cos_LAN_i_t)
    
    sin_LOP_i_t = h_i_t/e_i_t 
    cos_LOP_i_t = k_i_t/e_i_t 
    LOP_i_t = np.arctan2( sin_LOP_i_t, cos_LOP_i_t)
    AP_i_t = LOP_i_t - LAN_i_t

    return e_i_t, I_i_t, AP_i_t, LAN_i_t, h_i_t, k_i_t, p_i_t, q_i_t,e_max_i
    

def determine_stability_time(tend,Nsteps,N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii):

    particles = Tools.create_nested_multiple(N, masses,semimajor_axes,eccentricities,inclinations,APs,LANs,radii=radii) 
    orbits = [x for x in particles if x.is_binary==True]
    for o in orbits:
        o.check_for_physical_collision_or_orbit_crossing = True
    N_orbits = len(orbits)
    
    #binaries[0].include_1PN_terms = True
    code = SecularMultiple() ### initialize the code

    code.add_particles(particles)
    primary = code.particles[0]

    code.enable_tides = False
    code.enable_root_finding = True
    
    a_AU_print = [[] for x in range(N_orbits)]
    e_print = [[] for x in range(N_orbits)]
    INCL_print = [[] for x in range(N_orbits)]
    rel_INCL_print = [[] for x in range(N_orbits)]
    t_print = []
    
    t = 0.0
    dt = tend/float(Nsteps)
    import time
    start = time.time()
    while t<tend:
        t+=dt            
        code.evolve_model(t)
    
        #print 't',t,'es',[o.e for o in orbits]
        for i in range(N_orbits):
            rel_INCL_print[i].append(orbits[i].INCL_parent)
            a_AU_print[i].append(orbits[i].a)
            e_print[i].append(orbits[i].e)
            INCL_print[i].append(orbits[i].INCL)
        t_print.append(t)

        if code.flag == 2:
            t = code.model_time
            #print 'root found at t=',t
            break

        
    #print 'wall time',time.time()-start
    code.reset()

    return t
        

    
if __name__ == '__main__':
    import sys
    t=examples()
    if len(sys.argv)>1:
        i = int(sys.argv[1])
        if i<0:
            t.development_test()
        else:
            print('Running example %s'%i)
            function = getattr(t, 'example%s'%i)
            function()
