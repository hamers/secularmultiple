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
