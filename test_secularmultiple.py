import numpy as np
import argparse
import time
        
from secularmultiple import SecularMultiple,Particle,Tools

"""
Several routines for testing the code/installation. 
To run all tests, simply run `python test_secularmultiple.py'.
Specific tests can be run with the command line --t i, where i is the
number of the test.

Adrian Hamers, March 2020
"""

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

    parser.add_argument("--t",                           type=int,     dest="test",                        default=0,              help="Test number")
    
    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                         default=False,         help="verbose terminal output")
    add_bool_arg(parser, 'plot',                            default=False,         help="make plots")
    
    args = parser.parse_args()

    return args

class test_secularmultiple():

    def test1(self,args):
        print("Basic test using reference system of Naoz et al. (2009)")

        particles = Tools.create_nested_multiple(3, [1.0,1.0e-3,40.0e-3],[6.0,100.0],[0.001,0.6],[0.0,65.0*np.pi/180.0],[45.0*np.pi/180.0,0.0],[0.0,0.0])
        binaries = [x for x in particles if x.is_binary==True]
        
        code = SecularMultiple()
        code.add_particles(particles)
        inner_binary = binaries[0]
        outer_binary = binaries[1]
        
        e_print = []
        INCL_print = []
        rel_INCL_print = []
        t_print = []
        
        start = time.time()
    
        t = 0.0
        N = 1000
        tend = 3.0e7

        dt = tend/float(N)
        while t<=tend:
            code.evolve_model(t)
            t+=dt
        
            if args.verbose==True:
                print( 't',t,'es',[x.e for x in binaries],'INCL_parent',inner_binary.INCL_parent)
                        
            rel_INCL_print.append(inner_binary.INCL_parent)
            e_print.append(inner_binary.e)
            INCL_print.append(inner_binary.INCL)
            t_print.append(t)
        
        if args.verbose==True:
            print('wall time',time.time()-start)
        
        t_print = np.array(t_print)
        rel_INCL_print = np.array(rel_INCL_print)
        e_print = np.array(e_print)

        assert round(e_print[-1],3) == 0.204
        assert round(rel_INCL_print[-1],3) == 1.217
        
        print("Test passed")

        code.reset()
                
        if HAS_MATPLOTLIB==True and args.plot==True:
            fig=pyplot.figure()
            plot1=fig.add_subplot(2,1,1)
            plot2=fig.add_subplot(2,1,2,yscale="log")
            plot1.plot(1.0e-6*t_print,rel_INCL_print*180.0/np.pi)
            plot2.plot(1.0e-6*t_print,1.0-e_print)
            plot2.set_xlabel("$t/\mathrm{Myr}$")
            plot1.set_ylabel("$i_\mathrm{rel}/\mathrm{deg}$")
            plot2.set_ylabel("$1-e_\mathrm{in}$")
            pyplot.show()

    def test2(self,args):
        print("Test 1PN precession in 2-body system")

        particles = Tools.create_nested_multiple(2,[1.0, 1.0], [1.0], [0.99], [0.01], [0.01], [0.01])
        binaries = [x for x in particles if x.is_binary == True]

        for b in binaries:
            b.include_pairwise_1PN_terms = True

        code = SecularMultiple()
        code.add_particles(particles)

        t = 0.0
        N=1000
        tend = 1.0e7
        dt=tend/float(N)

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)

            if args.verbose==True:
                print( 't/Myr',t,'omega',binaries[0].AP)
            t_print_array.append(t)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)
            AP_print_array.append(binaries[0].AP)
        
        t_print_array = np.array(t_print_array)
        a_print_array = np.array(a_print_array)
        e_print_array = np.array(e_print_array)
        AP_print_array = np.array(AP_print_array)
        
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C

        # Theoretical prediction #
        a = binaries[0].a
        e = binaries[0].e
        M = binaries[0].mass
        rg = CONST_G*M/(CONST_C**2)
        P = 2.0*np.pi*np.sqrt(a**3/(CONST_G*M))
        t_1PN = (1.0/3.0)*P*(1.0-e**2)*(a/rg)

        AP = 0.01 +2.0*np.pi*tend/t_1PN
        AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
        
        assert round(AP_print_array[-1],5) == round(AP,5)
        print("Test passed")

        code.reset()
                
        if HAS_MATPLOTLIB == True and args.plot==True:
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array*1.0e-6,AP_print_array, color='r',label="$\mathrm{SecularMultiple}$")
            points = np.linspace(0.0,tend*1.0e-6,N)
            AP = 0.01 +2.0*np.pi*points/(t_1PN*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',label="$\mathrm{Analytic}$")

            plot2.plot(t_print_array*1.0e-6,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array*1.0e-6,np.fabs((a-a_print_array)/a), color='r')
            plot4.plot(t_print_array*1.0e-6,np.fabs((e-e_print_array)/e), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)
            
            handles,labels = plot1.get_legend_handles_labels()
            plot1.legend(handles,labels,loc="upper left",fontsize=0.6*fontsize)

            pyplot.show()

    def test3(self,args):
        print("Test GW emission in 2-body system + collision detection")

        code = SecularMultiple()
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C

        a0 = 1.0
        e0 = 0.999
        m1 = 1.0
        m2 = 1.0
        particles = Tools.create_nested_multiple(2,[m1, m2], [a0], [e0], [0.01], [0.01], [0.01])
        bodies = [x for x in particles if x.is_binary == False]
        binaries = [x for x in particles if x.is_binary == True]
        
        for b in binaries:
            b.include_pairwise_1PN_terms = False
            b.include_pairwise_25PN_terms = True
            b.check_for_physical_collision_or_orbit_crossing = True
        
        rg = (m1+m2)*CONST_G/(CONST_C**2)
        
        for b in bodies:
            b.radius = 100.0*rg

        binary = binaries[0]

        code.add_particles(particles)

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        tend = 1e8
        N = 1000
        dt = tend/float(N)
        t = 0.0
        while (t<tend):
            t+=dt

            code.evolve_model(t)
            flag = code.flag

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

            if args.verbose==True:
                print("t",t*1e-6,'a/AU',binary.a,'e',binary.e)
            
            if flag == 2:
                if args.verbose==True:
                    print( 'root found')
                break

        assert round(t*1e-6,1) == 97.3
        print("Test passed")
        
        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure(figsize=(16,10))
            plot1 = fig.add_subplot(2,1,1)
            plot2 = fig.add_subplot(2,1,2,yscale="log")

            plot1.plot(t_print_array,e_print_array, color='r')

            plot2.plot(t_print_array,a_print_array, color='r')

            ### Peters 1964 ###
            c0 = a0*(1.0-e0**2)/( pow(e0,12.0/19.0)*pow(1.0 + (121.0/304.0)*e0**2,870.0/2299.0))
            a_an = c0*pow(e_print_array,12.0/19.0)*pow(1.0+(121.0/304.0)*e_print_array**2,870.0/2299.0)/(1.0-e_print_array**2)
            beta = (64.0/5.0)*CONST_G**3*m1*m2*(m1+m2)/(CONST_C**5)
            #T = c0**4*pow(e0,48.0/19.0)/(4.0*beta)
            T_c = a0**4/(4.0*beta)
            T = (768.0/425.0)*T_c*pow(1.0-e0**2,7.0/2.0)
            print( 'T/Myr (approx)',T*1.0e-6)
            plot2.plot(t_print_array,a_an,color='g',linestyle='dashed',linewidth=2)

            fontsize = 15
            plot1.set_ylabel("$e$",fontsize=fontsize)
            plot2.set_ylabel("$a/\mathrm{AU}$",fontsize=fontsize)
            plot2.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    
    def test4(self,args):
        print("Test tidal friction in 2-body system")
        
        code = SecularMultiple()
        code.enable_tides = True

        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        M = 0.0009546386983890755 ### Jupiter mass
        R = 40.0*0.1027922358015816*CONST_R_SUN ### 40 R_J
        m_per = 1.0
        mu = m_per*M/(m_per+M)
        a0 = 0.1
        e0 = 0.3
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        aF = a0*(1.0-e0**2)
        nF = np.sqrt( CONST_G*(M+m_per)/(aF**3) )

        particles = Tools.create_nested_multiple(2, [m_per, M], [a0], [e0], [0.01], [0.01], [0.01])
        binary = particles[2]
        
        particles[0].radius = CONST_R_SUN
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0
        particles[1].spin_vec_y = 0.0
        particles[1].spin_vec_z = 4.0e-2/day

        k_L = 0.38
        k_AM = k_L/2.0
        rg = 0.25
        tau = 0.66*second

        I = rg*M*R**2
        alpha = I/(mu*a0**2)
        T = R**3/(CONST_G*M*tau)
        t_V = 3.0*(1.0 + 1.0/k_L)*T
        if args.verbose==True:
            print( 't_V',t_V,'M',M,'R',R)

        particles[0].include_tidal_friction_terms = False
        
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = True
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-8

        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_viscous_time_scale = t_V
        particles[1].tides_gyration_radius = rg

        tD = M*aF**8/(3.0*k_L*tau*CONST_G*m_per*(M+m_per)*R**5)
        particles[2].check_for_physical_collision_or_orbit_crossing = True

        code.add_particles(particles)
        binary = code.particles[2]

        t = 0.0
        N=100
        tend = 1.0e4
        dt = tend/float(N)

        t_print_array = []
        a_print_array = []
        n_print_array = []
        e_print_array = []
        AP_print_array = []
        spin_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            if args.verbose==True:
                print( 'flag',code.flag,'t/yr',t,'a/AU',binary.a,'e',binary.e)


            t_print_array.append(t)
            a_print_array.append(binary.a)
            n_print_array.append(np.sqrt(CONST_G*(M+m_per)/(binary.a**3)))
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)
            spin_print_array.append( np.sqrt( particles[1].spin_vec_x**2 + particles[1].spin_vec_y**2 + particles[1].spin_vec_z**2) )

            bodies = particles[0:2]
            if args.verbose==True:
                for body in bodies:
                    print( 'S_x',body.spin_vec_x)
                    print( 'S_y',body.spin_vec_y)
                    print( 'S_z',body.spin_vec_z)
                print( '='*50)

        N_r = 2
        assert round(spin_print_array[-1],N_r) == round(n_print_array[-1],N_r)
        assert round(a_print_array[-1],N_r) == round(aF,N_r)
        assert len([x for x in range(len(t_print_array)) if round(aF,N_r) not in [round(a*(1.0-e**2),N_r) for a,e in zip( a_print_array,e_print_array)] ] ) == 0
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            fig = pyplot.figure(figsize=(10,10))
            fontsize=12
    
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)
            spin_print_array = np.array(spin_print_array)
            n_print_array = np.array(n_print_array)
            
            N_p = 4
            plot1 = fig.add_subplot(N_p,1,1)
            plot1.plot(t_print_array*1.0e-6,a_print_array, color='r')
            plot1.set_ylabel("$a/\mathrm{AU}$",fontsize=fontsize)

            plot2 = fig.add_subplot(N_p,1,2)
            plot2.plot(t_print_array*1.0e-6,e_print_array,color='k')
            plot2.set_ylabel("$e$",fontsize=fontsize)

            plot3 = fig.add_subplot(N_p,1,3,yscale="log")

            plot3.plot(t_print_array*1.0e-6,a_print_array*(1.0-e_print_array**2),color='k')
            plot3.axhline(y = a0*(1.0 - e0**2), color='k')
            plot3.set_ylabel("$a(1-e^2)/\mathrm{AU}$",fontsize=fontsize)

            plot4 = fig.add_subplot(N_p,1,4)
            plot4.plot(t_print_array*1.0e-6,spin_print_array/n_print_array)
            plot4.set_ylabel("$\Omega/n$",fontsize=fontsize)

            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    def test5(self,args):
        print("Test precession due to tidal bulges")

        code = SecularMultiple()
        code.enable_tides = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)

        M = 0.0009546386983890755 ### Jupiter mass
        R = 1.0*0.1027922358015816*CONST_R_SUN ### Jupiter radius ~ 0.1 R_SUN

        m_per = 1.0
        a0 = 30.0
        e0 = 0.999
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        particles = Tools.create_nested_multiple(2, [m_per, M], [a0], [e0], [0.01], [0.01], [0.01])
        binary = particles[2]
        particles[0].radius = 1.0*CONST_R_SUN
        particles[1].radius = R

        k_L = 0.41
        k_AM = k_L/2.0

        particles[1].tides_method = 0
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = True
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_gyration_radius = 0.25
        
        code.add_particles(particles)

        t = 0.0
        dt = 1.0e6
        tend = 1.0e8

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        g_dot_TB = (15.0/8.0)*n0*(8.0+12.0*e0**2+e0**4)*(m_per/M)*k_AM*pow(R/a0,5.0)/pow(1.0-e0**2,5.0)
        t_TB = 2.0*np.pi/g_dot_TB

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            if args.verbose==True:
                print( 'flag',code.flag,'t',t,'a/AU',binary.a,'e',binary.e)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

        AP = 0.01 +2.0*np.pi*tend/(t_TB)
        AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
        
        N_r = 8
        assert round(AP,N_r) == round(AP_print_array[-1],N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)
            
            fig = pyplot.figure(figsize=(10,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,len(t_print_array))
            AP = 0.01 +2.0*np.pi*points/(t_TB*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',linestyle='dotted',linewidth=2)

            plot2.plot(t_print_array,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array,np.fabs((a0-a_print_array)/a0), color='r')
            plot4.plot(t_print_array,np.fabs((e0-e_print_array)/e0), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega/\mathrm{rad}$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)
            
            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)

            pyplot.show()

    def test6(self,args):
        print("Test precession due to rotation")

        code = SecularMultiple()
        code.enable_tides = True
        CONST_G = code.CONST_G
        CONST_C = code.CONST_C
        CONST_R_SUN = code.CONST_R_SUN
        day = 1.0/365.25
        second = day/(24.0*3600.0)
        
        M = 0.0009546386983890755 ### Jupiter mass
        R = 1.5*0.1027922358015816*CONST_R_SUN ### Jupiter radius ~ 0.1 R_SUN
        m_per = 1.0
        a0 = 30.0
        e0 = 0.999
        P0 = 2.0*np.pi*np.sqrt(a0**3/(CONST_G*(M+m_per)))
        n0 = 2.0*np.pi/P0

        aF = a0*(1.0-e0**2)
        nF = np.sqrt( CONST_G*(M+m_per)/(aF**3) )

        particles = Tools.create_nested_multiple(2, [m_per, M], [a0], [e0], [1.0e-5], [1.0e-5], [1.0e-5])
        binary = particles[2]
        particles[0].radius = 1.0
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0
        particles[1].spin_vec_y = 0.0
        Omega_PS0 = n0*(33.0/10.0)*pow(a0/aF,3.0/2.0)
        particles[1].spin_vec_z = Omega_PS0

        k_L = 0.51
        k_AM = k_L/2.0
        rg = 0.25
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = True
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_gyration_radius = rg

        code.add_particles(particles)

        t = 0.0
        dt = 1.0e6
        tend = 1.0e8

        t_print_array = []
        a_print_array = []
        e_print_array = []
        AP_print_array = []

        Omega_vec = [particles[1].spin_vec_x,particles[1].spin_vec_y,particles[1].spin_vec_z]
        Omega = np.sqrt(Omega_vec[0]**2 + Omega_vec[1]**2 + Omega_vec[2]**2)
        if args.verbose==True:
            print( 'Omega/n',Omega/n0)

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            if args.verbose==True:
                print( 'flag',code.flag,'t',t,'a',binary.a,'e',binary.e)


            t_print_array.append(t*1.0e-6)
            a_print_array.append(binary.a)
            e_print_array.append(binary.e)
            AP_print_array.append(binary.AP)

        g_dot_rot = n0*(1.0 + m_per/M)*k_AM*pow(R/a0,5.0)*(Omega/n0)**2/((1.0-e0**2)**2)
        t_rot = 2.0*np.pi/g_dot_rot

        AP = 0.01 + 2.0*np.pi*tend/(t_rot)
        AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi

        N_r = 1
        assert round(AP,N_r) == round(AP_print_array[-1],N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            t_print_array = np.array(t_print_array)
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            AP_print_array = np.array(AP_print_array)

            fig = pyplot.figure(figsize=(10,10))
            plot1 = fig.add_subplot(4,1,1)
            plot2 = fig.add_subplot(4,1,2,yscale="log")
            plot3 = fig.add_subplot(4,1,3,yscale="log")
            plot4 = fig.add_subplot(4,1,4,yscale="log")

            plot1.plot(t_print_array,AP_print_array, color='r')
            points = np.linspace(0.0,tend*1.0e-6,len(t_print_array))
            AP = 0.01 +2.0*np.pi*points/(t_rot*1.0e-6)
            AP = (AP+np.pi)%(2.0*np.pi) - np.pi ### -pi < AP < pi
            plot1.plot(points,AP,color='g',linestyle='dotted',linewidth=2)

            plot2.plot(t_print_array,np.fabs( (AP - AP_print_array)/AP ), color='r')
            plot3.plot(t_print_array,np.fabs((a0-a_print_array)/a0), color='r')
            plot4.plot(t_print_array,np.fabs((e0-e_print_array)/e0), color='r')

            fontsize = 15
            plot1.set_ylabel("$\omega/\mathrm{rad}$",fontsize=fontsize)
            plot2.set_ylabel("$|(\omega_p-\omega)/\omega_p|$",fontsize=fontsize)
            plot3.set_ylabel("$|(a_0-a)/a_0|$",fontsize=fontsize)
            plot4.set_ylabel("$|(e_0-e)/e_0|$",fontsize=fontsize)

            plot4.set_xlabel("$t/\mathrm{Myr}$",fontsize=fontsize)
            pyplot.show()

    def test7(self,args):
        print("Test collision detection in 3-body system")

        code = SecularMultiple()
        
        particles = Tools.create_nested_multiple(3,[1.0, 1.2, 0.9], [1.0, 100.0], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        binaries[0].check_for_physical_collision_or_orbit_crossing = True
        R = 0.03 ### AU
        for body in bodies:
            body.radius = R

        code.add_particles(particles)

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6
        t_root = 0.0
        
        t_print_array = []
        a_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag

            if args.verbose==True:
                print("="*50)
                print("t/Myr",t*1e-6,"rp/AU",binaries[0].a*(1.0 - binaries[0].e) )
                print( 'secular_breakdown_has_occurred',binaries[0].secular_breakdown_has_occurred)
                print( 'dynamical_instability_has_occurred',binaries[0].dynamical_instability_has_occurred)
                print( 'physical_collision_or_orbit_crossing_has_occurred',binaries[0].physical_collision_or_orbit_crossing_has_occurred)
                print( 'minimum_periapse_distance_has_occurred',binaries[0].minimum_periapse_distance_has_occurred)
                print( 'RLOF_at_pericentre_has_occurred',binaries[0].RLOF_at_pericentre_has_occurred)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                t_root = code.model_time
                if args.verbose==True:
                    print( 'root found at t=',t_root)
                break
        
        N_r = 5
        assert round(a_print_array[-1]*(1.0 - e_print_array[-1]),N_r) == round(2.0*R,N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure()
            plot = fig.add_subplot(1,1,1)
            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))
            plot.axhline(y = bodies[0].radius + bodies[1].radius,color='k')
            plot.set_ylabel("$r_\mathrm{p}/\mathrm{AU}$")
            plot.set_xlabel("$t/\mathrm{Myr}$")

            pyplot.show()

    def test8(self,args):
        print("Test minimum periapse occurrence")

        code = SecularMultiple()
        
        particles = Tools.create_nested_multiple(3,[1.0, 1.2, 0.9], [1.0, 100.0], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        binaries[0].check_for_minimum_periapse_distance = True
        rp_min = 0.1
        binaries[0].check_for_minimum_periapse_distance_value = rp_min
        
        code.add_particles(particles)

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6
        t_root = 0.0
        
        t_print_array = []
        a_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag

            if args.verbose==True:
                print("="*50)
                print("t/Myr",t*1e-6,"rp/AU",binaries[0].a*(1.0 - binaries[0].e) )
                print( 'secular_breakdown_has_occurred',binaries[0].secular_breakdown_has_occurred)
                print( 'dynamical_instability_has_occurred',binaries[0].dynamical_instability_has_occurred)
                print( 'physical_collision_or_orbit_crossing_has_occurred',binaries[0].physical_collision_or_orbit_crossing_has_occurred)
                print( 'minimum_periapse_distance_has_occurred',binaries[0].minimum_periapse_distance_has_occurred)
                print( 'RLOF_at_pericentre_has_occurred',binaries[0].RLOF_at_pericentre_has_occurred)

            t_print_array.append(t*1.0e-6)
            a_print_array.append(binaries[0].a)
            e_print_array.append(binaries[0].e)

            if flag == 2:
                t_root = code.model_time
                if args.verbose==True:
                    print( 'root found at t=',t_root)
                break
        
        N_r = 5
        assert round(a_print_array[-1]*(1.0 - e_print_array[-1]),N_r) == round(rp_min,N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            a_print_array = np.array(a_print_array)
            e_print_array = np.array(e_print_array)
            
            fig = pyplot.figure()
            plot = fig.add_subplot(1,1,1)
            plot.plot(t_print_array,a_print_array*(1.0-e_print_array))
            plot.axhline(y = rp_min,color='k')
            plot.set_ylabel("$r_\mathrm{p}/\mathrm{AU}$")
            plot.set_xlabel("$t/\mathrm{Myr}$")

            pyplot.show()

    def test9(self,args):
        print("Test adiabatic mass loss")

        code = SecularMultiple()
        
        m1i = 1.0
        m2i = 1.2
        m3i = 0.9
        a1i = 1.0
        a2i = 100.0
        particles = Tools.create_nested_multiple(3,[m1i,m2i,m3i], [a1i,a2i], [0.1, 0.5], [0.01, 80.0*np.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]

        m1dot = -1.0e-7
        m2dot = -1.0e-8
        m3dot = -1.0e-9
        mdots = [m1dot,m2dot,m3dot]
        for index,body in enumerate(bodies):
            body.mass_dot = mdots[index]

        code.add_particles(particles)

        t = 0.0
        dt = 1.0e4
        tend = 1.0e6

        t_print_array = []
        a1_print_array = []
        a2_print_array = []
        e_print_array = []

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag

            if args.verbose==True:
                print( 't/Myr',t*1e-6,'m1',bodies[0].mass,'m2',bodies[1].mass,'m3',bodies[2].mass,'a/AU',binaries[0].a)

            t_print_array.append(t*1.0e-6)
            a1_print_array.append(binaries[0].a)
            a2_print_array.append(binaries[1].a)
            e_print_array.append(binaries[0].e)

        m1 = m1i + m1dot*tend
        m2 = m2i + m2dot*tend
        m3 = m3i + m3dot*tend
            
        a1 = a1i*(m1i+m2i)/(m1+m2)
        a2 = a2i*(m1i+m2i+m3i)/(m1+m2+m3)

        N_r = 4
        assert round(a1,N_r) == round(a1_print_array[-1],N_r)
        assert round(a2,N_r) == round(a2_print_array[-1],N_r)
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB == True and args.plot==True:
            a1_print_array = np.array(a1_print_array)
            a2_print_array = np.array(a2_print_array)
            e_print_array = np.array(e_print_array)

            t = pow(10.0,np.linspace(4.0,6.0,100))
            m1 = m1i + m1dot*t
            m2 = m2i + m2dot*t
            m3 = m3i + m3dot*t
            
            a1 = a1i*(m1i+m2i)/(m1+m2)
            a2 = a2i*(m1i+m2i+m3i)/(m1+m2+m3)

            fig = pyplot.figure(figsize=(10,8))
            plot1 = fig.add_subplot(2,1,1,yscale="log")
            plot2 = fig.add_subplot(2,1,2,yscale="log")
            plot1.plot(t_print_array,a1_print_array,color='k',linestyle='solid',label='$a_1$',zorder=10)
            plot2.plot(t_print_array,a2_print_array,color='k',linestyle='solid',label='$a_2$',zorder=10)
            
            plot1.plot(t*1.0e-6,a1,color='g',linestyle='dashed',linewidth=3,label='$a_1\,(\mathrm{an})$')
            plot2.plot(t*1.0e-6,a2,color='g',linestyle='dashed',linewidth=3,label='$a_2\,(\mathrm{an})$')

            handles,labels = plot1.get_legend_handles_labels()
            plot1.legend(handles,labels,loc="lower left",fontsize=18)

            plot1.set_ylabel("$a_1$")
            plot2.set_ylabel("$a_2$")
            plot2.set_xlabel("$t/\mathrm{Myr}$")

            pyplot.show()


    def test10(self,args):
        print("Test hybrid integration")

        m1=1.0
        m2=1.0e-6
        m3=1.0
        e1=0
        e2=0.4
        a1=1.0
        a2=10.0

        i1=0.01
        i2=65.0*np.pi/180.0
        AP1=0
        AP2=0
        LAN1=0
        LAN2=0

        particles = Tools.create_nested_multiple(3, [m1,m2,m3],[a1,a2],[e1,e2],[i1,i2],[AP1,AP2],[LAN1,LAN2])
        
        bodies = [x for x in particles if x.is_binary==False]
        binaries = [x for x in particles if x.is_binary==True]
        N_binaries = len(binaries)
        
        code = SecularMultiple()
        code.add_particles(particles)

        integration_methods = [0,1]
        KS_use_V = [True,True]
        terms = [False,True,True,True,True,True]
        
        binaries[0].integration_method = integration_methods[0]
        binaries[1].integration_method = integration_methods[1]

        binaries[0].KS_use_perturbing_potential = KS_use_V[0]
        binaries[1].KS_use_perturbing_potential = KS_use_V[1]

        CONST_G = code.CONST_G
        P1=2.0*np.pi*np.sqrt(a1**3/(CONST_G*(m1+m2)))
        P2=2.0*np.pi*np.sqrt(a2**3/(CONST_G*(m1+m2+m3)))
        P_LK12 = (P2**2/P1)*((m1+m2+m3)/m3)*pow(1.0-e2**2,3.0/2.0)
        if args.verbose==True:
            print("P_orbs/yr",P1,P2)
            print("P_LK/yr",P_LK12)

        code.include_double_averaging_corrections = terms[0]
        code.include_quadrupole_order_terms = terms[1]
        code.include_octupole_order_binary_pair_terms = terms[2]
        code.include_octupole_order_binary_triplet_terms = terms[3]
        code.include_hexadecupole_order_binary_pair_terms = terms[4]
        code.include_dotriacontupole_order_binary_pair_terms = terms[5]
        
        a_print = [[] for x in range(N_binaries)]
        e_print = [[] for x in range(N_binaries)]
        t_print = []
        rel_INCL_print = []
        

        code.enable_root_finding = True
        binaries[0].check_for_physical_collision_or_orbit_crossing=True
        bodies[0].radius=1.0e-5
        bodies[1].radius=1.0e-5
        
        import time
        
        t = 0.0
        N = 1000
        tend = 1.e3
        dt = tend/float(N)
                
        start = time.time()
        while t<tend:
            t+=dt
            code.evolve_model(t)
            
            if args.verbose==True:
                print( 't/Myr',t*1e-6,'es',[x.e for x in binaries],'smas',[x.a for x in binaries])

            rel_INCL_print.append(binaries[0].INCL_parent)                        
            for i in range(N_binaries):
                a_print[i].append(binaries[i].a)
                e_print[i].append(binaries[i].e)
            t_print.append(t)
            
            flag=code.flag
            if flag==2:
                print("RF")
                break
        wall_time = time.time()-start

        t_print = np.array(t_print)
        rel_INCL_print = np.array(rel_INCL_print)
        for i in range(N_binaries):
            a_print[i] = np.array(a_print[i])
            e_print[i] = np.array(e_print[i])

        N_r = 4
        assert round(e_print[0][-1],N_r) == 0.5274
        print("Test passed")

        code.reset()
        
        if HAS_MATPLOTLIB==True and args.plot==True:
            if 1==0:
                pyplot.rc('text',usetex=True)
                pyplot.rc('legend',fancybox=True)

            
            fig=pyplot.figure(figsize=(10,8))
            plot1=fig.add_subplot(2,1,1)
            plot2=fig.add_subplot(2,1,2,yscale="log")

            plot1.plot(1.0e-6*t_print,rel_INCL_print*180.0/np.pi)
            plot2.plot(1.0e-6*t_print,1-e_print[0])

            fontsize=18

            plot1.set_ylabel("$i_\mathrm{rel}\,(\mathrm{deg})$",fontsize=fontsize)
            plot2.set_ylabel("$1-e$",fontsize=fontsize)

            plot2.set_xlabel("$t/\mathrm{Myr}$")
            
            pyplot.show()

    
if __name__ == '__main__':
    args = parse_arguments()
    
    N_tests = 10
    if args.test==0:
        tests = range(1,N_tests+1)
    else:
        tests = [args.test]

    t=test_secularmultiple()
    for i in tests:
        print( 'Running test number',i,'; verbose =',args.verbose,'plot =',args.plot)
        function = getattr(t, 'test%s'%i)
        function(args)
