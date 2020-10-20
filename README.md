# SecularMultiple
    
A code to compute the secular (orbit-averaged) gravitational dynamics of hierarchical multiple systems 
composed of nested binary orbits (simplex-type systems) with any configuration and any number of bodies.
A particle can repesent a binary (`is_binary = True`) or a body (`is_binary = False`).
The structure of the system is determined by linking to other particles with the attributes `child1` and `child2`.
Tidal interactions and relativistic corrections are included in an ad hoc fashion
(tides: treating the companion as a single body, even if it is not; relativistic terms:
only including binary-binary interactions). Hybrid integration (averaged, or direct integration) is also supported.
    
Includes routines for external perturbations (flybys & supernovae).

If you use this code for work in scientific publications, please cite:
https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2827H (the original paper)
https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.4139H (updates with external perturbations)
https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.5492H (updates with suborbital effects).

A C++ compiler is required, as well as Python (2/3) for the Python interface. Make sure to first compile the code using `make`. It will compile using the system's default C++ compiler (if you want to change this, you could modify the Makefile). 

The script `test_secularmultiple.py` can be used to test the
installation. See `examples.py` for some examples.

Adrian Hamers, February 2020
