Data and code to reproduce the results of 
Richter and Mastropieri, Muscle force distribution with co-activation: A game-theoretical optimality criterion
Journal of Biomechanics Open
https://doi.org/10.1016/j.jbmo.2026.100005

Hanz Richter, 2026
Center for Human-Machine Systems,  Cleveland State University
h.richter@csuohio.edu

The comparison with the method of Grimmelsmann et. al. is found in the Comparison folder, see README file there.

Scripts in top folder:

TorqueEstimation.m: 
Runs the muscle torques estimation per the extended quadratic or game-theoretic criterion, run 

fit_h.m:
Uses optimization (requires fmincon, Optimization Toolbox) to fit Hill-Katz muscle model parameters independently
to 3 data sets.

BetaSens.m:
Sweeps weight beta to show the effect on the peak and mean forces and co-contracting force.

 