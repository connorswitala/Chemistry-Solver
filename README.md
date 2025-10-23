# Chemical Equilibrium Solver

This repo is for solving for chemical equilibrium thermodynamics. It currently has everything necessary to compute mass/molar concentrations of air as well as thermodynamic properties at a given density and internal energy. 

## Executables

Three executables can be run. They are: ./gibbs_tp [T] [p] [mix_type], ./gibbs_re [rho] [e] [mix_type], and 

[rho] is the density, [e] is the internal energy, and [mix_type] is a string that says which mixture to use. Currently, there are 5 different air mixtures available:

- air5 - This uses a 5-species air mixture that contains N2, O2, NO, N, and O.
- air11 - This uses an 11-species air mixture that contains N2, O2, NO, N, O, N2+, O2+, NO+, N+, O+, and e-.
- air11_Ar - This uses an 11-species air mixture that contains N2, O2, NO, N, O, Ar, Ar+, NO+, N+, O+, and e-.
- air13 - This uses a 13-species air mixture that contains N2, O2, NO, N, O, Ar, Ar+, N2+, O2+, NO+, N+, O+, and e-.
- perf - This just returns perfect gas thermodynamics quantities for air (e.g. R = 287.0, gamma = 1.4, etc.) No minimization procedure is used.

An example of running this code looks like:

./gibbs 0.1 2e7 air13 

This will find the thermodynamic state of 13-species air in chemical equilibrium.

## Main plugin
