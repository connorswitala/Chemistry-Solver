# Chemical Equilibrium Solver

This repo is for solving for chemical equilibrium thermodynamics. It is currently under development to include Helmholtz and Gibbs energy minimization. 

## How to run.



#### Composition

Currently, there are four air compositions available and built into the code itself. These are accesible through the the use of the enum class `GasType::`. The four compositions are:

- `GasType::AIR5` Contains N2, O2, NO, N, and O.
- `GasType::AIR11` Contains N2, O2, NO, N, O, N2+, O2+, N+, O+, e-.
- `GasType::Air11_AR` Contains N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-.
- `GasType::Air13` Contains N2, O2, NO, N, O, Ar, Ar+, N2+, O2+, N+, O+, e-.

In a future update, `GasUpdate::CREATE` will allow the user to specify the composition that is needed for their case.

#### Minimization Procedure
This code is set up to minimize either Gibbs or Helmholtz energies using the `CESolver` class. The method used is specified through the used of the enum class `ConstraintType::`. The usual constraints are:

- `ConstrainType::TP` Minimization holds temperature (T) and pressure (P) constant.
- `ConstrainType::HP` Minimization holds enthalpy (H) and pressure (P) constant.
- `ConstrainType::SP` Minimization holds entropy (S) and pressure (P) constant.
- `ConstrainType::TV` Minimization holds temperature (T) and volume (V) constant.
- `ConstrainType::UV` Minimization holds internal energy (U) and volume (V) constant.
- `ConstrainType::SV` Minimization holds entropy (S) and volume (V) constant.


The `CESolver` constructor function automatically chooses which energy to minimize when it receives the constraint type. The constructor function creates a pointer to the correct
minimization function and is called with `CESolver CE(g, c)` Here, `g` is the `GasType`, and `c` is the `ConstraintType`. This pointer is used in the public member function of `CESolver` named `compute_equilibrium`.
This function takes in two arguments that are alligned with the contraint type, e.g. if you have chosen `ConstraintType::TP`, then calling `compute_equilibrium(1000, 101325)` will assume that temperature T = 1000 K, 
and pressure P = 101325 Pa. If you set c = `ConstraintType::UV`, then `compute_equilibrium(1000, 101325)` would set internal energy (U) to 1000, and volume (V) to 101315.