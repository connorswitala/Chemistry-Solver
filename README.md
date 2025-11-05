# Chemical Equilibrium Solver

This repo is for solving for chemical equilibrium thermodynamics. It is currently under development to include Helmholtz and Gibbs energy minimization. 

## How to run.



#### Composition

Currently, there are four air compositions available and built into the code itself. These are accesible through the the use of the enum class `GasType::`. The four compositions are:

- `GasType::AIR5` Contains N2, O2, NO, N, and O.
- `GasType::AIR11` Contains N2, O2, NO, N, O, N2+, O2+, N+, O+, e-.
- `GasType::Air11_AR` Contains N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-.
- `GasType::Air13` Contains N2, O2, NO, N, O, Ar, Ar+, N2+, O2+, N+, O+, e-.

In a future update, `GasUpdate::CREATE` will allow the user to specify the composition that is needed for their case. To create a gas mix, the following two lines are used:

```
GasType g           = GasType::AIR5; 
mix gas = common_air::create_air_mix(g);
```

This creates an instance of mix gas for you to access. It is also passed into the equilibrium solver where its values are used / edited during the minimization process.

#### Minimization Procedure
This code is set up to minimize either Gibbs or Helmholtz energies using the `CESolver` class. The method used is specified through the used of the enum class `ConstraintType::`. The usual constraints are:

- `ConstrainType::TP` Minimization holds temperature (T) and pressure (P) constant.
- `ConstrainType::HP` Minimization holds enthalpy (H) and pressure (P) constant.
- `ConstrainType::SP` Minimization holds entropy (S) and pressure (P) constant.
- `ConstrainType::TV` Minimization holds temperature (T) and volume (V) constant.
- `ConstrainType::UV` Minimization holds internal energy (U) and volume (V) constant.
- `ConstrainType::SV` Minimization holds entropy (S) and volume (V) constant.


The `CESolver` constructor function automatically chooses which energy to minimize when it receives the constraint type. The constructor function creates a pointer to the correct
minimization function and is called with `CESolver CE(gas, c)` Here, `gas` is the gas mix created aboce, and `c` is the `ConstraintType`. The minimization function pointer is used in the public member function of `CESolver` named `compute_equilibrium`.
This function takes in two arguments that are alligned with the contraint type, e.g. if you have chosen `ConstraintType::TP`, then calling `compute_equilibrium(1000, 101325)` will assume that temperature T = 1000 K, 
and pressure P = 101325 Pa. If you set c = `ConstraintType::UV`, then `compute_equilibrium(1000, 101325)` would set internal energy (U) to 1000, and volume (V) to 101325.

#### Access to results.

There are two ways to access the results of the minimization. Firstly, you can access each variable of the mix instance that you created. Assumed that for this README, it is created as `mix gas`. All results are stored in this variable, so you can use or print all of its member variables. Everything that is available to you is found in the header file "thermoObjs.h". You can access the variables you want by use of `gas.<variable name>` If you want temperature printed, you would write

`cout << gas.T << endl;`

Or for printing mass fractions of each species:

```
for (int j = 0; j < gas.NS; ++j) {
    cout << gas.Y[j] << endl;
}
```