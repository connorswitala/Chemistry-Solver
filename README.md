# Chemical Equilibrium Solver

This repository is for solving chemical equilibrium thermodynamics. It is currently under development to include Helmholtz and Gibbs energy minimization. 

## Completed Functions

Currently, the minimization methods that are completed are:

- Gibbs minimization holding T, P constant
- Helmholtz minimization holding T, V constant
- Helmholtz minimization holding U, V constant

## How to run

### Composition

Currently, there are five air mixtures and one Martian atmosphere available for easy access. These are accesible through the the use of the enum class `GasType::` and are created using the function `common_mixture()`. The five compositions included are:

- `GasType::AIR5` Contains N2, O2, NO, N, and O.
- `GasType::AIR7` Contains N2, O2, NO, N, O, NO+, e-.
- `GasType::AIR11` Contains N2, O2, NO, N, O, N2+, O2+, N+, O+, e-.
- `GasType::Air11_AR` Contains N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-.
- `GasType::Air13` Contains N2, O2, NO, N, O, Ar, Ar+, N2+, O2+, N+, O+, e-.
- `GasType::MARS8` Contains CO2, N2, O2, CO, O, C, NO, N.

Here is an example of creating a common air mixture:

```
GasType g = GasType::AIR5; 
mix gas = common_mixture(g);
```

This creates an instance of type `mix` named `gas` for you to access. It is also passed into the equilibrium solver where its values are used / manipulated during the minimization process. 

If you want to create your own mixture, you can define three vectors. The first is a vector of strings that contains the species names, the second is a vector of strings that contains a names of the elemeNts, and the third is a vector that contains the initial mass fractions of the elements. These are then passed into the `create_mixture()` function which returns the `mix` object. For example:

```
vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
vector<string> elements = {"C", "O", "N"};
vector<double> initial = {0.264, 0.7186, 0.0174}; 

mix gas = create_mixture(species, elements, initial);
```


### Minimization Procedure
This code is set up to minimize either Gibbs or Helmholtz energies using the `CESolver` class. The method used is specified through the used of the enum class `ConstraintType::`. The usual constraints are:

- `ConstrainType::TP` Minimization holds temperature (T) and pressure (P) constant.
- `ConstrainType::HP` Minimization holds enthalpy (H) and pressure (P) constant.
- `ConstrainType::SP` Minimization holds entropy (S) and pressure (P) constant.
- `ConstrainType::TV` Minimization holds temperature (T) and volume (V) constant.
- `ConstrainType::UV` Minimization holds internal energy (U) and volume (V) constant.
- `ConstrainType::SV` Minimization holds entropy (S) and volume (V) constant.


The `CESolver` constructor function automatically chooses which energy to minimize when it receives the constraint type. The constructor function creates a pointer to the correct
minimization function and is called with `CESolver CE(gas, constraint)`. Here, `gas` is the gas mix created above, and `constraint` is the `ConstraintType` you want to use. The minimization function pointer is used in the public member function of `CESolver` named `compute_equilibrium`.
This function takes in two arguments that are alligned with the contraint type, e.g. if you have chosen `ConstraintType::TP`, then calling `compute_equilibrium(1000, 101325)` will assume that temperature T = 1000 K, 
and pressure P = 101325 Pa. If you set 'constraint' = `ConstraintType::UV`, then `compute_equilibrium(1e6, 0.1)` would set internal energy (U) to 1e6 J/kg , and volume (V) to 0.1 m^3/kg.

### Access to results.

There are two ways to access the results of the minimization.

Firstly, you can access each variable of the `mix` instance that you created in your code. Assume that, for this README, it is created as `mix gas` as defined above. All results are stored in this variable, so you can use or print all of its member variables. Everything that is available to you is found in the header file 'thermoObjs.h'. You can access the variables you want by use of `gas.<variable name>`. If you want temperature printed, you would write

```
cout << gas.T << endl;
```

Or for printing mass fractions of each species:

```
for (int j = 0; j < gas.NS; ++j) {
    cout << gas.Y[j] << endl;
}
```

Secondly, if you only care about the results without putting it in other code, you can print the results to the terminal with:

```
print_properties(gas);
```

This function will display every thermodynamic quantity you may need in `mix gas`. The quantities that are printed in this function are denoted in the structure definition in 'thermoObjs.h'

## Sample Code

Code is already setup to test the functionality of this repository. To check it out, open [minimize.cpp](./source/).

## Documentation

There is scientific documentation on the formulation of the equations being used and in what functions they show up. It is written in Latex and can be found [here](./documents/EnergyMinimization.pdf).
