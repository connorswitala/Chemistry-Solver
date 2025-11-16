# Chemical Equilibrium Solver

This repository is for solving chemical equilibrium thermodynamics. It is currently under development to include Helmholtz and Gibbs energy minimization. 

# Programs

In the [bin](./bin/) directory, there are three programs ready to be run after you build. The first is an example program that has everything defined for you to run it and see the minimization results. You can toy with options to understand how to use correct `enum class` types, or `ConstraintTypes::` (both defined below).

The second program is the command-line driven code that lets you choose options and set mixtures for the chemical equilibrium process. The program has a 'help' command that will tell you what you need to run, but an brief overview is given here. When you start the program, it will give you some options right off the bat. In order to run minimization, all fields need to filled in that appear when you type `show`. Two fields are already initialized by default. 

The `--mode` option tells the program which mode it is running in. By inputting `--mode=sweep` you can change the programs from a standalone minimization given a single value for T and P, to a verson that sweeps through an user-input min and max temperature value and then plots the results to a user-named file.

The `--constraint` option tells the programs with minimization procedure to use. TP, TV, UV, CFD are the current ones available. It needs to be in caps when input as well. For example`--constraint=TP`.

The `--species` option tells the program which species you want to include in the minimization process. You can list your own with comma (and space) separated variables as such: `--species=N2, O2, NO, N, O`. If you choose to build species this way, you will need to specify a couple more things like `--elements` and `-Y`. If instead you use an already created mixture (available ones are air5, air7, air11, air13, mars8), you can do so with `--species=air11`, for example. This will fill in the `--elements` and `-Y` flags for you as well. You can always double check what needs to be filled in with the `show` command.

If you do need to specifiy elements and mass fractions, you can do so with `--elements=N, O` and also `-Y=0.7572, 0.2428`. Just make sure that the order that you input the mass fractions lines up with the elemental ordering.

When the configuration input has been filled in, `run` will run the program for you. You will still need to follow command line prompts for specifying filenames, sweeping values, etc.

Finally, the third program is created for timing the speed of the minimization processes and is really more for developing than for getting values.

# Codebase Information

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

- `ConstraintType::TP` Minimization holds temperature (T) and pressure (P) constant.
- `ConstraintType::HP` Minimization holds enthalpy (H) and pressure (P) constant.
- `ConstraintType::SP` Minimization holds entropy (S) and pressure (P) constant.
- `ConstraintType::TV` Minimization holds temperature (T) and volume (V) constant.
- `ConstraintType::UV` Minimization holds internal energy (U) and volume (V) constant.
- `ConstraintType::SV` Minimization holds entropy (S) and volume (V) constant.


The `CESolver` constructor function automatically chooses which energy to minimize when it receives the constraint type. The constructor function creates a pointer to the correct minimization function and is called with `CESolver CE(gas, constraint)`. Here, `gas` is the gas mix created above, and `constraint` is the `ConstraintType` you want to use. The minimization function pointer is used in the public member function of `CESolver` named `compute_equilibrium()`. This function takes in two arguments that are alligned with the contraint type, e.g. if you have chosen `ConstraintType::TP`, then calling `compute_equilibrium(1000, 101325)` will assume that temperature T = 1000 K, and pressure P = 101325 Pa. If you set 'constraint' = `ConstraintType::UV`, then `compute_equilibrium(1e6, 0.1)` would set internal energy (U) to 1e6 J/kg , and volume (V) to 0.1 m^3/kg.

Another constraint type is used when this code is utilized for plugging in to CFD code. This is `ConstraintType::CFD` which does not create a pointer to a certain function. It forced you to call the CFD version `CFD_equilibrium(e, rho)` which minimizes Helmholtz energy hold internal energy and density constant.

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

This function will display every thermodynamic quantity you may need in `mix gas`. The quantities that are printed in this function are denoted in the structure definition [here](./libraries/includes/thermoObjs.h).

## Sample Code

Code is already setup to test the functionality of this repository. To check it out, open [minimize.cpp](./source/).

## Documentation

There is scientific documentation on the formulation of the equations being used and in what functions they show up. It is written in Latex and can be found [here](./documents/EnergyMinimization.pdf).
