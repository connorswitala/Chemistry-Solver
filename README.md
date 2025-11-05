# Chemical Equilibrium Solver

This repo is for solving for chemical equilibrium thermodynamics. It is currently under development to include Helmholtz and Gibbs energy minimization. 

## How to run.

Currently, there are four air compositions available and built into the code itself. These are accesible through the the use of the enum class `GasType::`. The four compositions are:

- `GasType::AIR5` => Contains N2, O2, NO, N, and O.
- `GasType::AIR11` => Contains N2, O2, NO, N, O, N2+, O2+, N+, O+, e-.
- `GasType::Air11_AR` => Contains N2, O2, NO, N, O, Ar, Ar+, N+, O+, e-.
- `GasType::Air13` => Contains N2, O2, NO, N, O, Ar, Ar+, N2+, O2+, N+, O+, e-.

