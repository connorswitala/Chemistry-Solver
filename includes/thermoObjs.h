#ifndef THERMOOBJS_H
#define THERMOOBJS_H

#include "common.h"

constexpr int NCOEF = 9; // Number of coefficients for NASA polynomials
using namespace std;

enum class EnergyType {Helmholtz, Gibbs};
enum class ConstraintType {TP, HP, SP, TV, UV, SV};
enum class GasType {AIR5, AIR7, AIR11, AIR11_AR, AIR13, CREATE};

struct SpeciesInfo {
    string name;        // Name of species
    double mw;          // Molecular weight of species
    Vector poly;        // NASA polynomial coefficients
    int q;              // Charge of species
    double hf;
    double href;
};

struct mix {

    string name;
    // ======================: Displayed properties in this box :=======================
    //                                                                                  |
    vector<SpeciesInfo> species;            // Vector of Species                        |
    int NS, NE;                             // Sizes                                    |
    double R, gamma, cp, cv, MW, Pr, k, D;  // Thermodynamic properties of mixture.     |
    double e, rho, T, p, V;                 // Thermodynamic state variables            |
    Vector Y, X;                            // Mass and mole fractions.                 |
    //                                                                                  |
    // =================================================================================

    Vector H0_RT, S0_R, mu0_RT, CP0_R, U0_RT;  // NASA polynomial
    double up, hp, sp, uo, ho, so;  // Used for specification of enthalpy, internal energy, and entropy.
    double N_tot;
    double e_ref;
    Vector X0;                      // Initial Moles of air.
    Vector a, b;                    // Stoichiometric coefficients / number of moles of element i 
    Vector N, mu_RT;                // Current solution for number of moles as well as chemical potential.
    int J_SIZE;                     // Size of solution vector.
    bool HAS_IONS, NEEDS_T;         // Booleans used in solver, stores here for passing to function.
};

#endif