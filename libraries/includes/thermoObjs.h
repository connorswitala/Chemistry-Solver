#ifndef THERMOOBJS_H
#define THERMOOBJS_H

#include "common.h"
using namespace std;

constexpr int NCOEF = 9; // Number of coefficients for NASA polynomials

enum class ConstraintType {TP, HP, SP, TV, UV, SV, CFD};            // Enum class for minimization constraint
enum class GasType {AIR5, AIR7, AIR11, AIR11_AR, AIR13, MARS8};     // Enum class for common gas mixes


struct SpeciesInfo {
    string name;        // Name of species
    double mw;          // Molecular weight of species
    Vector poly;        // NASA polynomial coefficients
    int q;              // Charge of species
    double hf;          // Enthalpy of formation
    double href;        // Reference enthalpy 
};

struct mix {

    string name;
    // ======================: Displayed properties in this box :=======================
    //       
    vector<string> speciesnames;            // For access in command line version                                                                           
    vector<SpeciesInfo> species;            // Vector of SpeciesInfo                         
    vector<string> elements;                // Vector of element names                 
    int NS, NE;                             // NS = number of species, NE = number of atomic elements                                   
    double R, gamma, cp, cv, MW, c;         // Thermodynamic properties of mixture.    
    double e, rho, T, p, V;                 // Thermodynamic state variables           
    Vector Y, X;                            // Mass and mole fractions.                
    double dpdr, dtdr;                      // Dericatives wrt rho holding e const     
    //                                                                                 
    // =================================================================================

    // ==================== NASA data and solver related information ===================
    Vector elemental_mw;                        // Molecular weight of elements
    Vector H0_RT, S0_R, mu0_RT, CP0_R, U0_RT;   // NASA polynomial curve fit data
    double up, hp, sp, uo, ho, so;              // Used for specification of enthalpy, internal energy, and entropy.
    double N_tot;                               // Total N
    double e_ref;                               // Sum of weighted reference enthalpies H(298.15) - H(0)
    Vector X0;                                  // Initial mass fraction of elements.
    Vector a, b;                                // Stoichiometric coefficients / number of moles of element i 
    Vector N, mu_RT;                            // Current solution for number of moles as well as chemical potential.
    int J_SIZE;                                 // Size of solution vector.
    bool HAS_IONS, NEEDS_T, COND;               // Booleans used in solver, stores here for passing to function.
};

inline ConstraintType constraintFromString(string& s) {

    ConstraintType c;

    if (s == "TP")
        return ConstraintType::TP;
    else if (s == "UV")
        return ConstraintType::UV;
    else if (s == "TV")
        return ConstraintType::TV;
    else if (s == "CFD")
        return ConstraintType::CFD;    
}


#endif