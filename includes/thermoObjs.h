#ifndef THERMOOBJS_H
#define THERMOOBJS_H

#include "common.h"

constexpr int NCOEF = 9; // Number of coefficients for NASA polynomials
using namespace std;

enum class EnergyType {Helmholtz, Gibbs};
enum class ConstraintType {TP, HP, SP, TV, UV, SV};
enum class GasType {AIR5, AIR11, AIR11_AR, AIR13, CREATE};

struct SpeciesInfo {
    string name;
    double mw;
    Vector poly;
    double theta_v;
    int q;
    double R;
};

struct mix {

    vector<SpeciesInfo> species;

    int NS, NE, NI, N_RE, N_PROD;       // Sizes
    Vector H0_RT, S0_R, mu0_RT, CP0_R;  // NASA polynomial
    Vector N, mu_RT;                           // Current solution for number of moles.

    int J_SIZE; // Size of solution vector.

    Vector hf;
    double R, gamma, cp, cv, MW, Pr, k, D;
    double e, rho, T, p, V;
    Vector Y;
    Vector X, X0;
    Vector a, b;
    double initial_moles, N_tot;
    vector<int> diatomic_list;
    vector<int> mono_list;
    bool HAS_IONS;
    Vector reactant_idx;
    Vector product_idx;
    Vector reactions;
};

struct Context {
    bool HAS_IONS, NEEDS_T; 
    EnergyType energy_type;
    ConstraintType constraint_type; 
};

#endif