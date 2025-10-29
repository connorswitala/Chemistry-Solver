#ifndef THERMOOBJS_H
#define THERMOOBJS_H

#include "common.h"

constexpr int NCOEF = 9; // Number of coefficients for NASA polynomials
using namespace std;

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

    int N_SP, N_EL, N_ION, N_RE, N_PROD;

    Vector H0, S0, mu0, hf;
    double R, gamma, cp, cv, MW, Pr, mu, k, D;
    double e, rho, T, p;
    Vector Y;
    Vector X, X0;
    Vector a, b;
    double initial_moles, N_tot;
    vector<int> diatomic_list;
    vector<int> mono_list;
    bool perf_flag, HAS_IONS;
    Vector reactant_idx;
    Vector product_idx;
    Vector reactions;
};

#endif