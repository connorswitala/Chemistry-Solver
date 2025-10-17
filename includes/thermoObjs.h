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
    Vector H0, S0, mu0, Cp0;
    double R, gamma, cp, cv, gam;
    double e, rho, T, p;
    Vector Y;
    Vector X, X0;
    int N_SP, N_EL, N_ION;
    Vector a, b;
    double initial_moles;
};

#endif