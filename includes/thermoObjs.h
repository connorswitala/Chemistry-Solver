#ifndef THERMOOBJS_H
#define THERMOOBJS_H

#include "math.h"

constexpr int NCOEF = 9; // Number of coefficients for NASA polynomials

struct SpeciesInfo {

    string name;
    double mw;
    Vector poly;
    double theta_v;
    int q;
};

struct mix {

    vector<SpeciesInfo> species;
    Vector H0, S0, mu0;
    double R, gamma, cp, cv;
    double e, rho, T, p;
    Vector Y;
    Vector X;

};

#endif