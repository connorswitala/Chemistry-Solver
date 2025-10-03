#ifndef MATH_H
#define MATH_H

#include "common.h"




inline void addVectors(const double* v1, const double* v2,  double* r, const int n) {
    for (int i = 0; i < n; ++i) {
        r[i] = v1[i] + v2[i];
    }
}

inline void subVectors(const double* v1, const double* v2, Vector& r, const int n) {
    for (int i = 0; i < n; ++i) {
        r[i] = v1[i] - v2[i];
    }
}

inline void outerProd(const double* v1, const double* v2, double* r, const int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            r[i * n + j] = v1[i] * v2[j];
        }
    }
}

inline double dotProd(const double* v1, const double* v2, const int n) {
    double r = 0.0;
    for (int i = 0; i < n; ++i) {
        r += v1[i] * v2[i];
    }
    return r;
}

#endif

