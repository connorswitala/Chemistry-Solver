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

inline void LUSolve(double* A, const double* B, double* X, int n, int m) {
    // Copy A because LU is in-place and we don't want to destroy the original
    std::vector<double> LU(A, A + n * n);

    // In-place LU Decomposition (no pivoting)
    for (int k = 0; k < n; ++k) {
        for (int i = k + 1; i < n; ++i) {
            LU[i * n + k] /= LU[k * n + k];
            for (int j = k + 1; j < n; ++j) {
                LU[i * n + j] -= LU[i * n + k] * LU[k * n + j];
            }
        }
    }

    std::vector<double> y(n, 0.0);

    for (int col = 0; col < m; ++col) {
        // Forward Substitution: Ly = b
        for (int i = 0; i < n; ++i) {
            y[i] = B[i * m + col];  // row-major access
            for (int j = 0; j < i; ++j)
                y[i] -= LU[i * n + j] * y[j];
        }

        // Backward Substitution: Ux = y
        for (int i = n - 1; i >= 0; --i) {
            double sum = y[i];
            for (int j = i + 1; j < n; ++j)
                sum -= LU[i * n + j] * X[j * m + col];  // row-major
            X[i * m + col] = sum / LU[i * n + i];
        }
    }
}

inline double norm(double* v1, double* v2, int n) {
        double result = 0.0;
        for (int i = 0; i < n; ++i) {
                result += fabs(v1[i] - v2[i]) * fabs(v1[i] - v2[i]);
        }
        return sqrt(result); 
}



#endif

