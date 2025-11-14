#ifndef CESolver_H
#define CESolver_H

#include "../mixes/mixes.h"
#include "../includes/functions.h"
#include "../includes/math.h"
#include "../includes/readNASA.h"


class CESolver {
        
    public:

    mix& gas; // Main object.

    // Main constructor functions
    CESolver(mix& gas_in, ConstraintType& contrainttype); 

    // Display results.
    void print_properties(); 

    // Main equilibrium function
    void compute_equilibrium(double& a, double& b);

    // Main CFD function
    void CFD_equilibrium(double T, double V);
    
    private:

    // Used for compute_equilibrium function
    using FuncPtr = void (CESolver::*)(double&, double&);
    FuncPtr equilibrium = nullptr;

    int T_flag;     // Flags which set of NASA coefficients to use.
    int J_SIZE;     // Size of solution vector in Newton Method.
    int NS, NE, NI; // Equation sizes

    Vector J, F, DELTA;   // Newton method vectors

    inline void compute_equilibrium_TV(double& T, double& V);     // All in one Equilibirum solver function
    inline void compute_equilibrium_UV(double& U, double& V);     // All in one Equilibirum solver function
    inline void compute_equilibrium_SV(double& S, double& V);     // All in one Equilibirum solver function
    inline void compute_equilibrium_TP(double& T, double& P);     // All in one Equilibirum solver function
    inline void compute_equilibrium_HP(double& H, double& P);     // All in one Equilibirum solver function
    inline void compute_equilibrium_SP(double& S, double& P);     // All in one Equilibirum solver function


    inline void findTRange();               // Finds Temperature range to use proper NASA coeff.
    inline array<double, 7> temp_base();    // Calculates Temperature functions for NASA poly.
    inline void NASA_fits();                // Calculates H0, S0, H0, and MU0.

    inline void compute_mixture_properties();

    // Convergence functions
    inline bool check_convergence(double* dlnj, double& dln);
    inline double compute_damping(double* dlnj, double& dln, double& dlnT);
};

#endif

