#ifndef CESolver_H
#define CESolver_H

#include "../common_mixes/commonMixes.h"
#include "../includes/functions.h"
#include "../includes/math.h"


class CESolver {
        
    public:

    mix& gas;        // Main object.

    CESolver(mix& gas_in, ConstraintType& contrainttype);                       // Constructor.

    void print_properties();                      // Display results.

    void plot_concentrations_for_T_range();

    void compute_equilibrium(double& a, double& b);

    void CFD_equilibrium(double T, double V);
    
    private:

    using FuncPtr = void (CESolver::*)(double&, double&);
    FuncPtr equilibrium = nullptr;

    int T_flag; // Flags which set of NASA coefficients to use.

    int J_SIZE; // Size of solution vector in Newton Method.

    int NS, NE, NI;

    Vector J, F, DELTA;   // Old solution vector of molar concentrations.
    Vector XI, J_STAR;
    
    int XI_SIZE, XI_ROWS;

    inline void compute_equilibrium_TV(double& T, double& V);     // All in one Equilibirum solver function

    inline void compute_equilibrium_UV(double& U, double& V);     // All in one Equilibirum solver function

    inline void compute_equilibrium_SV(double& S, double& V);     // All in one Equilibirum solver function

    inline void compute_equilibrium_TP(double& T, double& P);     // All in one Equilibirum solver function

    inline void compute_equilibrium_HP(double& H, double& P);     // All in one Equilibirum solver function

    inline void compute_equilibrium_SP(double& S, double& P);     // All in one Equilibirum solver function


    inline void findTRange();                       // Finds Temperature range to use proper NASA coeff.

    inline array<double, 7> temp_base();    // Calculates Temperature functions for NASA poly.

    inline void NASA_fits();                               // Calculates H0, S0, and MU0.

    inline void compute_mass_fractions();                  // Converts molar fractions to mass fractions.

    inline bool check_convergence(double* dlnj, double& dln);

    inline double compute_damping(double* dlnj, double& dln, double& dlnT);

    inline void form_xi();

    inline void formJF();

    inline void XI_TEST();

    inline bool CFD_convergence(double* dlnj);

};

#endif

