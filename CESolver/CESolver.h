#ifndef CESolver_H
#define CESolver_H

#include "../common_mixes/commonMixes.h"
#include "../energy_functions/functions.h"
#include "../includes/math.h"

using Task = void(*)(double*, double*, mix&);

class CESolver {
        
    public:

    mix gas;    // Main object.
    Context cfg; // Configuration Struct

    CESolver(GasType& gastype, ConstraintType& contrainttype);                       // Constructor.

    void compute_equilibrium_T(double rho, double e);     // All in one Equilibirum solver function

    void display_gas_properties();                      // Display results.

    void plot_concentrations_for_T_range();
    
    private:


    int T_flag; // Flags which set of NASA coefficients to use.

    int J_SIZE; // Size of solution vector in Newton Method.

    vector<MTask> matrix_tasks;
    gasTask mu_task;
    NjTask Nj_task;

    int NS, NE, NI;

    Vector X, g, lnN_old, lnN_new, dln, N;   // Old solution vector of molar concentrations.

    inline void findTRange();                       // Finds Temperature range to use proper NASA coeff.

    inline array<double, 7> temp_base(double T);    // Calculates Temperature functions for NASA poly.

    void NASA_fits();                               // Calculates H0, S0, and MU0.

    void compute_molar_fractions();                 // Newton Iteration to solve for concentrations.

    void compute_mass_fractions();                  // Converts molar fractions to mass fractions.

    void compute_formation_enthalpies();            // Computes formation enthalpies.

    bool check_convergence();

    double find_damping();

    void form_system_ions(double* J, double* F);
    void form_system_neut(double* J, double* F);

};

#endif

