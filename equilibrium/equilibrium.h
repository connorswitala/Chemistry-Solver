#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "../common_mixes/commonMixes.h"
#include "../includes/math.h"

using namespace std;

class equilibrium {
        
    public:

    mix gas;    // Main object.

    equilibrium(string gas_type);                       // Constructor.

    void compute_equilibrium(double rho, double e);     // All in one Equilibirum solver function

    void display_gas_properties();                      // Display results.

    void plot_concentrations_for_T_range();
    
    private:

    int T_flag; // Flags which set of NASA coefficients to use.

    int J_SIZE; // Size of solution vector in Newton Method.

    int NCOEF;  // Number of NASA polynomial coefficients (Always 9).

    int NSP, NEL, NION;

    Vector X;   // Old solution vector of molar concentrations..


    inline void findTRange();                       // Finds Temperature range to use proper NASA coeff.

    inline array<double, 7> temp_base(double T);    // Calculates Temperature functions for NASA poly.

    void NASA_fits();                               // Calculates H0, S0, and MU0.

    void compute_molar_fractions();                 // Newton Iteration to solve for concentrations.

    void compute_mass_fractions();                  // Converts molar fractions to mass fractions.

    void compute_formation_enthalpies_ions();            // Computes formation enthalpies.
    void compute_formation_enthalpies_neut();            // Computes formation enthalpies.

    void form_system_ions(double* J, double* F);
    void form_system_neut(double* J, double* F);

    using hfptr = void (equilibrium::*)();
    using systemptr = void (equilibrium::*)(double*, double*);

    systemptr form_system = nullptr;
    hfptr compute_hf = nullptr;

};

#endif

