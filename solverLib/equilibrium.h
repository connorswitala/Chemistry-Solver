#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "../common_mixes/commonMixes.h"
#include "../includes/math.h"

using namespace std;



class equilibrium {
        
    public:

    mix gas;
    int T_flag; // Flags which set of NASA coefficients to use.
    int J_SIZE;
    int NCOEF;

    equilibrium(string gas_type);
    void display_gas_properties();   

    void compute_equilibrium(double rho, double e);

    private:

    Vector Ts;
    Vector X;

    inline void findTRange();
    inline array<double, 7> temp_base(double T);

    void NASA_fits();


    void compute_molar_fractions();
    void compute_mass_fractions();
    void redistribute_thermo(); 
    void compute_formation_enthalpies();


    double norm(double* v1, double* v2);
    // ====== Display Functions ======


};


#endif

