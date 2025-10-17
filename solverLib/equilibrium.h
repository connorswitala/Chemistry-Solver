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

    equilibrium(string gas_type, double rho, double e);
    void diplay_gas_properties();   

    void compute_equilibrium();

    private:

    Vector Ts;

    inline void findTRange();
    inline array<double, 7> temp_base(double T);

    void NASA_fits();


    void compute_molar_fractions();
    void compute_mass_fractions();
    void redistribute_thermo(); 


    double norm(double* v1, double* v2);
    // ====== Display Functions ======


};


#endif

