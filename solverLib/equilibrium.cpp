#include "equilibrium.hpp"


equilibrium::equilibrium() {

    int 
    cout << "Declare mixture. (Currently just type '5' or '11' for 5 or 11 species air.): ";
    cin >> 

}


int equilibrium::findTRange() {
    if      (gas.T >= 200.0   && gas.T <= 1000.0)  { T_flag = 0; return 0; }
    else if (gas.T >  1000.0  && gas.T <= 6000.0)  { T_flag = 1; return 1; }
    else if (gas.T >  6000.0  && gas.T <= 20000.0) { T_flag = 2; return 2; }
    T_flag = -1;
    return -1; // caller decides what to do
}

array<double, 7> equilibrium::temp_base(double T) {

    array<double, 7> Ts;
    double T_inverse = 1/T;

    Ts[0] = T_inverse * T_inverse;  // 1 / T^2
    Ts[1] = T_inverse;  // 1 / T
    Ts[2] = log(T);     // ln(T)
    Ts[3] = T;          // T
    Ts[4] = Ts[3] * T;  // T^2
    Ts[5] = Ts[4] * T;  // T^3
    Ts[6] = Ts[5] * T;  // T^4
    
    return Ts;
}

void equilibrium::NASA_fits() {

    // H(T) / RT = -a0 * T^(-2) + a1 * ln(T)/T + a2 + a3 * T/2 + a4 * T^2 / 3 + a5 * T^3/4 + a6 * T^4 / 5 + a7 / T
    // S(T) / R = -a0 * T^(-2)/2 - a1 / T + a2 * ln(T) + a3 * T + a4 * T^2 / 2 + a5 * T^3/3 + a6 * T^4 / 4 + a8

    const auto Ts = temp_base(gas.T);    

    for (int j = 0; j < n_sp; ++j) {

        const auto& poly = gas.species[j].poly;
        const double* coeff = poly.data() + T_flag * NCOEF;

        gas.H0[j] = gcon * gas.T * (
                   -coeff[0] * Ts[0] 
                   + coeff[1] * Ts[1] * Ts[2] 
                   + coeff[2] 
                   + 0.5 * coeff[3] * Ts[3] 
                   + 0.333 * coeff[4] * Ts[4] 
                   + 0.25 * coeff[5] * Ts[5] 
                   + 0.2 * coeff[6] * Ts[6] 
                   + coeff[7] * Ts[1]);

        gas.S0[j] = gcon * (
                   -0.5 * coeff[0] * Ts[0] 
                   - coeff[1] * Ts[1]
                   + coeff[2] * Ts[2]
                   + coeff[3] * Ts[3] 
                   + 0.5 * coeff[4] * Ts[4] 
                   + 0.333 * coeff[5] * Ts[5] 
                   + 0.25 * coeff[6] * Ts[6] 
                   + coeff[8]);

        gas.mu0[j] = gas.H0[j] - gas.T * gas.S0[j];
    }
}


void equilibrium::create_air_11() {

    


}