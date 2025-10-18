#include "equilibrium.h"


equilibrium::equilibrium(string gas_type) {

    if (gas_type == "air11") gas = common_air::make_air11();
    if (gas_type == "air5") gas = common_air::make_air5(); 
    NCOEF = 9;
}

inline void equilibrium::findTRange() {
    if      (gas.T >= 200.0   && gas.T <= 1000.0)   T_flag = 0; 
    else if (gas.T >  1000.0  && gas.T <= 6000.0)   T_flag = 1; 
    else if (gas.T >  6000.0  && gas.T <= 20000.0)  T_flag = 2;
    else T_flag = -1;
}

inline array<double, 7> equilibrium::temp_base(double T) {

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

    findTRange();                           // Find the temperature range to use for NASA Polynomials
    const auto Ts = temp_base(gas.T);       // Calculate Temperature variables.

    for (int j = 0; j < gas.N_SP; ++j) {

        const auto& poly = gas.species[j].poly;
        const double* coeff = poly.data() + T_flag * NCOEF;

        gas.H0[j] = gcon * gas.T * (
                   - coeff[0] * Ts[0] 
                   + coeff[1] * Ts[1] * Ts[2] 
                   + coeff[2] 
                   + 0.5 * coeff[3] * Ts[3] 
                   + 0.333 * coeff[4] * Ts[4] 
                   + 0.25 * coeff[5] * Ts[5] 
                   + 0.2 * coeff[6] * Ts[6] 
                   + coeff[7] * Ts[1]);

        gas.S0[j] = gcon * (
                   - 0.5 * coeff[0] * Ts[0] 
                   - coeff[1] * Ts[1]
                   + coeff[2] * Ts[2]
                   + coeff[3] * Ts[3] 
                   + 0.5 * coeff[4] * Ts[4] 
                   + 0.333 * coeff[5] * Ts[5] 
                   + 0.25 * coeff[6] * Ts[6] 
                   + coeff[8]);

        gas.Cp0[j] = gcon * (
                       coeff[0] * Ts[0] 
                     + coeff[1] * Ts[1] 
                     + coeff[2] 
                     + coeff[3] * Ts[3] 
                     + coeff[4] * Ts[4] 
                     + coeff[5] * Ts[5] 
                     + coeff[6] * Ts[6] );

        gas.mu0[j] = gas.H0[j] - gas.T * gas.S0[j];        
    }
}

void equilibrium::compute_equilibrium(double rho, double e) {
    
    double e_new = 0, cv_new = 0;     
    gas.rho = rho;
    gas.e = e;

    if (e > 1e7) {
        gas.T = 15000.0;
    }
    else gas.T = gas.e / gas.cv;

    gas.p = rho * gas.R * gas.T; 


    int iteration = 0;

    J_SIZE = gas.N_SP + gas.N_EL + 1;

    X = vector<double>(J_SIZE);
    findTRange();
    for (int i = 0; i < gas.N_SP; ++i) {
        X[i] = gas.guesses[T_flag * gas.N_SP + i];
    }   
    
    // double avg_x = gas.initial_moles / gas.N_SP;
    // for (int i = 0; i < gas.N_SP; ++i) X[i] = avg_x;

    while (fabs(e_new - e) >= 100) {

        findTRange();
        NASA_fits();
        compute_molar_fractions();
        compute_formation_enthalpies();

        e_new = 0.0;
        for (int i : gas.diatomic_list) {
            e_new += gas.Y[i] * (2.5 * gas.species[i].R * gas.T + gas.species[i].R * gas.species[i].theta_v / (exp(gas.species[i].theta_v / gas.T) - 1) + gas.hf[i]); //diatomic species
        }

        for (int i : gas.mono_list) {
            e_new += gas.Y[i] * (1.5 * gas.species[i].R * gas.T + gas.hf[i]); // atomic species
        }
        
        cv_new = e_new / gas.T;
        gas.T = gas.T - 0.6 * (e_new - e) / cv_new;

        gas.MW = 0.0;
        for (int i = 0; i < gas.N_SP; ++i) {    
            gas.MW += gas.X[i] * gas.species[i].mw;
        }
        
        gas.R = gcon * 1000.0 / gas.MW;
        gas.p = gas.rho * gas.R * gas.T;
        gas.cv = cv_new;
        iteration++;
        
    }

    gas.cp = gas.R + gas.cv;
    gas.gamma = gas.cp / gas.cv;
    cout << "-- Outer iterations: " << iteration << endl;
}



void equilibrium::compute_molar_fractions() {
    vector<double> X_new(J_SIZE, 0.0), dx(J_SIZE, 0.0), F(J_SIZE, 0.0), J(J_SIZE * J_SIZE, 0.0); 

    int nsp = gas.N_SP;
    int nel = gas.N_EL;
    int nion = gas.N_ION;

    int iteration = 0;
    
    double residual = norm(X_new.data(), X.data());

    while (residual > 1e-6) {

        for (int i = nsp; i < J_SIZE; ++i) {
            X[i] = 0.0;
        }

        for (int i = 0; i < nsp; ++i) {

            J[i * J_SIZE + i] = 1.0;  // Identity part

            for (int k = 0; k < nel; ++k) {
                    J[i * J_SIZE + nsp + k] = -gas.a[k * nsp + i];  // -a[k][i] 
            }

            J[i * J_SIZE + J_SIZE - 1] = -gas.species[i].q; 

        }       

        // Row 10–12: nitrogen, oxygen, argon atomic balance
        for (int k = 0; k < nel; ++k) {
                for (int j = 0; j < nsp; ++j) {
                        J[(nsp + k) * J_SIZE + j] = gas.a[k * nsp + j] * X[j];  // a[k][j] * X[j]
                }
        // last 4 entries (cols 10–13) remain 0
        }

        // Row 13: charge neutrality
        for (int j = 0; j < nsp; ++j) {
                J[(J_SIZE - 1) * J_SIZE + j] = gas.species[j].q * X[j];
        }
        // last 4 entries of row 13 are 0

        for (int i = 0; i < nsp; ++i) {
                double Xi_safe = max(X[i], 1e-10); // Protect log from zero  
                F[i] = -(gas.mu0[i] + gcon * gas.T * log(Xi_safe * gas.p / 101325.0)) / (gcon * gas.T);
        }

        for (int i = 0; i < nel; ++i) {
            double sum = 0.0;
            for (int j = 0; j < nsp; ++j) sum += gas.a[i * nsp + j] * X[j];
            F[nsp + i] = gas.b[i] - sum;
        }

        double sum = 0.0;
        for (int i = nsp - nion; i < nsp; ++i) sum += X[i] * gas.species[i].q;
        F[J_SIZE - 1] = -sum;
        
        matrix_divide(J.data(), F.data(), dx.data(), J_SIZE, 1);

        for (int i = 0; i < J_SIZE; ++i) {
                X_new[i] = X[i] * exp(dx[i]);

                if (isnan(X_new[i]) || isinf(X_new[i])) {
                        cout << "Nonphysical molar fraction computed (NaN or Inf)" << endl;
                }

                if (X_new[i] < 0) {
                        cout << "Negative molar fraction computed" << endl;
                }
        }

        residual = norm(X_new.data(), X.data());
        X = X_new;
        iteration++;
    }

    double sum = 0.0;
    for (int i = 0; i < nsp; ++i) sum += X[i];

    for (int i = 0; i < nsp; ++i) {
            gas.X[i] = X_new[i]/sum;
    }

    compute_mass_fractions();

}

void equilibrium::compute_mass_fractions() {

        gas.MW = 0.0;
        for (int i = 0; i < gas.N_SP; ++i) {
                gas.MW += gas.X[i] * gas.species[i].mw;
        }

        for (int i = 0; i < gas.N_SP; ++i) {
                gas.Y[i] = (gas.X[i] * gas.species[i].mw) / gas.MW;
        }
}

void equilibrium::compute_formation_enthalpies() {
    
    gas.hf[0] = 0.0;    // N2
    gas.hf[1] = 0.0;    // O2
    gas.hf[2] = (gas.H0[2] - 0.5 * gas.H0[0] - 0.5 * gas.H0[1]) * 1000 / gas.species[2].mw;  // NO
    gas.hf[3] = (gas.H0[3] - 0.5 * gas.H0[0]) * 1000 / gas.species[3].mw;    // N
    gas.hf[4] = (gas.H0[4] - 0.5 * gas.H0[1]) * 1000 / gas.species[4].mw;    // O
    gas.hf[5] = 0.0;    // Ar
    gas.hf[6] = (gas.H0[6] - gas.H0[5]) * 1000 / gas.species[6].mw;  // Ar+
    gas.hf[7] = (gas.H0[7] - 0.5 * gas.H0[0]) * 1000 / gas.species[7].mw;    // N+
    gas.hf[8] = (gas.H0[8] - 0.5 * gas.H0[1]) * 1000 / gas.species[8].mw;    // O+
    gas.hf[9] = (gas.H0[9] - 0.5 * gas.H0[0] - 0.5 * gas.H0[1]) * 1000 / gas.species[9].mw;   // NO+
    gas.hf[10] = 0.0;

}


double equilibrium::norm(double* v1, double* v2) {
        double result = 0.0;
        for (int i = 0; i < 14; ++i) {
                result += fabs(v1[i] - v2[i]) * fabs(v1[i] - v2[i]);
        }
        return sqrt(result); 
}

void equilibrium::display_gas_properties() {
    cout << endl << setw(30) << "Mixture Temperature: " << gas.T << " [K]" << endl;
    cout << setw(30) << "Mixture Pressure: " << gas.p / 1000.0 << " [kPa]" << endl;
    cout << setw(30) << "Mixture Density: " << gas.rho << " [kg/m^3]" << endl;
    cout << setw(30) << "Mixture Internal Energy: " << gas.e / 1000.0 << " [kJ/kg]" << endl;
    cout << setw(30) << "Mixture Gas Constant: " << gas.R << " [J/kg-K]" << endl;
    cout << setw(30) << "Mixture cp: " << gas.cp << " [J/kg-K]" << endl;
    cout << setw(30) << "Mixture cv: " << gas.cv << " [J/kg-K]" << endl;
    cout << setw(30) << "Mixture gamma: " << gas.gamma << endl;
    cout << setw(30) << "Mixture MW: " << gas.MW << " [g/mol]" << endl;


    cout << endl << "-- Species Molar Fractions -- " << endl;
    for (int i = 0; i < gas.N_SP; ++i) {
        cout << setw(3) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.X[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }

    cout << endl << endl << "-- Species Mass Fractions -- " << endl;
    for (int i = 0; i < gas.N_SP; ++i) {
        cout << setw(3) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.Y[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }


}