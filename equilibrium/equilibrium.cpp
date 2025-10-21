#include "equilibrium.h"


equilibrium::equilibrium(string gas_type) {

    if (gas_type == "air11") {
        gas = common_air::make_air11();
        compute_hf = &equilibrium::compute_formation_enthalpies_ions;
        form_system = &equilibrium::form_system_ions;
    }
    if (gas_type == "air5") {
        gas = common_air::make_air5(); 
        compute_hf = &equilibrium::compute_formation_enthalpies_neut;
        form_system  = &equilibrium::form_system_neut;
    }
    if (gas_type == "air11_Ar") {
        gas = common_air::make_air11_Ar();
        compute_hf = &equilibrium::compute_formation_enthalpies_ions;
        form_system = &equilibrium::form_system_ions;
    }

    if (gas_type == "air13") {
        gas = common_air::make_air13();
        compute_hf = &equilibrium::compute_formation_enthalpies_ions;
        form_system = &equilibrium::form_system_ions;
    }

    NCOEF = 9;
    NION = gas.N_ION;
    NSP = gas.N_SP;
    NEL = gas.N_EL;

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

    J_SIZE = NSP + NEL + 1;    
    X = vector<double>(J_SIZE);

    findTRange();

    X[J_SIZE - 1] = gas.initial_moles;
    double avg_x = gas.initial_moles / NSP;
    for (int i = 0; i < NSP; ++i) X[i] = avg_x;
    while (fabs(e_new - e) >= 100) {

        NASA_fits();
        compute_molar_fractions();

        e_new = 0.0;
        for (int i : gas.diatomic_list) {
            e_new += gas.Y[i] * (2.5 * gas.species[i].R * gas.T 
                     + gas.species[i].R * gas.species[i].theta_v / (exp(gas.species[i].theta_v / gas.T) - 1)    
                     + gas.hf[i]); //diatomic species
        }

        for (int i : gas.mono_list) {
            e_new += gas.Y[i] * (1.5 * gas.species[i].R * gas.T 
                     + gas.hf[i]); // atomic species
        }
        
        cv_new = e_new / gas.T;
        gas.T = gas.T - 0.5 * (e_new - e) / cv_new;

        gas.MW = 0.0;
        for (int i = 0; i < NSP; ++i) {    
            gas.MW += gas.X[i] * gas.species[i].mw;
        }
        
        gas.R = gcon * 1000.0 / gas.MW;
        gas.p = gas.rho * gas.R * gas.T;
        gas.cv = cv_new;
        iteration++;   
    }

    gas.cp = gas.R + gas.cv;
    gas.gamma = gas.cp / gas.cv;
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

    for (int j = 0; j < NSP; ++j) {

        const auto& poly = gas.species[j].poly;
        const double* coeff = &poly.at(T_flag * NCOEF);

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

        gas.mu0[j] = gas.H0[j] - gas.T * gas.S0[j];        
    }

    (this->*compute_hf)();
}

void equilibrium::compute_formation_enthalpies_ions() {
    
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

void equilibrium::compute_formation_enthalpies_neut() {
    
    gas.hf[0] = 0.0;    // N2
    gas.hf[1] = 0.0;    // O2
    gas.hf[2] = (gas.H0[2] - 0.5 * gas.H0[0] - 0.5 * gas.H0[1]) * 1000 / gas.species[2].mw;  // NO
    gas.hf[3] = (gas.H0[3] - 0.5 * gas.H0[0]) * 1000 / gas.species[3].mw;    // N
    gas.hf[4] = (gas.H0[4] - 0.5 * gas.H0[1]) * 1000 / gas.species[4].mw;    // O

}

void equilibrium::compute_molar_fractions() {

    vector<double> X_new(J_SIZE, 0.0), dx(J_SIZE, 0.0), F(J_SIZE, 0.0), J(J_SIZE * J_SIZE, 0.0); 

    int iteration = 0;    
    double residual = norm(X_new.data(), X.data(), NSP);

    while (residual > 1e-8) {

        fill(J.begin(), J.end(), 0.0);
        (this->*form_system)(J.data(), F.data());        
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

        residual = norm(X_new.data(), X.data(), NSP);
        X = X_new;
        iteration++;
    }

    double sum = 0.0;
    for (int i = 0; i < NSP; ++i) sum += X[i];

    for (int i = 0; i < NSP; ++i) {
            gas.X[i] = X_new[i]/sum;
    }

    compute_mass_fractions();
}

void equilibrium::compute_mass_fractions() {

        gas.MW = 0.0;
        for (int i = 0; i < NSP; ++i) {
                gas.MW += gas.X[i] * gas.species[i].mw;
        }

        for (int i = 0; i < NSP; ++i) {
                gas.Y[i] = (gas.X[i] * gas.species[i].mw) / gas.MW;
        }
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


    cout << endl << setw(67) << "====: Species Molar Fractions :==== " << endl;
    for (int i = 0; i < NSP; ++i) {
        cout << setw(10) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.X[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }

    cout << endl << endl << setw(67) <<  "====: Species Mass Fractions :==== " << endl;
    for (int i = 0; i < NSP; ++i) {
        cout << setw(10) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.Y[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }


}

void equilibrium::plot_concentrations_for_T_range() {
    
    double rho = 0.1, e;
    double T_eq = 1300;

    ofstream mass("mass_fractions.csv");
    ofstream molar("molar_fractions.csv");

    mass << "T, N2, O2, NO, N, O, Ar, Ar+, N+, O+, NO+, e-" << endl;
    molar << "T, N2, O2, NO, N, O, Ar, Ar+, N+, O+, NO+, e-" << endl;

    int counter = 0;
    while (T_eq < 20000) {
            e = 717 * T_eq + 20000 * counter;
            compute_equilibrium(rho, e);

            mass << gas.T; // Temperature
            molar << gas.T; // Temperature
            for (int i = 0; i < NSP; ++i) mass << ", " << gas.Y[i];  // Mass fractions 
            for (int i = 0; i < NSP; ++i) molar << ", " << gas.X[i];  // Molar fractions 
            mass << "\n";
            molar << "\n";

            counter++; 
            T_eq = gas.T;
            if (counter % 100 == 0) cout << T_eq << endl;
    }

    mass.close();
    molar.close();
}

void equilibrium::form_system_ions(double* J, double* F) {

        for (int i = NSP; i < J_SIZE; ++i) {
            X[i] = 0.0;
        } 

        for (int i = 0; i < NSP; ++i) {

            J[i * J_SIZE + i] = 1.0;  // Identity part

            for (int k = 0; k < NEL; ++k) {
                J[i * J_SIZE + NSP + k] = -gas.a[k * NSP + i];  // -a[k][i] 
            }

            J[i * J_SIZE + J_SIZE - 1] = -gas.species[i].q; 

        }       

        // Row 10–12: nitrogen, oxygen, argon atomic balance
        for (int k = 0; k < NEL; ++k) {
            for (int j = 0; j < NSP; ++j) {
                J[(NSP + k) * J_SIZE + j] = gas.a[k * NSP + j] * X[j];  // a[k][j] * X[j]
            }
        // last 4 entries (cols 10–13) remain 0
        }

        // Row 13: charge neutrality
        for (int j = 0; j < NSP; ++j) {
            J[(J_SIZE - 1) * J_SIZE + j] = gas.species[j].q * X[j];
        }
        // last 4 entries of row 13 are 0

        for (int i = 0; i < NSP; ++i) {
            double Xi_safe = max(X[i], 1e-10); // Protect log from zero  
            F[i] = -(gas.mu0[i] + gcon * gas.T * log(Xi_safe * gas.p / 101325.0)) / (gcon * gas.T);
        }

        for (int i = 0; i < NEL; ++i) {
            double sum = 0.0;
            for (int j = 0; j < NSP; ++j) sum += gas.a[i * NSP + j] * X[j];
            F[NSP + i] = gas.b[i] - sum;
        }

        double sum = 0.0;
        for (int i = NSP - NION; i < NSP; ++i) sum += X[i] * gas.species[i].q;
        F[J_SIZE - 1] = -sum;
}

void equilibrium::form_system_neut(double* J, double* F) {

        for (int i = NSP; i < J_SIZE - 1; ++i) {
            X[i] = 0.0;
        } 

        for (int i = 0; i < NSP; ++i) {

            J[i * J_SIZE + i] = 1.0;  // Identity part

            for (int k = 0; k < NEL; ++k) {
                J[i * J_SIZE + NSP + k] = -gas.a[k * NSP + i];  // -a[k][i] 
            }

            J[i * J_SIZE + J_SIZE - 1] = -1.0;
        }       

        for (int k = 0; k < NEL; ++k) {
            for (int j = 0; j < NSP; ++j) {
                J[(NSP + k) * J_SIZE + j] = gas.a[k * NSP + j] * X[j];  // a[k][j] * X[j]   
            }
        // last 4 entries (cols 10–13) remain 0
        }

        for (int j = 0; j < NSP; ++j) {
            J[(J_SIZE - 1) * J_SIZE + j] = X[j];
        }

        J[J_SIZE * J_SIZE - 1] = -X[J_SIZE - 1];

        for (int i = 0; i < NSP; ++i) {
            double Xi_safe = max(X[i], 1e-10); // Protect log from zero
            F[i] = -(gas.mu0[i] + gcon * gas.T * log(Xi_safe * gas.p / 101325.0)) / (gcon * gas.T);
        }

        for (int i = 0; i < NEL; ++i) {
            double sum = 0.0;
            for (int j = 0; j < NSP; ++j) sum += gas.a[i * NSP + j] * X[j];
            F[NSP + i] = gas.b[i] - sum;
        }

        double sum = 0.0;
        for (int i = 0; i < NSP; ++i) sum += X[i];
        F[J_SIZE - 1] = gas.initial_moles - sum;
}























