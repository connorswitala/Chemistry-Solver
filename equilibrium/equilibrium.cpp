#include "equilibrium.h"


equilibrium::equilibrium(string gas_type) {

    if (gas_type == "air11") {
        gas = common_air::make_air11();
        form_system = &equilibrium::form_system_neut;
    }
    if (gas_type == "air5") {
        gas = common_air::make_air5(); 
        form_system  = &equilibrium::form_system_neut;
    }
    if (gas_type == "air11_Ar") {
        gas = common_air::make_air11_Ar();
        form_system = &equilibrium::form_system_ions;
    }

    if (gas_type == "air13") {
        gas = common_air::make_air13();
        form_system = &equilibrium::form_system_ions;
    }

    if (gas_type == "perf") {
        gas = common_air::make_perf();
        return;
    }

    NCOEF = 9;
    NION = gas.N_ION;
    NSP = gas.N_SP;
    NEL = gas.N_EL;

}

void equilibrium::compute_equilibrium(double rho, double e) {
    
    if (gas.perf_flag) {
        gas.T = e / gas.cv;
        gas.rho = rho;
        gas.p = rho * gas.R * gas.T;
        gas.e = e;
        compute_mass_fractions();
        return; 
    }

    double e_new = 0, cv_new = 0;     
    gas.rho = rho;
    gas.e = e;

    if (e > 1e7) {
        gas.T = 15000.0;
    }
    else gas.T = gas.e / gas.cv;

    gas.p = rho * gas.R * gas.T; 


    int iteration = 0;

    J_SIZE = NEL + 1;
    X = vector<double>(J_SIZE);

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

void equilibrium::compute_equilibrium_TP(double T, double P) {
    
    if (gas.perf_flag) {
        gas.T = T;
        gas.p = P;
        gas.rho = gas.p / (gas.R * gas.T);
        gas.e = gas.cv * gas.T;
        return; 
    }

    gas.T = T;
    gas.p = P;
    gas.rho = gas.p / (gas.R * gas.T);
    gas.e = T * gas.cv;

    int iteration = 0;
    
    J_SIZE = NEL + 1;    
    
    g = vector<double>(NSP);
    lnN_old = vector<double>(NSP + 1, 0.0);
    lnN_new = vector<double>(NSP + 1, 0.0);
    dln = vector<double>(NSP + 1, 0.0);
    N = vector<double>(NSP + 1, 0.0);

    gas.N_tot = 0.1;
    double avg_N = gas.N_tot / NSP;

    for (int i = 0; i < NSP; ++i) {
        N[i] = avg_N;
        lnN_old[i] = log(avg_N);
    }

    lnN_old[NSP] = gas.N_tot;
    N[NSP] = log(gas.N_tot);

    NASA_fits();
    compute_molar_fractions();

    gas.e = 0.0;
    for (int i : gas.diatomic_list) {
        gas.e += gas.Y[i] * (2.5 * gas.species[i].R * gas.T 
                    + gas.species[i].R * gas.species[i].theta_v / (exp(gas.species[i].theta_v / gas.T) - 1)    
                    + gas.hf[i]); //diatomic species
    }

    for (int i : gas.mono_list) {
        gas.e += gas.Y[i] * (1.5 * gas.species[i].R * gas.T 
                    + gas.hf[i]); // atomic species
    }
    
    gas.cv = gas.e / gas.T;

    gas.MW = 0.0;
    for (int i = 0; i < NSP; ++i) {    
        gas.MW += gas.X[i] * gas.species[i].mw;
    }
    
    gas.R = gcon * 1000.0 / gas.MW;
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

        gas.H0[j] = gcon * gas.T * (- coeff[0] * Ts[0] 
                   + coeff[1] * Ts[1] * Ts[2] 
                   + coeff[2] 
                   + 0.5 * coeff[3] * Ts[3] 
                   + 0.333 * coeff[4] * Ts[4] 
                   + 0.25 * coeff[5] * Ts[5] 
                   + 0.2 * coeff[6] * Ts[6] 
                   + coeff[7] * Ts[1]);

        gas.S0[j] = gcon * (- 0.5 * coeff[0] * Ts[0] 
                   - coeff[1] * Ts[1]
                   + coeff[2] * Ts[2]
                   + coeff[3] * Ts[3] 
                   + 0.5 * coeff[4] * Ts[4] 
                   + 0.333 * coeff[5] * Ts[5] 
                   + 0.25 * coeff[6] * Ts[6] 
                   + coeff[8]);

        gas.mu0[j] = gas.H0[j] - gas.T * gas.S0[j];        
    }

    compute_formation_enthalpies();
}

void equilibrium::compute_formation_enthalpies() {
    
    for (int i : gas.reactant_idx) {
        gas.hf[i] = 0.0;
    }
    
    int k = 0;
    for (int i : gas.product_idx) {

        double sum = 0.0;

        for (int j = 0; j < gas.N_RE; ++j) {
            sum += gas.reactions[k * gas.N_RE + j] * gas.H0[gas.reactant_idx[j]];
        }
        gas.hf[i] = (gas.H0[i] - sum) * 1000.0 / gas.species[i].mw;
        k++;
    }

    gas.hf.back() = 0.0;     
}

void equilibrium::compute_molar_fractions() {

    vector<double> dx(J_SIZE, 0.0), F(J_SIZE, 0.0), J(J_SIZE * J_SIZE, 0.0); 

    int iteration = 0;    
    double RH_sum, e;
    bool converged = false;

    while (!converged) {


        (this->*form_system)(J.data(), F.data());        
        matrix_divide(J.data(), F.data(), dx.data(), J_SIZE, 1);

        // ---- compute delta ln(N_j) for each species
        for (int j = 0; j < NSP; ++j) {

            RH_sum = 0.0;
            for (int i = 0; i < NEL; ++i){
                RH_sum += gas.a[i * NSP + j] * dx[i];
            }

            dln[j] = dx[NEL] + RH_sum - g[j];
        }

        dln[NSP] = dx[NEL];

        // ---- find damping coefficient for update
        e = find_damping();
       

        // ---- apply damped update
        for (int j = 0; j < NSP; ++j) {
            lnN_new[j] = lnN_old[j] + e * dln[j];
            N[j]       = exp(lnN_new[j]);
        }
        lnN_new[NSP] = lnN_old[NSP] + e * dln[NSP];   // total
        N[NSP]       = exp(lnN_new[NSP]);

        // ---- check convergence and set old = new    
        converged = check_convergence(); 
        lnN_old = lnN_new;
        iteration++;
    }

    cout << "-- Iterations: " << iteration << endl;

    for (int i = 0; i < NSP; ++i) {
            gas.X[i] = N[i]/N[NSP];
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

    if (gas.perf_flag) return;

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


    // Buffers to hold (T, fractions...) so we can write Tecplot zones with correct I count
    std::vector<std::vector<double>> mass_rows;   // each row: [T, Y1, Y2, ...]
    std::vector<std::vector<double>> molar_rows;  // each row: [T, X1, X2, ...]

    int counter = 0;

    int ni = 194;

    for (int i = 0; i < ni; ++i) {
        double T = 500.0 + static_cast<double>(i) * 100.0;
        
        compute_equilibrium_TP(T, 101325.0);

        // Build a row for mass fractions
        std::vector<double> mass_row;
        mass_row.reserve(NSP + 1);
        mass_row.push_back(gas.T);
        for (int i = 0; i < NSP; ++i) mass_row.push_back(gas.Y[i]);
        mass_rows.push_back(std::move(mass_row));

        // Build a row for molar fractions
        std::vector<double> molar_row;
        molar_row.reserve(NSP + 1);
        molar_row.push_back(gas.T);
        for (int i = 0; i < NSP; ++i) molar_row.push_back(gas.X[i]);
        molar_rows.push_back(std::move(molar_row));

        counter++;
    }

    // Write Tecplot file with two zones
    std::ofstream tec("concentrations_tecplot.dat");
    tec << "TITLE = \"Equilibrium concentrations vs Temperature\"\n";

    // VARIABLES line: T plus each species (same for both zones)
    tec << "VARIABLES = \"T\"";
    for (int i = 0; i < NSP; ++i) tec << ", \"" << gas.species[i].name << "\"";
    tec << "\n";

    tec << std::scientific << std::setprecision(8);

    // ---- ZONE 1: Mass fractions ----
    tec << "ZONE T=\"Mass Fractions\", I=" << ni
        << ", F=POINT\n";
    for (const auto& row : mass_rows) {
        for (size_t j = 0; j < row.size(); ++j) {
            if (j) tec << ' ';
            tec << row[j];
        }
        tec << '\n';
    }

    // ---- ZONE 2: Molar fractions ----
    tec << "ZONE T=\"Molar Fractions\", I=" << static_cast<int>(molar_rows.size())
        << ", F=POINT\n";
    for (const auto& row : molar_rows) {
        for (size_t j = 0; j < row.size(); ++j) {
            if (j) tec << ' ';
            tec << row[j];
        }
        tec << '\n';
    }

    tec.close();
}

void equilibrium::add_system_ions(double* J, double* F) {

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
            F[i] = -(gas.mu0[i] + gcon * gas.T * log(Xi_safe * gas.p / 100000.0)) / (gcon * gas.T);
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

        vector<double> akj_N(NSP);
        double asum, gsum, stoich_sum, Nsum;
        double N_safe;

        double pref = 100000.0;

        for (int j = 0; j < NSP; ++j) {
            N_safe = max(N[j], 1e-10); // Protect log from zero
            g[j] = log(gas.p/pref * N_safe) + gas.mu0[j] / (gcon * gas.T);
        }

        // Loop through rows of Jacobian matrix
        for (int k = 0; k < NEL; ++k) {

            asum = 0.0;
            gsum = 0.0;

            for (int j = 0; j < NSP; ++j) {
                akj_N[j] = gas.a[k * NSP + j] * N[j];
                asum += akj_N[j];
            }

            // Loop through columns of row k for pi_i solution vector.
            for (int i = 0; i < NEL; ++i) {

                stoich_sum = 0.0;
                for (int j = 0; j < NSP; ++j) {
                    stoich_sum += akj_N[j] * gas.a[i * NSP + j];
                }

                J[k * J_SIZE + i] = stoich_sum;
            }

            J[k * J_SIZE + NEL] = asum;

            // Fill RHS vector.
            for (int j = 0; j < NSP; ++j) 
                gsum += akj_N[j] * g[j];
            
            F[k] = gas.b[k] - asum + gsum;
        }

        gsum = 0.0;
        Nsum = 0.0;

        for (int j = 0; j < NSP; ++j) {
            Nsum += N[j];
            gsum += N[j] * g[j];
        }

        for (int i = 0; i < NEL; ++i) {
            
            asum = 0.0;

            for (int j = 0; j < NSP; ++j) {
                asum += gas.a[i * NSP + j] * N[j];
            }

            J[NEL * J_SIZE + i] = asum;
        }

        J[NEL * J_SIZE + NEL] = Nsum - N[NSP];
        F[NEL] = N[NSP] - Nsum + gsum;

}

bool equilibrium::check_convergence() {

    double sum = 0.0, check, tol = 0.5e-5;

    for (int j = 0; j < NSP; ++j) 
        sum += N[j];
    

    for (int j = 0; j < NSP; ++j) {
        check = N[j] * fabs(dln[j]) / sum;
        if (check > tol) return false;
    }

    if (N[NSP] * fabs(dln[NSP]) / sum > tol)
        return false;

    return true;
}

double equilibrium::find_damping() {

    const double SIZE = -log(1e-8), eps = 1e-14;

     // ---- e1
    double maxAbs_dlnj = 0.0;
    for (int j = 0; j < NSP; ++j)
        maxAbs_dlnj = max(maxAbs_dlnj, abs(dln[j]));
    double e1 = 2.0 / max(5.0 * abs(dln[NSP]), maxAbs_dlnj);

    // ---- e2 (min over eligible species)
    double e2 = 1.0;
    bool any = false;
    for (int j = 0; j < NSP; ++j) {
        double frac_log = log(max(N[j],1e-300) / max(N[NSP],1e-300)); // ln(Nj/N)
        if (frac_log <= -SIZE && dln[j] >= 0.0) {
            double den = dln[j] - dln[NSP];
            if (abs(den) > eps) {
                double num = -frac_log - 9.2103404; // ln(1e4)
                double cand = abs(num / den);
                e2 = min(e2, cand);
                any = true;
            }
        }
    }
    if (!any) e2 = 1.0;

    // ---- final damping
    return min(1.0, min(e1, e2));
}



















