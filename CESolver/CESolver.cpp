#include "CESolver.h"


CESolver::CESolver(GasType& gastype, ConstraintType& constrainttype) {

    switch (gastype) {
        case::GasType::AIR5:
            gas = common_air::make_air5();
            break;

        case::GasType::AIR11:
            gas = common_air::make_air11();
            break;

        case::GasType::AIR11_AR:
            gas = common_air::make_air11_Ar();
            break;

        case::GasType::AIR13:
            gas = common_air::make_air13();
            break;

        case::GasType::CREATE:
            // Will need to be made.
            break;
    }

    switch (constrainttype) {
        case::ConstraintType::TP:
            cfg.energy_type = EnergyType::Gibbs;
            cfg.constraint_type = ConstraintType::TP;
            cfg.NEEDS_T = false;
            break;

        case::ConstraintType::HP:
            cfg.energy_type = EnergyType::Gibbs;
            cfg.constraint_type = ConstraintType::HP;
            cfg.NEEDS_T = true;
            break;

        case::ConstraintType::SP:
            cfg.energy_type = EnergyType::Gibbs;
            cfg.constraint_type = ConstraintType::SP;
            cfg.NEEDS_T = true;
            break;

        case::ConstraintType::TV:
            cfg.energy_type = EnergyType::Helmholtz;
            cfg.constraint_type = ConstraintType::TV;
            cfg.NEEDS_T = false;
            break;

        case::ConstraintType::UV:
            cfg.energy_type = EnergyType::Helmholtz;
            cfg.constraint_type = ConstraintType::UV;
            cfg.NEEDS_T = true;
            break;

        case::ConstraintType::SV:
            cfg.energy_type = EnergyType::Helmholtz;
            cfg.constraint_type = ConstraintType::SV;
            cfg.NEEDS_T = true;
            break;
    }

    if (gas.HAS_IONS) cfg.HAS_IONS = true;
    else cfg.HAS_IONS = false;

    matrix_tasks = energy::build_matrix_tasks(cfg); 
    mu_task = energy::build_mu_task(cfg);
    Nj_task = energy::build_Nj_task(cfg);

    NI = gas.NI;
    NS = gas.NS;
    NE = gas.NE;
    J_SIZE = NE + gas.HAS_IONS + cfg.NEEDS_T;
    gas.J_SIZE = J_SIZE;

    cout << "-- Configuration complete" << endl;

}

void CESolver::compute_equilibrium_T(double e, double rho) {
    
    gas.V = 1.0 / rho;
    gas.T = 3000.0;

    int iteration = 0;

    vector<double> J(J_SIZE * J_SIZE), F(J_SIZE), DELTA(J_SIZE), corr_vars(J_SIZE); 

    double N_sum  = 0.0;
    for (int j = 0; j < NS; ++j) N_sum += gas.X0[j];
    for (int j = 0; j < NS; ++j) {
        gas.N[j] = N_sum / NS;
        corr_vars[j] = 
    }

    while (iteration < 25) {

        NASA_fits();
        mu_task;

        for (MTask FN : matrix_tasks) 
            FN(J.data(), F.data(), gas);

        matrix_divide(J.data(), F.data(), DELTA.data(), J_SIZE, 1);

        iteration++;
    }



    cout << "-- Equilibrium computed" << endl;
}


inline void CESolver::findTRange() {
    if      (gas.T >= 200.0   && gas.T <= 1000.0)   T_flag = 0; 
    else if (gas.T >  1000.0  && gas.T <= 6000.0)   T_flag = 1; 
    else if (gas.T >  6000.0  && gas.T <= 20000.0)  T_flag = 2;
    else T_flag = -1;
}

inline array<double, 7> CESolver::temp_base(double T) {

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

void CESolver::NASA_fits() {

    // H(T) / RT = -a0 * T^(-2) + a1 * ln(T)/T + a2 + a3 * T/2 + a4 * T^2 / 3 + a5 * T^3/4 + a6 * T^4 / 5 + a7 / T
    // S(T) / R = -a0 * T^(-2)/2 - a1 / T + a2 * ln(T) + a3 * T + a4 * T^2 / 2 + a5 * T^3/3 + a6 * T^4 / 4 + a8

    findTRange();                           // Find the temperature range to use for NASA Polynomials
    const auto Ts = temp_base(gas.T);       // Calculate Temperature variables.

    for (int j = 0; j < NS; ++j) {

        const auto& poly = gas.species[j].poly;
        const double* coeff = &poly.at(T_flag * NCOEF);

        gas.H0_RT[j] = - coeff[0] * Ts[0] 
                   + coeff[1] * Ts[1] * Ts[2] 
                   + coeff[2] 
                   + 0.5 * coeff[3] * Ts[3] 
                   + 0.333 * coeff[4] * Ts[4] 
                   + 0.25 * coeff[5] * Ts[5] 
                   + 0.2 * coeff[6] * Ts[6] 
                   + coeff[7] * Ts[1];

        gas.S0_R[j] = - 0.5 * coeff[0] * Ts[0] 
                   - coeff[1] * Ts[1]
                   + coeff[2] * Ts[2]
                   + coeff[3] * Ts[3] 
                   + 0.5 * coeff[4] * Ts[4] 
                   + 0.333 * coeff[5] * Ts[5] 
                   + 0.25 * coeff[6] * Ts[6] 
                   + coeff[8];

        gas.CP0_R[j] = coeff[0] * Ts[0]
                     + coeff[1] * Ts[1]
                     + coeff[2]
                     + coeff[3] * Ts[3]
                     + coeff[4] * Ts[4]
                     + coeff[5] * Ts[5]
                     + coeff[6] * Ts[6];

        gas.mu0_RT[j] = gas.H0_RT[j] - gas.S0_R[j];        
    }

    compute_formation_enthalpies();
}

void CESolver::compute_formation_enthalpies() {
    
    for (int i : gas.reactant_idx) {
        gas.hf[i] = 0.0;
    }
    
    int k = 0;
    for (int i : gas.product_idx) {

        double sum = 0.0;

        for (int j = 0; j < gas.N_RE; ++j) {
            sum += gas.reactions[k * gas.N_RE + j] * gas.H0_RT[gas.reactant_idx[j]];
        }
        gas.hf[i] = (gas.H0_RT[i] - sum) * 1000.0 / gas.species[i].mw;
        k++;
    }

    gas.hf.back() = 0.0;     
}

void CESolver::compute_molar_fractions() {

    vector<double> dx(J_SIZE, 0.0), F(J_SIZE, 0.0), J(J_SIZE * J_SIZE, 0.0); 

    int iteration = 0;    
    double RH_sum, e;
    bool converged = false;

    while (!converged) {
     
        matrix_divide(J.data(), F.data(), dx.data(), J_SIZE, 1);

        // ---- compute delta ln(N_j) for each species
        for (int j = 0; j < NS; ++j) {

            RH_sum = 0.0;
            for (int i = 0; i < NE; ++i){
                RH_sum += gas.a[i * NS + j] * dx[i];
            }

            dln[j] = dx[NE] + RH_sum - g[j];
        }

        dln[NS] = dx[NE];

        // ---- find damping coefficient for update
        e = find_damping();
       

        // ---- apply damped update
        for (int j = 0; j < NS; ++j) {
            lnN_new[j] = lnN_old[j] + e * dln[j];
            N[j]       = exp(lnN_new[j]);
        }
        lnN_new[NS] = lnN_old[NS] + e * dln[NS];   // total
        N[NS]       = exp(lnN_new[NS]);

        // ---- check convergence and set old = new    
        converged = check_convergence(); 
        lnN_old = lnN_new;
        iteration++;
    }

    cout << "-- Iterations: " << iteration << endl;

    for (int i = 0; i < NS; ++i) {
            gas.X[i] = N[i]/N[NS];
    }

    compute_mass_fractions();
}

void CESolver::compute_mass_fractions() {

        gas.MW = 0.0;
        for (int i = 0; i < NS; ++i) {
                gas.MW += gas.X[i] * gas.species[i].mw;
        }

        for (int i = 0; i < NS; ++i) {
                gas.Y[i] = (gas.X[i] * gas.species[i].mw) / gas.MW;
        }
}

void CESolver::display_gas_properties() {


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
    for (int i = 0; i < NS; ++i) {
        cout << setw(10) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.X[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }

    cout << endl << endl << setw(67) <<  "====: Species Mass Fractions :==== " << endl;
    for (int i = 0; i < NS; ++i) {
        cout << setw(10) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.Y[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }
}

bool CESolver::check_convergence() {

    double sum = 0.0, check, tol = 0.5e-5;

    for (int j = 0; j < NS; ++j) 
        sum += N[j];
    

    for (int j = 0; j < NS; ++j) {
        check = N[j] * fabs(dln[j]) / sum;
        if (check > tol) return false;
    }

    if (N[NS] * fabs(dln[NS]) / sum > tol)
        return false;

    return true;
}

double CESolver::find_damping() {

    const double SIZE = -log(1e-8), eps = 1e-14;

     // ---- e1
    double maxAbs_dlnj = 0.0;
    for (int j = 0; j < NS; ++j)
        maxAbs_dlnj = max(maxAbs_dlnj, abs(dln[j]));
    double e1 = 2.0 / max(5.0 * abs(dln[NS]), maxAbs_dlnj);

    // ---- e2 (min over eligible species)
    double e2 = 1.0;
    bool any = false;
    for (int j = 0; j < NS; ++j) {
        double frac_log = log(max(N[j],1e-300) / max(N[NS],1e-300)); // ln(Nj/N)
        if (frac_log <= -SIZE && dln[j] >= 0.0) {
            double den = dln[j] - dln[NS];
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



















