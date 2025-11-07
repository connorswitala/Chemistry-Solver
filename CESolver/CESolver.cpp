#include "CESolver.h"


CESolver::CESolver(mix& gas_in, ConstraintType& constrainttype) : gas(gas_in) {

    string energy;
    string contraint;

    switch (constrainttype) {
        case::ConstraintType::TP:
            gas.NEEDS_T = false;
            equilibrium = &CESolver::compute_equilibrium_TP;
            energy = "Gibbs";
            contraint = "temperature and pressure";
            NS = gas.NS;
            NE = gas.NE;
            J_SIZE = NE + 1 + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;

        case::ConstraintType::HP:
            gas.NEEDS_T = true;
            equilibrium = &CESolver::compute_equilibrium_HP;
            energy = "Gibbs";
            contraint = "enthalpy and pressure";
            NS = gas.NS;    
            NE = gas.NE;
            J_SIZE = NE + 1 + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;

        case::ConstraintType::SP:
            gas.NEEDS_T = true;
            equilibrium = &CESolver::compute_equilibrium_SP;
            energy = "Gibbs";
            contraint = "entropy and pressure";
            NS = gas.NS;
            NE = gas.NE;
            J_SIZE = NE + 1 + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;

        case::ConstraintType::TV:
            gas.NEEDS_T = false;
            equilibrium = &CESolver::compute_equilibrium_TV;
            energy = "Helmholtz";
            contraint = "temperature and volume";
            NS = gas.NS;
            NE = gas.NE;
            J_SIZE = NE + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;

        case::ConstraintType::UV:
            gas.NEEDS_T = true;
            equilibrium = &CESolver::compute_equilibrium_UV;
            energy = "Helmholtz";
            contraint = "internal energy and volume";
            NS = gas.NS;
            NE = gas.NE;
            J_SIZE = NE + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;

        case::ConstraintType::SV:
            gas.NEEDS_T = true;
            equilibrium = &CESolver::compute_equilibrium_SV;
            energy = "Helmholtz";
            contraint = "entropy and volume";
            NS = gas.NS;
            NE = gas.NE;
            J_SIZE = NE + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;
    }

    string ions;
    if (gas.HAS_IONS) ions = " and charge constraint enforced";
    else ions = " and no charge constraint";
    cout << "-- Configuration complete. Using " << energy << " minimization with " << contraint << " held constant" << ions << endl << endl;
}

void CESolver::compute_equilibrium(double& a, double& b) {
    (this->*equilibrium)(a, b);
}

inline void CESolver::compute_equilibrium_TV(double& T, double& V) {
    
    bool converged = false;
    static Vector J(J_SIZE * J_SIZE), F(J_SIZE), DELTA(J_SIZE), DlnNj(NS); 

    gas.V = V;
    gas.T = T;

    double sum;

    sum = 0.0;
    for (int j = 0; j < NS; ++j) 
        sum += gas.X0[j];

    for (int j = 0; j < NS; ++j) {
        gas.N[j] = max(1e-12, gas.X0[j]) / sum;
    }

    NASA_fits();
    // auto start = NOW;
    
    while (!converged) {        

        helm::compute_mu(gas);
        helm::form_elemental(J.data(), F.data(), gas);
        if (gas.HAS_IONS) helm::form_charge(J.data(), F.data(), gas);

        LUSolve(J.data(), F.data(), DELTA.data(), J_SIZE, 1);

        for (int j = 0; j < NS; ++j) {
            sum = 0.0;
            for (int i = 0; i < NE; ++i) {
                sum += gas.a[i * NS + j] * DELTA[i];
            }
            DlnNj[j] = sum - gas.mu_RT[j];
            if (gas.HAS_IONS) DlnNj[j] += gas.species[j].q * DELTA[NE];
        }
    
        for (int j = 0; j < NS; ++j) {
            gas.N[j] *= exp(DlnNj[j]);
        }
        
        // converged = check_convergence(DlnNj.data());   
    }


        sum = 0.0;
        for (int j = 0; j < NS; ++j) 
            sum += gas.N[j];

        for (int j = 0; j < NS; ++j) {
            gas.X[j] = gas.N[j] / sum;
            gas.X0[j] = gas.N[j];
        }

    // auto end = NOW;
    // auto duration = chrono::duration<double>(end - start);
    // cout << "Average iteration takes " << duration.count()/iteration << " seconds" << endl;

    // cout << "-- Equilibrium computed in " << iteration << " iterations" << endl;
}

inline void CESolver::compute_equilibrium_UV(double& U, double& V) {
    
    bool converged = false;

    gas.V = V;
    gas.uo = U;
    
    double sum;
    sum = 0.0;
    for (int j = 0; j < NS; ++j) 
        sum += gas.X0[j];

    for (int j = 0; j < NS; ++j) {
        gas.N[j] = max(1e-12, gas.X0[j]) / sum;
    }

    vector<double> J(J_SIZE * J_SIZE), F(J_SIZE), DELTA(J_SIZE), DlnNj(NS); 

    while (!converged) {   

 
        NASA_fits();

        // Compute u'
        gas.up = 0.0;
        for (int j = 0; j < NS; ++j) 
            gas.up += gas.N[j] * (gas.H0_RT[j] - 1.0) * gcon * gas.T;

        // Form matrix system
        helm::compute_mu(gas);
        helm::form_elemental(J.data(), F.data(), gas);
        if (gas.HAS_IONS) helm::form_charge(J.data(), F.data(), gas);
        helm::form_U(J.data(), F.data(), gas);

        // Solve matrix system
        LUSolve(J.data(), F.data(), DELTA.data(), J_SIZE, 1); 

        // Reform Dln(N_j)
        for (int j = 0; j < NS; ++j) {
            sum = 0.0;
            for (int i = 0; i < NE; ++i) {
                sum += gas.a[i * NS + j] * DELTA[i];
            }
            DlnNj[j] = sum + (gas.H0_RT[j] - 1.0) * DELTA[J_SIZE - 1] - gas.mu_RT[j];
            if (gas.HAS_IONS) DlnNj[j] += gas.species[j].q * DELTA[NE];
        }

        // Get new temperature
        gas.T *= exp(DELTA[J_SIZE - 1]);
    
        // Get new moles
        for (int j = 0; j < NS; ++j) {
            gas.N[j] *= exp(DlnNj[j]);
        }
        
        // Check convergence
        // converged = check_convergence(DlnNj.data());    
        if (fabs(DELTA[J_SIZE - 1]) > 1.0e-4) converged = false;    
    }

    // Get final molar concentrations
    sum = 0.0;
    for (int j = 0; j < NS; ++j) 
        sum += gas.N[j];

    for (int j = 0; j < NS; ++j) {
        gas.X[j] = gas.N[j] / sum;
        gas.X0[j] = gas.N[j];
    }
}

inline void CESolver::compute_equilibrium_SV(double& S, double& V) {
    cout << "Not made yet" << endl;
}   

inline void CESolver::compute_equilibrium_TP(double& T, double& P) {
   
    bool converged = false;
    static Vector J(J_SIZE * J_SIZE), F(J_SIZE), DELTA(J_SIZE), DlnNj(NS); 

    gas.T = T;
    gas.p = P;

    double sum, e, dlnT = 0.0;

    sum = 0.0;
    for (int j = 0; j < NS; ++j) 
        sum += max(1e-12, gas.X0[j]);

    gas.N_tot = sum;

    for (int j = 0; j < NS; ++j) {
        gas.N[j] = max(1e-6, gas.X0[j]) / sum;
    }

    NASA_fits();

    while (!converged) {        

        gibbs::compute_mu(gas);
        gibbs::form_elemental(J.data(), F.data(), gas);
        if (gas.HAS_IONS) gibbs::form_charge(J.data(), F.data(), gas);

        LUSolve(J.data(), F.data(), DELTA.data(), J_SIZE, 1);
        
        // cout << "========" << endl;
        // for (int i = 0; i < J_SIZE; ++i) {
        //     for (int j = 0; j < J_SIZE; ++j) {
        //         cout << J[i * J_SIZE + j] << "\t";
        //     }
        //     cout << "|" << F[i] << "\t|" << DELTA[i] << endl;
        // }
        // cout << "========" << endl;

        for (int j = 0; j < NS; ++j) {
            sum = 0.0;
            for (int i = 0; i < NE; ++i) {
                sum += gas.a[i * NS + j] * DELTA[i];
            }
            DlnNj[j] = DELTA[NE] + sum - gas.mu_RT[j];
            if (gas.HAS_IONS) DlnNj[j] += gas.species[j].q * DELTA[NE + 1];
        }

        e = compute_damping(DlnNj.data(), DELTA[NE], dlnT);
     
        gas.N_tot *=  exp(e * DELTA[NE]);
    
        for (int j = 0; j < NS; ++j)
            gas.N[j] *= exp(e * DlnNj[j]);

        converged = check_convergence(DlnNj.data(), DELTA[NE]);           
    }

    sum = 0.0;
    for (int j = 0; j < NS; ++j) 
        sum += gas.N[j];

    for (int j = 0; j < NS; ++j) {
        gas.X[j] = gas.N[j] / sum;
        gas.X0[j] = gas.N[j];
    }
    compute_mass_fractions();
}

inline void CESolver::compute_equilibrium_HP(double& H, double& P) {
    cout << "Not made yet" << endl;
}

inline void CESolver::compute_equilibrium_SP(double& S, double& P) {
    cout << "Not made yet" << endl;
}

inline void CESolver::findTRange() {
    if      (gas.T >= 200.0   && gas.T <= 1000.0)   T_flag = 0; 
    else if (gas.T >  1000.0  && gas.T <= 6000.0)   T_flag = 1; 
    else if (gas.T >  6000.0  && gas.T <= 20000.0)  T_flag = 2;
    else {
        T_flag = -1;
        cout << "T = " << gas.T << ", outside polynomial range" << endl;
    }
}

inline array<double, 7> CESolver::temp_base() {

    double T = gas.T;
    array<double, 7> Ts;
    double Ti = 1/T;

    Ts[0] = Ti * Ti;  // 1 / T^2
    Ts[1] = Ti;  // 1 / T
    Ts[2] = log(T);     // ln(T)
    Ts[3] = T;          // T
    Ts[4] = Ts[3] * T;  // T^2
    Ts[5] = Ts[4] * T;  // T^3
    Ts[6] = Ts[5] * T;  // T^4
    
    return Ts;
}

inline void CESolver::NASA_fits() {

    // H(T) / RT = -a0 * T^(-2) + a1 * ln(T)/T + a2 + a3 * T/2 + a4 * T^2 / 3 + a5 * T^3/4 + a6 * T^4 / 5 + a7 / T
    // S(T) / R = -a0 * T^(-2)/2 - a1 / T + a2 * ln(T) + a3 * T + a4 * T^2 / 2 + a5 * T^3/3 + a6 * T^4 / 4 + a8

    findTRange();                           // Find the temperature range to use for NASA Polynomials
    const auto Ts = temp_base();       // Calculate Temperature variables.

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

}

inline void CESolver::compute_mass_fractions() {

        gas.MW = 0.0;
        for (int i = 0; i < NS; ++i) {
                gas.MW += gas.X[i] * gas.species[i].mw;
        }

        for (int i = 0; i < NS; ++i) {
                gas.Y[i] = (gas.X[i] * gas.species[i].mw) / gas.MW;
        }
}

void CESolver::print_properties() {


    cout << endl << setw(30) << "Mixture Temperature: " << gas.T << " [K]" << endl;
    cout << setw(30) << "Mixture Pressure: " << gas.p / 1000.0 << " [kPa]" << endl;
    // cout << setw(30) << "Mixture Density: " << gas.rho << " [kg/m^3]" << endl;
    // cout << setw(30) << "Mixture Internal Energy: " << gas.e / 1000.0 << " [kJ/kg]" << endl;
    // cout << setw(30) << "Mixture Gas Constant: " << gas.R << " [J/kg-K]" << endl;
    // cout << setw(30) << "Mixture cp: " << gas.cp << " [J/kg-K]" << endl;
    // cout << setw(30) << "Mixture cv: " << gas.cv << " [J/kg-K]" << endl;
    // cout << setw(30) << "Mixture gamma: " << gas.gamma << endl;
    // cout << setw(30) << "Mixture MW: " << gas.MW << " [g/mol]" << endl;


    cout << endl << setw(67) << "====: Species Molar Fractions :==== " << endl;
    for (int i = 0; i < NS; ++i) {
        cout << setw(10) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.X[i] << "\t";
        if ((i + 1) % 4 == 0) cout << endl;
    }

    // cout << endl << endl << setw(67) <<  "====: Species Mass Fractions :==== " << endl;
    // for (int i = 0; i < NS; ++i) {
    //     cout << setw(10) << gas.species[i].name << ": " << fixed << setprecision(4) << gas.Y[i] << "\t";
    //     if ((i + 1) % 4 == 0) cout << endl;
    // }
}

inline bool CESolver::check_convergence(double* dlnj, double& dln) {

    double sum = 0.0, check, tol = 0.5e-5;

    for (int j = 0; j < NS; ++j) 
        sum += gas.N[j];
    

    for (int j = 0; j < NS; ++j) {
        check = gas.N[j] * fabs(dlnj[j]) / sum;
        if (check > tol) return false;
    }


    check = gas.N_tot * fabs(dln) / sum;
    if (check > tol) return false;
    
    return true;
}

inline double CESolver::compute_damping(double* dlnj, double& dlnN, double& dlnT) {
    // Constants from the text
    constexpr double SIZE = 18.420681;          // ~ ln(1e8)
    constexpr double C1   =  9.2103404;         // ~ ln(1e4)

    // ---- λ1 = 2 / max( 5|ΔlnT|, 5|ΔlnN|, max_j |Δln n_j| ) ----
    double max_dlnj = 0.0;
    for (int j = 0; j < NS; ++j) {               // gaseous species only
        double a = fabs(dlnj[j]);
        if (a > max_dlnj) max_dlnj = a;
    }
    double denom1 = max( 5.0 * abs(dlnT), max(5.0 * abs(dlnN), max_dlnj));
    double lambda1 = (denom1 > 0.0) ? std::min(1.0, 2.0 / denom1) : 1.0;

    // ---- λ2 from (3.2) for species with ln(nj/n) <= -SIZE and Δln nj >= 0 ----
    double lambda2 = 1.0;     // no constraint unless triggered
    bool triggered = false;


    for (int j = 0; j < NS; ++j) {

        const double ratio = log(gas.N[j] / gas.N_tot);  // ln(nj/n)

        if (ratio <= -SIZE && dlnj[j] >= 0.0) {

            double denom2 = dlnj[j] - dlnN;

            if (abs(denom2) > 0.0) {

                double frac = abs((-ratio - C1) / denom2);
                lambda2 = min(lambda2, frac);
                triggered = true;

            }
        }
    }

    if (!triggered) lambda2 = 1.0;

    // ---- Final damping factor λ = min(1, λ1, λ2) ----
    double lambda = min(1.0, min(lambda1, lambda2));
    if (lambda < 0.0) lambda = 0.0;  // safety clamp

    return lambda;
}

















