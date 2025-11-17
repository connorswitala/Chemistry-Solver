#include "CESolver.h"

// Constructor function
CESolver::CESolver(mix& gas_in, ConstraintType constrainttype) : gas(gas_in) {

    // For terminial output
    string energy;
    string contraint;

    // Create Newton solver Matrix/Vector sizes and set booleans/function pointers.
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

        case::ConstraintType::CFD:
            gas.NEEDS_T = true;
            energy = "Helmholtz";
            contraint = "internal energy and density";
            NS = gas.NS;
            NE = gas.NE;
            J_SIZE = NE + gas.HAS_IONS + gas.NEEDS_T;
            gas.J_SIZE = J_SIZE;
            break;
    }

    J = vector<double>(J_SIZE * J_SIZE);
    F = vector<double>(J_SIZE);
    DELTA = vector<double>(J_SIZE);

    string ions;
    if (gas.HAS_IONS) ions = " and charge constraint enforced";
    else ions = " and no charge constraint";
    cout << "-- Chemical equilibrium configuration complete." << endl;
    cout << "-- Using " << energy << " minimization with " << contraint << " held constant" << ions << endl << endl;
}

void CESolver::CFD_equilibrium(double& e, double& rho) {
    bool converged = false;

    gas.rho = rho;
    gas.uo = e;

    double eps;     // Damping coefficient
    double rlntot;  // Total ln(N)

    int iteration = 0;
    int maxiter = 200;
    vector<double> DlnNj(NS);

    // Initial conditions
    gas.N_tot = 0.0;
    for (int j = 0; j < NS; ++j) {
        gas.N[j] = gas.X0[j];
        gas.N_tot += gas.N[j];
    }

    // Being iterations
    while (iteration < maxiter) {

        NASA_fits();    // Call NASA fits to update temperature dependant variables

        // Compute u'
        gas.up = 0.0;
        for (int j = 0; j < NS; ++j) {
            gas.up += gas.N[j] * gas.U0_RT[j];
        }

        helm::compute_mu(gas);                          // Compute chemical potential
        helm::form_elemental(J.data(), F.data(), gas);  // Form elemental constrains part of Jacobian and RHS
        if (gas.HAS_IONS) 
            helm::form_charge(J.data(), F.data(), gas); // If ions present, form charge constraint section
        helm::form_U(J.data(), F.data(), gas);          // Form ln(T) row

        LUSolve(J.data(), F.data(), DELTA.data(), J_SIZE, 1);   // Solve matrix system

        rlntot = 0.0; 

        // Reform Dln(N_j)
        for (int j = 0; j < NS; ++j) {
            DlnNj[j] =  gas.U0_RT[j] * DELTA[J_SIZE - 1] - gas.mu_RT[j];
            for (int i = 0; i < NE; ++i) {
                DlnNj[j] += gas.a[i * NS + j] * DELTA[i];
            }
            if (gas.HAS_IONS) DlnNj[j] += gas.species[j].q * DELTA[NE];
            rlntot += gas.N[j] * DlnNj[j];
        }

        rlntot /= gas.N_tot;

        eps = compute_damping(DlnNj, rlntot, DELTA[J_SIZE - 1]);  // 

        // Update N[j]
        gas.N_tot = 0.0;
        for (int j = 0; j < NS; ++j) {
            gas.N[j] = max(gas.N[j] * exp(eps *DlnNj[j]), 1e-18);
            gas.N_tot += gas.N[j];
        }

        // Get new temperature
        gas.T *= exp(eps * DELTA[J_SIZE - 1]);

        iteration++;

        // Check convergence
        converged = check_convergence(DlnNj.data(), rlntot, DELTA[J_SIZE - 1]);
        if (converged) iteration = maxiter;
    }

    // Get final molar concentrations
    for (int j = 0; j < NS; ++j) {
        gas.X[j] = gas.N[j] / gas.N_tot;
        gas.X0[j] = gas.N[j];
    }

    compute_mixture_properties();   // Compute molecular weights 
    compute_derivativesCFD();          // Find important derivatives
}

// Universal function caller
void CESolver::compute_equilibrium(double a, double b) {
    (this->*equilibrium)(a, b);
}


// ===== Mostly the same as CFD_equilibrium(), not commented =====
inline void CESolver::compute_equilibrium_TV(double T, double V) {

    bool converged = false;
    static Vector DlnNj(NS);

    gas.V = V;
    gas.T = T;
    gas.rho = 1 / V;

    double sum, e;
    double dlnt = 0.0;
    double dlnn = 0.0;

    gas.N_tot = 0.0;
    for (int j = 0; j < NS; ++j) {
        gas.N[j] = gas.X0[j];
        gas.N_tot += gas.N[j];
    }

    NASA_fits();

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

        e = compute_damping(DlnNj, dlnn, dlnt);

        for (int j = 0; j < NS; ++j) {
            gas.N[j] *= exp(e * DlnNj[j]);
        }

        converged = check_convergence(DlnNj.data(), dlnn, dlnt);
    }

        gas.N_tot = 0.0;
        for (int j = 0; j < NS; ++j)
            gas.N_tot += gas.N[j];

        for (int j = 0; j < NS; ++j) {
            gas.X[j] = gas.N[j] / gas.N_tot;
            gas.X0[j] = gas.N[j];
        }
        compute_mixture_properties();
        compute_derivatives();
}

inline void CESolver::compute_equilibrium_UV(double U, double V) {

    bool converged = false;

    gas.V = V;
    gas.rho = 1 / V;
    gas.uo = U;

    double sum;
    double dlnn = 0.0;
    double e, rlntot;

    int iteration = 0;
    int maxiter = 200;
    vector<double> DlnNj(NS);

    // Initial conditions
    gas.N_tot = 0.0;
    for (int j = 0; j < NS; ++j) {
        gas.N[j] = gas.X0[j];
        gas.N_tot += gas.N[j];
    }

    while (iteration < maxiter) {

        NASA_fits();

        // Compute u'
        gas.up = 0.0;
        // gas.e_ref = 0.0;
        for (int j = 0; j < NS; ++j) {
            gas.up += gas.N[j] * gas.U0_RT[j];
            gas.e_ref += gas.N[j] * gas.species[j].href;
        }

        // Form matrix system
        helm::compute_mu(gas);
        helm::form_elemental(J.data(), F.data(), gas);
        if (gas.HAS_IONS) helm::form_charge(J.data(), F.data(), gas);
        helm::form_U(J.data(), F.data(), gas);

        // Solve matrix system
        LUSolve(J.data(), F.data(), DELTA.data(), J_SIZE, 1);

        rlntot = 0.0;

        // Reform Dln(N_j)
        for (int j = 0; j < NS; ++j) {
            DlnNj[j] =  gas.U0_RT[j] * DELTA[J_SIZE - 1] - gas.mu_RT[j];
            for (int i = 0; i < NE; ++i) {
                DlnNj[j] += gas.a[i * NS + j] * DELTA[i];
            }
            if (gas.HAS_IONS) DlnNj[j] += gas.species[j].q * DELTA[NE];
            rlntot += gas.N[j] * DlnNj[j];
        }

        rlntot /= gas.N_tot;
        e = compute_damping(DlnNj, rlntot, DELTA[J_SIZE - 1]);

        gas.N_tot = 0.0;
        for (int j = 0; j < NS; ++j) {
            gas.N[j] = max(gas.N[j] * exp(e *DlnNj[j]), 1e-18);
            gas.N_tot += gas.N[j];
        }

        // Get new temperature
        gas.T *= exp(e * DELTA[J_SIZE - 1]);

        iteration++;

        // Check convergence
        converged = check_convergence(DlnNj.data(), dlnn, DELTA[J_SIZE - 1]);
        if (fabs(DELTA[J_SIZE - 1]) > 1.0e-5) converged = false;
        if (converged) iteration = maxiter;        
    }

    // Get final molar concentrations
    for (int j = 0; j < NS; ++j) {
        gas.X[j] = gas.N[j] / gas.N_tot;
        gas.X0[j] = gas.N[j];
    }

    // gas.up = 0.0;
    // for (int j = 0; j < NS; ++j) {
    //     gas.up += gas.N[j] * (gas.U0_RT[j] * gcon * gas.T + gas.species[j].href);
    // }

    compute_mixture_properties();
    compute_derivativesCFD();

}

inline void CESolver::compute_equilibrium_SV(double S, double V) {
    cout << "Not made yet" << endl;
}

inline void CESolver::compute_equilibrium_TP(double T, double P) {

    bool converged = false;
    static Vector DlnNj(NS);

    gas.T = T;
    gas.p = P;

    double sum, e, dlnT = 0.0;
    int iteration = 0;

    // Initial conditions
    gas.N_tot = 0.0;
    for (int j = 0; j < NS; ++j) {
        gas.N[j] = max(gas.X0[j], 1e-12);
        gas.N_tot += gas.N[j];
    }

    NASA_fits();

    while (!converged) {

        gibbs::compute_mu(gas);
        gibbs::form_elemental(J.data(), F.data(), gas);
        if (gas.HAS_IONS) gibbs::form_charge(J.data(), F.data(), gas);

        LUSolve(J.data(), F.data(), DELTA.data(), J_SIZE, 1);

        for (int j = 0; j < NS; ++j) {
            sum = 0.0;
            for (int i = 0; i < NE; ++i) {
                sum += gas.a[i * NS + j] * DELTA[i];
            }
            DlnNj[j] = DELTA[NE] + sum - gas.mu_RT[j];
            if (gas.HAS_IONS) DlnNj[j] += gas.species[j].q * DELTA[NE + 1];
        }

        e = compute_damping(DlnNj, DELTA[NE], dlnT);

        gas.N_tot *= exp(e * DELTA[NE]);

        for (int j = 0; j < NS; ++j)
            gas.N[j] *= exp(e * DlnNj[j]);

        converged = check_convergence(DlnNj.data(), DELTA[NE], dlnT);
    }

    for (int j = 0; j < NS; ++j) {
        gas.X[j] = gas.N[j] / gas.N_tot;
        gas.X0[j] = gas.N[j];
    }

    compute_mixture_properties();
    gas.rho = gas.p / (gas.R * gas.T);
    gas.V = 1.0 / gas.rho;
    compute_derivatives();

    // gas.up = 0.0;
    // for (int j = 0; j < NS; ++j) {
    //     gas.up += gas.N[j] * (gas.U0_RT[j] * gcon * gas.T + gas.species[j].href);
    // }    
}

inline void CESolver::compute_equilibrium_HP(double H, double P) {
    cout << "Not made yet" << endl;
}

inline void CESolver::compute_equilibrium_SP(double S, double P) {
    cout << "Not made yet" << endl;
}

// ===== Helper functions =====

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


    findTRange();                   // Find the temperature range to use for NASA Polynomials
    const auto Ts = temp_base();    // Calculate temperature variables.

    for (int j = 0; j < NS; ++j) {

        const auto& poly = gas.species[j].poly;
        const double* coeff = &poly.at(T_flag * NCOEF);

        // Standard state enthlapy
        gas.H0_RT[j] = - coeff[0] * Ts[0]
                   + coeff[1] * Ts[1] * Ts[2]
                   + coeff[2]
                   + 0.5 * coeff[3] * Ts[3]
                   + coeff[4] * Ts[4] / 3.0
                   + 0.25 * coeff[5] * Ts[5]
                   + 0.2 * coeff[6] * Ts[6]
                   + coeff[7] * Ts[1];                  

        // Standard state entropy                   
        gas.S0_R[j] = - 0.5 * coeff[0] * Ts[0]
                   - coeff[1] * Ts[1]
                   + coeff[2] * Ts[2]
                   + coeff[3] * Ts[3]
                   + 0.5 * coeff[4] * Ts[4]
                   + coeff[5] * Ts[5] / 3.0
                   + 0.25 * coeff[6] * Ts[6]
                   + coeff[8];

        // Standard state specific heat
        gas.CP0_R[j] = coeff[0] * Ts[0]
                     + coeff[1] * Ts[1]
                     + coeff[2]
                     + coeff[3] * Ts[3]
                     + coeff[4] * Ts[4]
                     + coeff[5] * Ts[5]
                     + coeff[6] * Ts[6];

        // Standard state internal energy
        gas.U0_RT[j] = gas.H0_RT[j] - 1.0;

        // Standard state chemical potential
        gas.mu0_RT[j] = gas.H0_RT[j] - gas.S0_R[j];
    }

}

// Compute molecular weights, gas constant, mass fractions
inline void CESolver::compute_mixture_properties() {
    gas.MW = 0.0;
    for (int i = 0; i < NS; ++i) {
            gas.MW += gas.X[i] * gas.species[i].mw;
    }

    gas.R = gcon * gas.N_tot;
    for (int i = 0; i < NS; ++i) {
            gas.Y[i] = (gas.X[i] * gas.species[i].mw) / gas.MW;
    }
}

// Convergence checker
inline bool CESolver::check_convergence(double* dlnj, double& dln, double& dlnt) {

    double sum = 0.0, check, tol = 0.5e-5;

    for (int j = 0; j < NS; ++j)
        sum += gas.N[j];

    for (int j = 0; j < NS; ++j) {
        check = gas.N[j] * fabs(dlnj[j]) / sum;
        if (check > tol) return false;
    }

    if (fabs(dlnt) > 1.0e-5) 
        return false;

    check = gas.N_tot * fabs(dln) / sum;
    if (check > tol) return false;

    return true;
}

// Compute damping coefficient for update to solution variables
inline double CESolver::compute_damping(vector<double>& DlnNj, double& dlnN, double& dlnT) {

    double lam1 = 5.0 * max(max( fabs(*max_element(DlnNj.begin(), DlnNj.end())), fabs(dlnN)), dlnT);
    double lam2 = 1.0;
    double ratio;

    for (int j = 0 ; j < NS; ++j) {
        ratio = log(gas.N[j] / gas.N_tot);
        if (ratio > -18.420681) 
            lam1 = max(lam1, fabs(DlnNj[j]));
        if (ratio < -18.420681 && DlnNj[j] > 0.0) 
            lam2 = min( lam2, fabs((-DlnNj[j] - 9.2103404)/(fabs(DlnNj[j] - dlnN) + 1.0e-16)) );
    }
    lam1 = 2.0 / (lam1 + 1.0e-16);
    return min(1.0, min(lam1, lam2));
}

// Compute derivatices for A_rho Jacobian
inline void CESolver::compute_derivativesCFD() {

    int size = NE + 1;
    vector<double> A(size * size, 0.0), RHS(size, 0.0), x(size, 0.0);
    int offset, idx = NE + gas.HAS_IONS;


    vector<double> an(NE, 0.0);

    double sum = 0.0, sum1 = 0.0;

    // ===== Derivatives wrt T =====
    
    for (int k = 0; k < NE; ++k) {

        offset = k * size;

        for (int i = 0; i < NE; ++i) {
            for (int j = 0; j < NS; ++j)
                A[offset + i] = J[k * J_SIZE + i];
        }

        for (int j = 0; j < NS; ++j) {
            an[k] += gas.a[k * NS + j] * gas.N[j];       
        }

        A[offset + NE] = an[k];
        RHS[k] = -(J[k * J_SIZE + idx] + an[k]);
    }

    for (int i = 0; i < NE; ++i) 
        A[NE * size + i] = an[i];

    for (int j = 0; j < NS; ++j)
        RHS[size - 1] -= gas.N[j] * gas.H0_RT[j];


    LUSolve(A.data(), RHS.data(), x.data(), size, 1);

    double dlndlt = x[size - 1];

    // ===== Compute cp =====
    gas.cp = 0.0;

    for (int j = 0; j < NS; ++j) 
        gas.cp += gas.N[j] * gas.CP0_R[j];


    for (int i = 0; i < NE; ++i) {        
        gas.cp += (J[i * J_SIZE + idx] + an[i]) * x[i];        
    }

    for (int j = 0; j < NS; ++j)
        gas.cp += gas.N[j] * gas.H0_RT[j] * (dlndlt + gas.H0_RT[j]);

    // ===== Derivatives wrt P =====
    
    for (int k = 0; k < NE; ++k) {
        RHS[k] = an[k];
    }   

    RHS[size - 1] = 0.0;
    for (int j = 0; j < NS; ++j)
        RHS[size - 1] += gas.N[j];

    LUSolve(A.data(), RHS.data(), x.data(), size, 1);

    double dlvdlp = x[size - 1] - 1.0;
    double dlvdlt = 1.0 + dlndlt;

    gas.p = gas.rho * gas.N_tot * gcon * gas.T;
    gas.cv = gas.cp + (gas.p * dlvdlt * dlvdlt) / (gas.rho * gas.T * dlvdlp) / gcon;
    gas.cp *= gcon;
    gas.cv *= gcon;
    gas.gamma = gas.cp / gas.cv;
    double gammas = -gas.gamma / dlvdlp;

    double deno = dlvdlt * dlvdlt + gas.cp * dlvdlp / gas.R;

    gas.dpdr = gas.p / gas.rho * (dlvdlt - gas.cp/gas.R) / deno;
    gas.dtdr = -gas.T / gas.rho * (dlvdlt + dlvdlp) / deno;
    gas.c = sqrt(gas.N_tot * gcon * gas.T * gammas);
}

// Compute derivatices for A_rho Jacobian
inline void CESolver::compute_derivatives() {

    int size = NE + 1;
    vector<double> A(size * size, 0.0), RHS(size, 0.0), x(size, 0.0);
    int offset, idx = NE + gas.HAS_IONS;

    vector<double> an(NE, 0.0);
    vector<double> anh(NE, 0.0);

    double sum = 0.0, sum1 = 0.0;
    double a;

    // ===== Derivatives wrt T =====
    
    for (int k = 0; k < NE; ++k) {

        offset = k * size;

        for (int i = 0; i < NE; ++i) {
            for (int j = 0; j < NS; ++j)
                A[offset + i] = J[k * J_SIZE + i];
        }

        for (int j = 0; j < NS; ++j) {
            a = gas.a[k * NS + j] * gas.N[j];
            an[k] += a;
            anh[k] += a * gas.H0_RT[j];
        }

        A[offset + NE] = an[k];
        RHS[k] = -anh[k];
    }

    for (int i = 0; i < NE; ++i) 
        A[NE * size + i] = an[i];

    for (int j = 0; j < NS; ++j)
        RHS[size - 1] -= gas.N[j] * gas.H0_RT[j];


    LUSolve(A.data(), RHS.data(), x.data(), size, 1);

    double dlndlt = x[size - 1];

    // ===== Compute cp =====
    gas.cp = 0.0;

    for (int j = 0; j < NS; ++j) 
        gas.cp += gas.N[j] * gas.CP0_R[j];


    for (int i = 0; i < NE; ++i)     
        gas.cp += anh[i] * x[i];        
    

    for (int j = 0; j < NS; ++j)
        gas.cp += gas.N[j] * gas.H0_RT[j] * (dlndlt + gas.H0_RT[j]);

    // ===== Derivatives wrt P =====
    
    for (int k = 0; k < NE; ++k) {
        RHS[k] = an[k];
    }   

    RHS[size - 1] = 0.0;
    for (int j = 0; j < NS; ++j)
        RHS[size - 1] += gas.N[j];

    LUSolve(A.data(), RHS.data(), x.data(), size, 1);

    double dlvdlp = x[size - 1] - 1.0;
    double dlvdlt = 1.0 + dlndlt;

    gas.p = gas.rho * gas.N_tot * gcon * gas.T;
    gas.cv = gas.cp + (gas.p * dlvdlt * dlvdlt) / (gas.rho * gas.T * dlvdlp) / gcon;
    gas.cp *= gcon;
    gas.cv *= gcon;
    gas.gamma = gas.cp / gas.cv;
    double gammas = -gas.gamma / dlvdlp;

    double deno = dlvdlt * dlvdlt + gas.cp * dlvdlp / gas.R;

    gas.dpdr = gas.p / gas.rho * (dlvdlt - gas.cp/gas.R) / deno;
    gas.dtdr = -gas.T / gas.rho * (dlvdlt + dlvdlp) / deno;
    gas.c = sqrt(gas.N_tot * gcon * gas.T * gammas);
}













