#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "thermoObjs.h"

using namespace std;

// Namespace for Gibbs minimization functions
namespace gibbs {

    // Compute chemical potential
    inline void compute_mu(mix& gas) {
        for (int j = 0; j < gas.NS; ++j) {
            gas.mu_RT[j] = gas.mu0_RT[j] + log(gas.N[j] / gas.N_tot) + log(gas.p / 1.0e5);
        }
    }        

    // This function forms the entries in J and F for the elemental contraints. 
    // It is called for every minimization procedure
    inline void form_elemental(double* J, double* F, mix& gas) {

        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        int offset;

        vector<double> b(NE, 0.0);

        // ===== Elemental rows in elemental columns [NE * NE]  =====
        for (int k = 0; k < NE; ++k) {

            offset = k * J_SIZE;
            
            // Used in both J and F
            b[k] = 0.0;
            for (int j = 0; j < NS; ++j) {
                b[k] += gas.a[k * NS + j] * gas.N[j];
            }

            // Elemental columns in elemental rows
            for (int i = 0; i < NE; ++i) {

                J[offset + i] = 0.0;
                for (int j = 0; j < NS; ++j)
                    J[offset + i] += gas.a[k * NS + j] * gas.a[i * NS + j] * gas.N[j];
            }

            // ln(N) column in elemental rows
            J[offset + NE] = b[k];


            // RHS
            F[k] = gas.b[k] - b[k];
            for (int j = 0; j < NS; ++j) 
                F[k] += gas.a[k * NS + j] * gas.N[j] * gas.mu_RT[j]; 
        }

        // ====== ln(N) row ======
        offset = NE * J_SIZE;

        // Elemental columns in ln(N) row
        for (int i = 0; i < NE; ++i)
            J[offset + i] = b[i];
        

        // ln(N) column in ln(N) row
        J[offset + NE] = -gas.N_tot; 
        for (int j = 0; j < NS; ++j) 
            J[offset + NE] += gas.N[j];
        
        // RHS
        F[NE] = -J[offset + NE];  
        for (int j = 0; j < NS; ++j)
            F[NE] += gas.N[j] * gas.mu_RT[j];
    }

    // This functions forms the entries in J and F for the charge constraint. 
    // It is called if the boolean gas.HAS_IONS is true
    inline void form_charge(double* J, double* F, mix& gas) {
        
        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        int offset = (NE + 1) * J_SIZE;

        vector<double> qa(NE, 0.0);
        double qnsum = 0.0;

        // ===== Elemental rows =====
        for (int k = 0; k < NE; ++k) {

            // Charge column in elemental rows
            for (int j = 0; j < NS; ++j) 
                qa[k] += gas.a[k * NS + j] * gas.species[j].q * gas.N[j];

            J[k * J_SIZE + NE + 1] = qa[k];
        }

        for (int j = 0; j < NS; ++j) 
            qnsum += gas.species[j].q * gas.N[j];

        // Charge column in ln(N) row
        J[NE * J_SIZE + NE + 1] = qnsum;
        
        // Elemental columns in charge row
        for (int i = 0; i < NE; ++i) 
            J[offset + i] = qa[i];
        
        // ln(N) column in charge row
        J[offset + NE] = qnsum;

        // Charge column in charge row;
        J[offset + NE + 1] = 0.0;
        for (int j = 0; j < NS; ++j) 
            J[offset + NE + 1] += gas.species[j].q * gas.species[j].q * gas.N[j];

        // RHS
        F[NE + 1] = -qnsum;
        for (int j = 0; j < NS; ++j) 
            F[NE + 1] += gas.species[j].q * gas.N[j] * gas.mu_RT[j];
    }

    // void form_H(double* J, double* F, mix& gas) {

    // }

    // void form_S(double* J, double* F, mix& gas) {

    // }

}

// Namespace for Helmholtz minimization functions
namespace helm {

    // Compute chemical potential
    inline void compute_mu(mix& gas) {
        double pb = gcon * gas.rho * gas.T * (1.0e-5);

        for (int j = 0; j < gas.NS; ++j) {
            gas.mu_RT[j] = gas.mu0_RT[j] + log(gas.N[j] * pb);
        }
    }

    // This function forms the entries in J and F for the elemental contraints.
    inline void form_elemental(double* J, double* F, mix& gas) {

        int NS = gas.NS;
        int NE = gas.NE; 
        int J_SIZE = gas.J_SIZE;
        
        double sum;
        int offset;

        // ===== Elemental rows =====
        for (int k = 0; k < NE; ++k) {

            offset = k * J_SIZE;

            sum = 0.0;
            for (int j = 0; j < NS; ++j) {                
                sum += gas.a[k * NS + j] * gas.N[j];
            }

            // Elemental columns in elemental rows
            for (int i = 0; i < NE; ++i) {                
                J[offset + i] = 0.0;
                for (int j = 0; j < NS; ++j)
                    J[offset + i] += gas.N[j] * gas.a[i * NS + j] * gas.a[k * NS + j];            
            }            

            // RHS
            F[k] = gas.b[k] - sum;
            for (int j = 0; j < NS; ++j) {
                F[k] += gas.N[j] * gas.a[k * NS + j] * gas.mu_RT[j];
            }
        }
    }

    // This function forms the entries in J and F for the charge contraint.
    inline void form_charge(double* J, double* F, mix& gas) {
        
        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        vector<double> qa(NE);

        // ===== Elemental rows =====
        for (int k = 0; k < NE; ++k) {
            qa[k] = 0.0;

            // Charge columns in elemental rows
            for (int j = 0; j < NS; ++j) 
                qa[k] += gas.species[j].q * gas.a[k * NS + j] * gas.N[j];

            J[k * J_SIZE + NE] = qa[k];    
        }

        // Elemental contraint columns in charge row.
        for (int i = 0; i < NE; ++i) 
            J[NE * J_SIZE + i] = qa[i];
        

        // Charge constraint column in charge row.
        J[NE * J_SIZE + NE] = 0.0;
        for (int j = 0; j < NS; ++j) 
            J[NE * J_SIZE + NE] += gas.species[j].q * gas.species[j].q * gas.N[j];

        // RHS
        F[NE] = 0.0;
        for (int j = 0; j < NS; ++j) 
            F[NE] += gas.species[j].q * gas.N[j] * (gas.mu_RT[j] - 1.0);  
    }

    // This function forms the entries in J and F for the U contraint.
    inline void form_U(double* J, double* F, mix& gas) {

        int NE = gas.NE;
        int NS = gas.NS;
        int J_SIZE = gas.J_SIZE;

        int idx = NE + gas.HAS_IONS;
        double sum, eref = 0.0;;

       // ===== Elemental rows
        for (int k = 0; k < NE; ++k) {
            sum = 0.0;
            for (int j = 0; j < NS; ++j) 
                sum += gas.a[k * NS + j] * gas.N[j] * gas.U0_RT[j];

            J[k * J_SIZE + idx] = sum;  // ln(T) columns in elemental rows
            J[idx * J_SIZE + k] = sum;  // elemental columns in ln(T) row
        }

        // ln(T) term in charge row and charge constraint in ln(T) row
        if (gas.HAS_IONS) {

            sum = 0.0;
            for (int j = 0; j < NS; ++j) {
                sum += gas.species[j].q * gas.N[j] * gas.U0_RT[j];
            }  

            J[NE * J_SIZE + NE + 1] = sum;          // Charge row
            J[(J_SIZE - 1) * J_SIZE + NE] = sum;    // ln(T) row
        }

        // ln(T) column in ln(T) row
        J[J_SIZE * J_SIZE - 1] = 0.0;    
        for (int j = 0; j < NS; ++j) 
            J[J_SIZE * J_SIZE - 1] += gas.N[j] * (gas.U0_RT[j] * gas.U0_RT[j] + gas.CP0_R[j] - 1.0) ;
        
        // RHS
        F[J_SIZE - 1] = 0.0;
        for (int j = 0; j < NS; ++j) 
            F[J_SIZE - 1] += gas.N[j] * gas.U0_RT[j] * gas.mu_RT[j];

        for (int j = 0; j < NS; ++j) 
            eref += gas.N[j] * gas.species[j].href;        

        F[J_SIZE - 1] += (gas.uo - eref) / (gcon * gas.T) - gas.up;
    }


    // inline void form_S(double* J, double* F, mix& gas) {

    // }

}



#endif