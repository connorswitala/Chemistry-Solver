#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "thermoObjs.h"

using namespace std;

namespace gibbs {

    inline void compute_mu(mix& gas) {
        for (int j = 0; j < gas.NS; ++j) {
            gas.mu_RT[j] = gas.mu0_RT[j] + log(gas.N[j] / gas.N_tot) + log(gas.p / 1.0e5);
        }
    }        

    // This function forms the entries in J and F for the elemental contraints.
    // ==== COMPLETED ====
    inline void form_elemental(double* J, double* F, mix& gas) {

        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        int offset;

        vector<double> b(NE, 0.0);

        // === Elemental rows ===
        for (int k = 0; k < NE; ++k) {

            offset = k * J_SIZE;
            
            b[k] = 0.0;
            for (int j = 0; j < NS; ++j) {
                b[k] += gas.a[k * NS + j] * gas.N[j];
            }

            for (int i = 0; i < NE; ++i) {

                J[offset + i] = 0.0;
                for (int j = 0; j < NS; ++j)
                    J[offset + i] += gas.a[k * NS + j] * gas.a[i * NS + j] * gas.N[j];
            }

            // N column
            J[offset + NE] = b[k];


            // RHS
            F[k] = gas.b[k] - b[k];
            for (int j = 0; j < NS; ++j) 
                F[k] += gas.a[k * NS + j] * gas.N[j] * gas.mu_RT[j]; 
        }

        // === N row ===
        offset = NE * J_SIZE;

        // Elemental columns in N row
        for (int i = 0; i < NE; ++i)
            J[offset + i] = b[i];
        

        // N column in N row
        J[offset + NE] = -gas.N_tot; 
        for (int j = 0; j < NS; ++j) 
            J[offset + NE] += gas.N[j];
        
        // RHS

        F[NE] = -J[offset + NE];  
        for (int j = 0; j < NS; ++j)
            F[NE] += gas.N[j] * gas.mu_RT[j];
    }

    inline void form_charge(double* J, double* F, mix& gas) {
        
        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        int offset = (NE + 1) * J_SIZE;

        vector<double> qa(NE, 0.0);
        double qnsum = 0.0;

        // Charge column in elemental rows
        for (int k = 0; k < NE; ++k) {

            for (int j = 0; j < NS; ++j) 
                qa[k] += gas.a[k * NS + j] * gas.species[j].q * gas.N[j];

            J[k * J_SIZE + NE + 1] = qa[k];
        }

        for (int j = 0; j < NS; ++j) 
            qnsum += gas.species[j].q * gas.N[j];

        // Charge column in N row
        J[NE * J_SIZE + NE + 1] = qnsum;
        
        // Elemental columns in charge row
        for (int i = 0; i < NE; ++i) 
            J[offset + i] = qa[i];
        

        // N column in charge row
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


namespace helm {

    inline void compute_mu(mix& gas) {
        double pb = gcon * gas.rho * gas.T * (1.0e-5);

        for (int j = 0; j < gas.NS; ++j) {
            gas.mu_RT[j] = gas.mu0_RT[j] + log(gas.N[j] * pb);
        }
    }

    // This function forms the entries in J and F for the elemental contraints.
    // ==== COMPLETED ====
    inline void form_elemental(double* J, double* F, mix& gas) {

        int NS = gas.NS;
        int NE = gas.NE; 
        int J_SIZE = gas.J_SIZE;
        
        double sum;
        int offset;

        for (int k = 0; k < NE; ++k) {

            offset = k * J_SIZE;

            sum = 0.0;
            for (int j = 0; j < NS; ++j) {                
                sum += gas.a[k * NS + j] * gas.N[j];
            }

            for (int i = 0; i < NE; ++i) {
                
                J[offset + i] = 0.0;
                for (int j = 0; j < NS; ++j)
                    J[offset + i] += gas.N[j] * gas.a[i * NS + j] * gas.a[k * NS + j];            
            }
            

            F[k] = gas.b[k] - sum;
            for (int j = 0; j < NS; ++j) {
                F[k] += gas.N[j] * gas.a[k * NS + j] * gas.mu_RT[j];
            }
        }
    }

    // This function forms the entries in J and F for the charge contraint.
    // It currently uses one if statement to check if it needs a T column.
    // ==== COMPLETED ====
    inline void form_charge(double* J, double* F, mix& gas) {
        
        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        vector<double> qa(NE);

        // Charge column in elemental constraint rows
        for (int k = 0; k < NE; ++k) {
            qa[k] = 0.0;
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
    // It currently uses one if statement to check if it needs a chage column.
    // ==== COMPLETED ====
    inline void form_U(double* J, double* F, mix& gas) {

        int NE = gas.NE;
        int NS = gas.NS;
        int J_SIZE = gas.J_SIZE;

        int idx = NE + gas.HAS_IONS;
        double sum, eref = 0.0;;

       // Elemental rows/columns
        for (int k = 0; k < NE; ++k) {
            sum = 0.0;
            for (int j = 0; j < NS; ++j) 
                sum += gas.a[k * NS + j] * gas.N[j] * gas.U0_RT[j];

            J[k * J_SIZE + idx] = sum;  // ln(T) columns in elemental rows
            J[idx * J_SIZE + k] = sum;  // elemental columns in ln(T) row
        }

        // ln(T) term in charge constraint row and charge constrain in ln(T) row
        if (gas.HAS_IONS) {

            sum = 0.0;
            for (int j = 0; j < NS; ++j) {
                sum += gas.species[j].q * gas.N[j] * gas.U0_RT[j];
            }  

            J[NE * J_SIZE + NE + 1] = sum;
            J[(J_SIZE - 1) * J_SIZE + NE] = sum;
        }

        // ln(T) column in ln(T) row
        J[J_SIZE * J_SIZE - 1] = 0.0;    
        for (int j = 0; j < NS; ++j) 
            J[J_SIZE * J_SIZE - 1] += gas.N[j] * (gas.U0_RT[j] * gas.U0_RT[j] + gas.CP0_R[j] - 1.0) ;
        

        F[J_SIZE - 1] = 0.0;
        for (int j = 0; j < NS; ++j) 
            F[J_SIZE - 1] += gas.N[j] * gas.U0_RT[j] * gas.mu_RT[j];

        for (int j = 0; j < NS; ++j) 
            eref += gas.N[j] * gas.species[j].href;        

        F[J_SIZE - 1] += (gas.uo - eref) / (gcon * gas.T) - gas.up;
    }


    inline void form_S(double* J, double* F, mix& gas) {

    }

}



#endif