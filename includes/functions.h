#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "thermoObjs.h"

using namespace std;




// namespace gibbs {

//     void compute_mu(mix& gas) {
        
//     }        

//     void form_elemental(double* J, double* F, mix& gas) {

//     }

//     void form_charge(double* J, double* F, mix& gas) {

//     }

//     void form_H(double* J, double* F, mix& gas) {

//     }

//     void form_S(double* J, double* F, mix& gas) {

//     }

// }


namespace helm {

    inline void compute_mu(mix& gas) {
        double Rp = gcon*(1.0e-5);

        for (int j = 0; j < gas.NS; ++j) {
            gas.mu_RT[j] = gas.mu0_RT[j] + log(gas.N[j] * Rp * gas.T / gas.V);
        }
    }

    // This function forms the entries in J and F for the elemental contrainys.
    // It currently uses two if statements to check if it needs ions and T columns.
    // ==== COMPLETED ====
    inline void form_elemental(double* J, double* F, mix& gas) {

        int NS = gas.NS;
        int NE = gas.NE; 
        int J_SIZE = gas.J_SIZE;
        
        double sum;
        double akj_N_sum;

        int offset;

        vector<double> akj_N(NS);
        for (int k = 0; k < NE; ++k) {

            offset = k * J_SIZE;

            akj_N_sum = 0.0;
            for (int j = 0; j < NS; ++j) {
                akj_N[j] = gas.a[k * NS + j] * gas.N[j];
                akj_N_sum += akj_N[j];
            }

            for (int i = 0; i < NE; ++i) {
                
                sum = 0.0;
                for (int j = 0; j < NS; ++j) {
                    sum += akj_N[j] * gas.a[i * NS + j];
                };
                J[offset + i] = sum;
            }
            

            sum = 0.0;
            for (int j = 0; j < NS; ++j) {
                sum += akj_N[j] * gas.mu_RT[j];
            }

            F[k] = gas.b[k] - akj_N_sum + sum;
        }
    }

    // This function forms the entries in J and F for the charge contraint.
    // It currently uses one if statement to check if it needs a T column.
    // ==== COMPLETED ====
    inline void form_charge(double* J, double* F, mix& gas) {
        
        int NS = gas.NS;
        int NE = gas.NE;
        int J_SIZE = gas.J_SIZE;

        vector<double> qj_Nj(NS);
        double qjNj_sum = 0.0, sum;

        // Sum terms used a lot.
        for (int j = 0; j < NS; ++j) {
            qj_Nj[j] = gas.species[j].q * max(1e-12, gas.N[j]); 
            qjNj_sum += qj_Nj[j];
        }

        // Charge column in elemental constraint rows
        for (int k = 0; k < NE; ++k) {
            sum = 0.0;
            for (int j = 0; j < NS; ++j) 
                sum += qj_Nj[j] * gas.a[k * NS + j];

            J[k * J_SIZE + NE] = sum;    
        }

        // Elemental contraint columns in charge row.
        for (int i = 0; i < NE; ++i) {

            sum = 0.0;
            for (int j = 0; j < NS; ++j) 
                sum += qj_Nj[j] * gas.a[i * NS + j];
            
            J[NE * J_SIZE + i] = sum;
        }

        // Charge constraint column in charge row.
        sum = 0.0;
        for (int j = 0; j < NS; ++j) 
            sum += qj_Nj[j] * gas.species[j].q;

        J[NE * J_SIZE + NE] = sum;

        // RHS
        sum = 0.0;
        for (int j = 0; j < NS; ++j) 
            sum += qj_Nj[j] * gas.mu_RT[j];

        F[NE] = sum - qjNj_sum;

        // Add charge column in ln(T) row is present.
        if (gas.NEEDS_T) {

            sum = 0.0;
            for (int j = 0; j < NS; ++j)
                sum += qj_Nj[j] * (gas.H0_RT[j] - 1.0);
                
            J[(J_SIZE - 1) * J_SIZE + NE] = sum;
        }        
    }

    // This function forms the entries in J and F for the U contraint.
    // It currently uses one if statement to check if it needs a chage column.
    // ==== COMPLETED ====
    inline void form_U(double* J, double* F, mix& gas) {

        int NE = gas.NE;
        int NS = gas.NS;
        int J_SIZE = gas.J_SIZE;

        int idx = NE + gas.HAS_IONS;

        double sum, sum1;
        vector<double> Uj_Nj(NS);

        // Importans sum terms
        for (int j = 0; j < NS; ++j) {
            Uj_Nj[j] = gas.N[j] * (gas.H0_RT[j] - 1.0);
        }

        // ln(T) column in elemental rows
        for (int k = 0; k < NE; ++k) {
            sum = 0.0;
            for (int j = 0; j < NS; ++j) 
                sum += gas.a[k * NS + j] * Uj_Nj[j];

            J[k * J_SIZE + idx] = sum;
        }

        // ln(T) term in charge constraint row
        if (gas.HAS_IONS) {

            sum = 0.0;
            for (int j = 0; j < NS; ++j) {
                sum += gas.species[j].q * Uj_Nj[j];
            }  

            J[NE * J_SIZE + NE + 1] = sum;
        }


        // Elemental constraint columns in ln(T) row
        for (int i = 0; i < NE; ++i) {

            sum = 0.0;
            for (int j = 0; j < NS; ++j) 
                sum += gas.a[i * NS + j] * Uj_Nj[j];

            J[idx * J_SIZE + i] = sum;
        }

        // ln(T) column in ln(T) row
        sum = 0.0;
        sum1 = 0.0;
        for (int j = 0; j < NS; ++j) {
            sum += Uj_Nj[j] * (gas.H0_RT[j] - 1.0);
            sum1 += gas.N[j] * (gas.CP0_R[j] - 1.0);
        }

        J[J_SIZE * J_SIZE - 1] = sum + sum1;

        sum = 0.0;
        for (int j = 0; j < NS; ++j) 
            sum += Uj_Nj[j] * gas.mu_RT[j];

        F[J_SIZE - 1] = (gas.uo - gas.up) / (gcon * gas.T) + sum;
    }


    inline void form_S(double* J, double* F, mix& gas) {

    }

}



#endif