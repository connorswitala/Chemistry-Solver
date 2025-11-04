#include "functions.h"


namespace energy{

    vector<MTask> build_matrix_tasks(Context& cfg) {

        int N = 1 + cfg.HAS_IONS + cfg.NEEDS_T;

        vector<MTask> tasks(N);

        MTask elem = nullptr, charge = nullptr, temp = nullptr;

        switch(cfg.energy_type) {
            case::EnergyType::Helmholtz:
                elem   = &helm::form_elemental;
                charge = &helm::form_charge;

                switch(cfg.constraint_type) {
                    case::ConstraintType::TV:
                    break;

                    case::ConstraintType::UV:
                        temp = &helm::form_U;
                    break;

                    case::ConstraintType::SV:
                        temp = &helm::form_S;
                    break;
                }
            break;

            case::EnergyType::Gibbs:
                elem   = &gibbs::form_elemental;
                charge = &gibbs::form_charge;

                switch(cfg.constraint_type) {
                    case::ConstraintType::TP:
                    break;

                    case::ConstraintType::HP:
                        temp = &gibbs::form_H;
                    break;

                    case::ConstraintType::SP:
                        temp = &gibbs::form_S;
                    break;
                }
            break;
        }

        tasks[0] = elem;
        if (cfg.HAS_IONS) tasks[1] = charge;
        if (cfg.NEEDS_T) tasks[N - 1] = temp;

        return tasks;   
    }

    gasTask build_mu_task(Context& cfg) {
        
        gasTask mu = nullptr;

        switch(cfg.energy_type) {

            case::EnergyType::Helmholtz:
                mu   = &helm::compute_mu;
            break;

            case::EnergyType::Gibbs:
                mu   = &gibbs::compute_mu;
            break;
        }

        return mu;
    }

    NjTask build_Nj_task(Context& cfg) {
        
        NjTask Nj = nullptr;

        switch(cfg.energy_type) {

            case::EnergyType::Helmholtz:
                Nj   = &helm::return_Nj;
            break;

            case::EnergyType::Gibbs:
                Nj   = &gibbs::return_Nj;
            break;
        }

        return Nj;
    }


    namespace gibbs {

        void compute_mu(mix& gas) {
            
        }        

        void form_elemental(double* J, double* F, mix& gas) {

        }

        void form_charge(double* J, double* F, mix& gas) {

        }

        void form_H(double* J, double* F, mix& gas) {

        }

        void form_S(double* J, double* F, mix& gas) {

        }

        void return_Nj(mix& gas){


        }

    }


    namespace helm {

        inline void compute_mu(mix& gas) {
            double Rp = gcon/(10e-5);

            for (int j = 0; j < gas.NS; ++j) {
                gas.mu_RT[j] = gas.mu0_RT[j] + gcon * gas.T * log(gas.N[j] * Rp * gas.T / gas.V);
            }
        }

        // This function forms the entries in J and F for the elemental contrains.
        // It currently uses two if statements to check if it needs ions and T columns.
        // ==== TENTATIVELY COMPLETED ====
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
                    akj_N[j] = gas.a[k * NCOEF + j] * gas.N[j];
                    akj_N_sum += akj_N[j];
                }

                for (int i = 0; i < NE; ++i) {
                    
                    sum = 0.0;
                    for (int j = 0; j < NS; ++j) {
                        sum += akj_N[j] * gas.a[i * NCOEF + j];
                    };

                    J[offset + i] = sum;
                }

                sum = 0.0;
                for (int j = 0; j < NS; ++j) {
                    sum += akj_N[j] * gas.mu_RT[j];
                }

                F[k] = gas.b[k] - akj_N_sum + sum;

                if (gas.HAS_IONS) {

                    sum = 0.0;
                    for (int j = 0; j < NS; ++j) 
                        sum += akj_N[j] * gas.species[j].q;

                    J[offset + NE] = sum;
                }

                if (gas.NEEDS_T) {

                    int idx = NE + gas.HAS_IONS;

                    sum = 0.0;
                    for (int j = 0; j < NS; ++j) 
                        sum += akj_N[j] * (gas.H0_RT[j] - 1.0);

                    J[offset + idx] = sum;
                }

            }
        }

        void form_charge(double* J, double* F, mix& gas) {
            
            vector<double> qj_Nj(gas.NS);
            double qjNj_sum = 0.0, sum;

            for (int j = 0; j < gas.NS; ++j) {
                qj_Nj[j] = gas.species[j].q * gas.N[j]; 
                qjNj_sum += qj_Nj[j];
            }

            for (int i = 0; i < gas.NE; ++i) {

                sum = 0.0;
                for (int j = 0; j < gas.NS; ++j) 
                    sum += qj_Nj[j] * gas.a[i * gas.NE + j];
                
                J[gas.NE * gas.J_SIZE + i] = sum;
            }

            sum = 0.0;
            for (int j = 0; j < gas.NS; ++j) 
                sum += qj_Nj[j] * gas.species[j].q;

            J[gas.NE * gas.J_SIZE + gas.NE] = sum;

            sum = 0.0;
            for (int j = 0; j < gas.NS; ++j) 
                sum += qj_Nj[j] * (gas.mu_RT[j]);

            J[gas.NE * gas.J_SIZE + gas.NE + 1] = sum - qjNj_sum;

            if (gas.NEEDS_T) 
        }

        void form_U(double* J, double* F, mix& gas) {

        }

        void form_S(double* J, double* F, mix& gas) {

        }

        void return_Nj(){

            

        }

    }

}