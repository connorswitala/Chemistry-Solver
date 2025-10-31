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

    muTask build_mu_task(Context& cfg) {
        
        muTask mu = nullptr;

        switch(cfg.energy_type) {

            case::EnergyType::Helmholtz:
                mu   = &helm::compute_mu;
            break;

            case::EnergyType::Gibbs:
                mu   = &gibbs::compute_mu;
            break;
        }
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

    }


    namespace helm {

        void compute_mu(mix& gas) {
            double Rp = gcon/(10e-5);

            for (int j = 0; j < gas.NS; ++j) {
                gas.mu_RT[j] = gas.mu0_RT[j] + gcon * gas.T * log(gas.N[j] * Rp * gas.T / gas.V);
            }
        }

        void form_elemental(double* J, double* F, mix& gas) {

            int NS = gas.NS;
            int NE = gas.NE; 
            int J_SIZE = gas.J_SIZE;
            
            double sum;
            double akj_N_sum;

            vector<double> akj_N(NS);

            for (int k = 0; k < NE; ++k) {

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

                    J[k * J_SIZE + i] = sum;
                }

                sum = 0.0;
                for (int j = 0; j < NS; ++j) {
                    sum += akj_N[j] * gas.mu_RT[j];
                }

                F[k] = gas.b[k] - akj_N_sum + sum;
            }
        }

        void form_charge(double* J, double* F, mix& gas) {

        }

        void form_U(double* J, double* F, mix& gas) {

        }

        void form_S(double* J, double* F, mix& gas) {

        }

    }

}