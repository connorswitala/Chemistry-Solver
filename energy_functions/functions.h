#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "../includes/thermoObjs.h"

using namespace std;

using MTask = void(*)(double*, double*, mix&);
using gasTask = void(*)(mix&);
using NjTask = void(*)(mix&, double*);


namespace energy {

    vector<MTask> build_matrix_tasks(Context& cfg); 
    gasTask build_mu_task(Context& cfg);
    NjTask build_Nj_task(Context& cfg);


    namespace gibbs {

        void compute_mu(mix& gas);

        void form_charge(double* J, double* F, mix& gas);
        void form_H(double* J, double* F, mix& gas);
        void form_S(double* J, double* F, mix& gas);
        void form_elemental(double* J, double* F, mix& gas);

        void return_Nj(mix& gas, double* DELTA);

    }


    namespace helm {

        inline void compute_mu(mix& gas);

        inline void form_charge(double* J, double* F, mix& gas);
        inline void form_U(double* J, double* F, mix& gas);
        inline void form_S(double* J, double* F, mix& gas);
        inline void form_elemental(double* J, double* F, mix& gas);

        inline void return_Nj(mix& gas, double* DELTA);

    }
}






#endif