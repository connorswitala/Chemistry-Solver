#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "../includes/thermoObjs.h"


using MTask = void(*)(double*, double*, mix&);
using muTask = void(*)(mix&);


namespace energy {

    vector<MTask> build_matrix_tasks(Context& cfg); 
    muTask build_mu_task(Context& cfg);


    namespace gibbs {

        void compute_mu(mix& gas);
        void form_charge(double* J, double* F, mix& gas);
        void form_temp(double* J, double* F, mix& gas);
        void form_elemental(double* J, double* F, mix& gas);

    }


    namespace helm {

        void compute_mu(mix& gas);
        void form_charge(double* J, double* F, mix& gas);
        void form_temp(double* J, double* F, mix& gas);
        void form_elemental(double* J, double* F, mix& gas);

    }
}






#endif