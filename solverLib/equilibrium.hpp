#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "../includes/common.h"
#include "../includes/thermoObjs.h"


using namespace std;



class equilibrium {
        
    public:

    const int n_sp;
    const int n_el;
    const int n_coeff;

    mix gas;
    int T_flag;


    private:

    Vector Ts;

    void create_air_5();
    void create_air_11();
    int findTRange();
    array<double, 7> temp_base(double T);

    void NASA_fits();
    void SfromT();
    void MUfromT();
   

};


#endif

