#ifndef COMMONMIXES_H
#define COMMONMIXES_H


#include "../includes/thermoObjs.h"
#include "../includes/readNASA.h"

using namespace std;


inline void print_species_table(const mix& gas);

namespace common_air {

    mix create_mix(GasType gastype);
    mix make_air13();
    mix make_air11_Ar();
    mix make_air11();
    mix make_air7();
    mix make_air5();
    mix make_perf();

}

namespace create_mix {
    mix mixture(vector<string>& speciesNames, vector<string>& elementNames, vector<double>& initial_Y); 
}


#endif