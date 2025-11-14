#ifndef MIXES_H
#define MIXES_H


#include "../includes/thermoObjs.h"
#include "../includes/readNASA.h"

using namespace std;


// Prints NASA thermodynamic data for gas mix
inline void print_NASA_mix(const mix& gas);

// Prints NASA thermodynamic data for a vector of strings holding species names
void print_NASA(vector<string> species_names);

namespace common {

    // Creates a mixture using the GasType enum parameter for creating a common air mixture.
    mix air_mixture(GasType gastype);

    // Basic funtions that create specific air mixtures. Called in above function;
    mix make_air13();
    mix make_air11_Ar();
    mix make_air11();
    mix make_air7();
    mix make_air5();
    mix make_perf();

}

// Used for user-specification of mixture
mix create_mixture(vector<string>& speciesNames, 
                   vector<string>& elementNames, 
                   vector<double>& initial_Y); 

#endif