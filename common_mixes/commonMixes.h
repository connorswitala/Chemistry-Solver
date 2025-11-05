#ifndef COMMONMIXES_H
#define COMMONMIXES_H


#include "../includes/thermoObjs.h"

using namespace std;

namespace common_air {

    mix create_air_mix(GasType gastype); 

    mix make_air13();
    mix make_air11_Ar();
    mix make_air11();
    mix make_air5();
    mix make_perf();

}


#endif