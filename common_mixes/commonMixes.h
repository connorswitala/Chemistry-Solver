#ifndef COMMONMIXES_H
#define COMMONMIXES_H

#include "../includes/common.h"
#include "../includes/thermoObjs.h"

using namespace std;

namespace common_air {

    const vector<SpeciesInfo>& air11();
    const vector<SpeciesInfo>& air5(); 

    mix make_air11();
    mix make_air5();


}


#endif