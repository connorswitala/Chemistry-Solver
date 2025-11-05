
#include "../CESolver/CESolver.h"
#include<chrono>

using namespace std;

int main() {


    GasType g           = GasType::AIR5; 
    ConstraintType c    = ConstraintType::TV;

    mix gas = common_air::create_air_mix(g);

    CESolver CE(gas, c);

    auto start = NOW;
    CE.compute_equilibrium(1e4, 1.0);

    cout << gas.T << endl;

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}