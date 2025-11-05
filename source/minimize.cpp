
#include "../CESolver/CESolver.h"
#include<chrono>

using namespace std;

int main() {


    GasType g = GasType::AIR5; // Set gas type
    ConstraintType constraint = ConstraintType::TV; // Set minimization procedure

    mix gas = common_air::create_air_mix(g); // Create gas mix
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    auto start = NOW;
    CE.compute_equilibrium(1e4, 1.0); // Solve minimization

    cout << gas.T << endl; // Print gas mix temperature.

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();
    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}