
#include "../CESolver/CESolver.h"
#include<chrono>

using namespace std;

int main() {


    GasType g = GasType::AIR11; // Set gas type
    ConstraintType constraint = ConstraintType::TV; // Set minimization procedure

    mix gas = common_air::create_air_mix(g); // Create gas mix
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    int N = 200;
    double U_base = 300, U;
    double V = 1.0;

    auto start = NOW;
    for (int i = 0; i < N; ++i) {
        U =  U_base + i * 50;
        CE.compute_equilibrium(U, V); // Solve minimization
    }
    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();


    cout << endl << "-- Time per call = " << fixed << setprecision(8) << duration / N << " s." << endl << endl;

    return 0;
}