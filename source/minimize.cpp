#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {


    GasType g = GasType::AIR7; // Set gas type
    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

    mix gas = common::air_mixture(g);
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double T = 6000;
    double P = 101325.0;

    auto start = NOW;

    CE.compute_equilibrium(T,P); // Solve minimization 
    
    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    print_properties(gas);

    cout << endl << "Old time: " << setprecision(8) << duration << endl;

    return 0;
}
