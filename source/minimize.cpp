#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {


    GasType g = GasType::AIR11; // Set gas type
    ConstraintType constraint = ConstraintType::TV; // Set minimization procedure

    mix gas = common_air::create_air_mix(g); // Create gas mix
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double T = 5000.0;
    double V = 1.0;

    auto start = NOW;

    CE.compute_equilibrium(T, V); // Solve minimization 
    
    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    CE.print_properties();

    cout << endl << "Old time: " << setprecision(8) << duration << endl;

    start = NOW;

    CE.CFD_equilibrium(T, V); // Solve minimization 
    
    end = NOW;
    duration = chrono::duration<double>(end - start).count();

    CE.print_properties();

    cout << endl << "New time: " << setprecision(8) << duration << endl;
    return 0;
}
