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

    CE.print_properties();

    cout << endl << "Old time: " << setprecision(8) << duration << endl;


    constraint = ConstraintType::UV; // Set minimization procedure
    double U = 1e7;
    double V = 1 / 1.225;

    CESolver CE1(gas, constraint);   // Construct CESolver class for minimization.

    start = NOW;

    CE1.compute_equilibrium(U,V); // Solve minimization 
    
    end = NOW;
    duration = chrono::duration<double>(end - start).count();

    CE1.print_properties();

    cout << endl << "Old time: " << setprecision(8) << duration << endl;
    return 0;
}
