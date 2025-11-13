#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {
    
    /**
     * Required inputs: species names, and initial mass fractions of each species.
    */
    // vector<string> species = {"N2", "O2", "NO", "N", "O", "N2+", "O2+", "N+", "O+", "NO+", "e-"};        
    // vector<string> elements = {"N", "O"};
    // vector<double> initial = {0.7572, 0.2428}; 

    vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
    vector<string> elements = {"C", "O", "N"};
    vector<double> initial = {0.264, 0.7186, 0.0174}; 

     mix gas = create_mix::mixture(species, elements, initial);

    // GasType g = GasType::AIR5;
    // mix gas = common_air::create_mix(g);

    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double T = 500;
    double P = 100000;

    CE.compute_equilibrium(T,P); // Solve minimization 
    CE.print_properties();

    return 0;
}
