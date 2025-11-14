#include "../libraries/CESolver/CESolver.h"

using namespace std;

int main() {

   
    // GasType g = GasType::AIR13;
    // mix gas = common::air_mixture(g);

    vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
    vector<string> elements = {"C", "O", "N"};
    vector<double> initial = {0.264, 0.7186, 0.0174}; 

     mix gas = create_mixture(species, elements, initial);


    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double T = 210;
    double P = ATM;
    double U = 522131;
    double V = 62;

    CE.compute_equilibrium(T,P); // Solve minimization 
    CE.print_properties();

    return 0;
}
