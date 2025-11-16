#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {

    // ====== Create gas mix using enum class for common mixes

    // GasType g = GasType::MARS8;                         // Set gas type
    // mix gas = common_mixture(g);                        // Create gas mixture from common species


    // ======= Create user-defined mix ======

    vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
    vector<string> elements = {"C", "O", "N"};
    vector<double> initial = {0.264, 0.7186, 0.0174}; 

    mix gas = create_mixture(species, elements, initial);


    ConstraintType constraint = ConstraintType::TP;     // Set minimization procedure (temperature, pressure held constant in this case).
    CESolver CE(gas, constraint);                       // Construct CESolver class for minimization.

    // Temperature and pressure to solve for
    double T = 210;
    double P = 610;

    CE.compute_equilibrium(T,P);    // Solve
    print_properties(gas);          // Print mixture properties
    print_NASA_mix(gas);            // Print NASA data takesn from thermodynamic tables

    return 0;
}
