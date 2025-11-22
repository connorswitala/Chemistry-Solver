#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {

    // ======: Create gas mix using enum class for common mixes :======

    mix gas = common_mixture(GasType::AIR11);                        


    // =======: Create user-defined mix (Mars 8-species) :======

    // vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
    // vector<string> elements = {"C", "O", "N"};
    // vector<double> initial = {0.264, 0.7186, 0.0174}; 

    // mix gas = create_mixture(species, elements, initial);

    CESolver CE(gas, ConstraintType::TP);                       // Construct CESolver class for minimization, using TP constraint.

    CE.compute_equilibrium(5000.0, 101325.0);  // Solve for equilibrium at T = 5000.0 K, and P = 101325.0 Pa

    print_properties(gas);  // Print mixture properties
    print_NASA_mix(gas);    // Print NASA thermodynamic table data for mixture (used for debugging)

    cout << gas.R << endl;      // Print just mixture specific gas constant
    cout << gas.Y[0] << endl;   // Print just species [0] mass fraction.
    
    return 0;
}
