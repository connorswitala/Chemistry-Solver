#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {

    // ====== Create gas mix using enum class for common mixes

    mix gas = common_mixture(GasType::AIR5);                        


    // ======= Create user-defined mix ======

    // vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
    // vector<string> elements = {"C", "O", "N"};
    // vector<double> initial = {0.264, 0.7186, 0.0174}; 

    // mix gas = create_mixture(species, elements, initial);

    ConstraintType constraint = ConstraintType::TP;     // Set minimization procedure (temperature, pressure held constant in this case).
    CESolver CE(gas, constraint);                       // Construct CESolver class for minimization.


    string filename = "../user-files/plot.dat";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open " << filename << " for writing.\n";
        exit;
    }

    // Buffers for BLOCK output
    int N = 5000;
    double T;

        // ===== Tecplot ASCII Header (BLOCK format) =====
    // Title
    write << "TITLE = \"Equilibrium TP Sweep: " << gas.name << "\"\n";

    // Variables: T, P, and species mass fractions
    write << "VARIABLES = \"T [K]\"";
    for (int j = 0; j < gas.NS; ++j) {
        write << ", \"X[" << gas.species[j].name << "]\"";
    }
    write << "\n";

    // One 1D zone with N points
    write << "ZONE T=\"Sweep\", I=" << N << ", F=POINT\n";

    for (int i = 0; i < N; ++i) {
        T = 300.0 + (20000.0 - 300.0) / (N - 1) * double(i);
        CE.compute_equilibrium(T, 101325.0);

        write << T;
        for (int j = 0; j < gas.NS; ++j)
            write << ", " << gas.X[j];
        write << endl;
    }

    return 0;
}
