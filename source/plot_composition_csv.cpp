#include "../libraries/CESolver/CESolver.h"

using namespace std;

int main() {
    GasType g = GasType::AIR11;                     // Set gas type
    ConstraintType constraint = ConstraintType::UV; // Set minimization procedure

    mix gas = common::air_mixture(g);
    CESolver CE(gas, constraint);                   // Construct CESolver

    // Output CSV file
    string filename = "../misc/files/equilibrium_" + gas.name + ".csv";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open " << filename << " for writing.\n";
        return 1;
    }

    // Sweep settings
    double U_min = 3e5;       // K
    double V = 1.0 / 1.225;
    int i = 0;

    // Write CSV header
    write << "T [K]";
    for (int j = 0; j < gas.NS; ++j) {
        write << ",Y(" << gas.species[j].name << ")";
    }
    write << "\n";

    // Compute equilibrium across the sweep
    auto start = chrono::high_resolution_clock::now();

    while (gas.T < 19500.0) {
        double  U = U_min + 175.0 * i;
        CE.compute_equilibrium(U, V);

        // Write one row per temperature
        write << fixed << setprecision(8) << gas.T;
        for (int j = 0; j < gas.NS; ++j) {
            write << "," << gas.Y[j];
        }
        write << "\n";
        i++;
    }

    auto end = chrono::high_resolution_clock::now();
    const double duration = chrono::duration<double>(end - start).count();

    write.close();

    cout << "\nTotal time: " << fixed << setprecision(8) << duration
         << " s  -- Time per call = " << duration / i << " s.\n"
         << "-- CSV file saved to " << filename << "\n";
    cout << "Iterations: " << i << endl;

    return 0;
}
