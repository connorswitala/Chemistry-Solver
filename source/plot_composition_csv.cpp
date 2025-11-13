#include "../CESolver/CESolver.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main() {
    GasType g = GasType::AIR5;                     // Set gas type
    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

    mix gas = common_air::create_mix(g);
    CESolver CE(gas, constraint);                   // Construct CESolver

    // Output CSV file
    string filename = "../files/equilibrium_" + gas.name + ".csv";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open " << filename << " for writing.\n";
        return 1;
    }

    // Sweep settings
    double P = 101325.0;
    double T_min = 300.0;       // K
    double T_max = 20000.0;     // K
    int    N     = 2000;

    // Write CSV header
    write << "T [K]";
    for (int j = 0; j < gas.NS; ++j) {
        write << ",Y(" << gas.species[j].name << ")";
    }
    write << "\n";

    // Compute equilibrium across the sweep
    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < N; ++i) {
        double T = T_min + (T_max - T_min) * static_cast<double>(i) / static_cast<double>(N - 1);
        CE.compute_equilibrium(T, P);

        // Write one row per temperature
        write << fixed << setprecision(8) << T;
        for (int j = 0; j < gas.NS; ++j) {
            write << "," << gas.Y[j];
        }
        write << "\n";
    }

    auto end = chrono::high_resolution_clock::now();
    const double duration = chrono::duration<double>(end - start).count();

    write.close();

    cout << "\nTotal time: " << fixed << setprecision(8) << duration
         << " s  -- Time per call = " << duration / N << " s.\n"
         << "-- CSV file saved to " << filename << "\n";

    return 0;
}
