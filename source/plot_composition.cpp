#include "../CESolver/CESolver.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main() {
    GasType g = GasType::AIR11;                     // Set gas type
    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

    mix gas = common_air::create_air_mix(g);        // Create gas mix
    CESolver CE(gas, constraint);                   // Construct CESolver

    // Output Tecplot file
    string filename = "../files/equilibrium_" + gas.name + ".dat";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open " << filename << " for writing.\n";
        return 1;
    }

    // Sweep settings
    double P = 101325.0;   // Pa
    double T_min = 300.0;  // K
    double T_max = 20000.0;// K
    int    N     = 2000;

    // Buffers for BLOCK output
    vector<double> Tvals(N, 0.0);
    vector<double> Pvals(N, P);
    vector<vector<double>> Y(gas.NS, vector<double>(N, 0.0));

    // Compute equilibrium across the sweep
    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < N; ++i) {
        double T = T_min + (T_max - T_min) * static_cast<double>(i) / static_cast<double>(N - 1);
        CE.compute_equilibrium(T, P);

        Tvals[i] = T;
        for (int j = 0; j < gas.NS; ++j) {
            Y[j][i] = gas.Y[j];
        }
    }

    auto end = chrono::high_resolution_clock::now();
    const double duration = chrono::duration<double>(end - start).count();

    // ===== Tecplot ASCII Header (BLOCK format) =====
    // Title
    write << "TITLE = \"Equilibrium TP Sweep: " << gas.name << "\"\n";

    // Variables: T, P, and species mass fractions
    write << "VARIABLES = \"T [K]\", \"P [Pa]\"";
    for (int j = 0; j < gas.NS; ++j) {
        write << ", \"Y(" << gas.species[j].name << ")\"";
    }
    write << "\n";

    // One 1D zone with N points
    write << "ZONE T=\"TP_Sweep\", I=" << N << ", F=BLOCK\n";

    // Helper lambda to dump a vector with consistent formatting
    auto dump_vec = [&](const vector<double>& v) {
        write << scientific << setprecision(8);
        int count = 0;
        for (double val : v) {
            write << val << " ";
            // keep lines a bit shorter for readability
            if (++count % 5 == 0) write << "\n";
        }
        if (count % 5 != 0) write << "\n";
    };

    // BLOCK order: all T, then all P, then each Y_j
    dump_vec(Tvals);
    dump_vec(Pvals);
    for (int j = 0; j < gas.NS; ++j) dump_vec(Y[j]);

    write.close();

    cout << "\nTotal time: " << fixed << setprecision(8) << duration
         << " s  -- Time per call = " << duration / N << " s.\n"
         << "-- Tecplot file saved to " << filename << "\n";

    return 0;
}
