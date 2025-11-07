#include "../CESolver/CESolver.h"
#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

int main() {
    GasType g = GasType::AIR11;                  // Set gas type
    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

    mix gas = common_air::create_air_mix(g);    // Create gas mix
    CESolver CE(gas, constraint);               // Construct CESolver

    string filename = "../files/equilibrium_" + gas.name + ".csv";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open equilibrium.csv for writing.\n";
        return 1;
    }

    // Header: T, P, then species names from gas.species[j].name
    write << "T,P";
    for (int j = 0; j < gas.NS; ++j) {
        write << ',' << gas.species[j].name;    
    }
    write << '\n';

    double T, P = 101325.0, T_base = 300.0;
    int N = 2000;

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        T = T_base + i * (20000.0 - 300.0)/N;

        CE.compute_equilibrium(T, P);           

        write << fixed << setprecision(6) << T << ',' << P;
        for (int j = 0; j < gas.NS; ++j) {
            write << ',' << setprecision(8) << gas.Y[j];
        }
        write << '\n';
    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double>(end - start).count();

    write.close();

    cout << "\nTotal time: " << setprecision(8) << duration
         << " -- Time per call = " << fixed << setprecision(8) << duration / N << " s.\n\n";
         cout << "-- File saved to " << filename << endl;

    return 0;
}
