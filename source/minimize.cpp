
// #include "../CESolver/CESolver.h"
// #include<chrono>
// #include<random> 

// using namespace std;

// int main() {

//     // // Random Number Generators for testing time.
//     // std::random_device rd;                      // non-deterministic seed
//     // std::mt19937 gen(rd());                     // Mersenne Twister RNG
//     // std::uniform_int_distribution<> distP(20000.0, 100000.0);
//     // std::uniform_int_distribution<> distT(300.0, 20000.0);


//     GasType g = GasType::AIR7; // Set gas type
//     ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

//     mix gas = common_air::create_air_mix(g); // Create gas mix
//     CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

//     double T, P = 101325.0, T_base = 300.0;
//     int N = 2000;

//     auto start = NOW;
//     for (int i = 0; i < N; ++i) {
//         T = T_base + i * 8.0;
//         CE.compute_equilibrium(T, P); // Solve minimization

//         for (int j = 0; j < NS; ++j) 
//         write << gas.X[j];

//     }
//     auto end = NOW;
//     auto duration = chrono::duration<double>(end - start).count();


//     cout << endl << "Total time: " << setprecision(8) << duration << "-- Time per call = " << fixed << setprecision(8) << duration / N << " s." << endl << endl;

//     return 0;
// }


#include "../CESolver/CESolver.h"
#include <chrono>
#include <random>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

int main() {
    GasType g = GasType::AIR7;                  // Set gas type
    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

    mix gas = common_air::create_air_mix(g);    // Create gas mix
    CESolver CE(gas, constraint);               // Construct CESolver

    // === Open CSV and write header ===
    ofstream write("equilibrium.csv");
    if (!write) {
        cerr << "Failed to open equilibrium.csv for writing.\n";
        return 1;
    }

    // Header: T, P, then species names from gas.species[j].name
    write << "T,P";
    for (int j = 0; j < gas.NS; ++j) {
        write << ',' << gas.species[j].name;    // assumes .name is a simple identifier
    }
    write << '\n';

    double T, P = 101325.0, T_base = 300.0;
    int N = 2000;

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        T = T_base + i * 8.0;

        CE.compute_equilibrium(T, P);           // Solve minimization (updates gas via CE)

        // One CSV row: T, P, then X[j]
        write << fixed << setprecision(6) << T << ',' << P;
        for (int j = 0; j < gas.NS; ++j) {
            write << ',' << scientific << setprecision(8) << gas.X[j];
        }
        write << '\n';
    }
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double>(end - start).count();

    write.close(); // optional (also closes in destructor)

    cout << "\nTotal time: " << setprecision(8) << duration
         << " -- Time per call = " << fixed << setprecision(8) << duration / N << " s.\n\n";

    return 0;
}
