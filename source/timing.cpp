#include "../CESolver/CESolver.h"
#include<chrono>
#include<random> 

using namespace std;

int main() {

    // Random Number Generators for testing time.
    std::random_device rd;                      // non-deterministic seed
    std::mt19937 gen(rd());                     // Mersenne Twister RNG
    std::uniform_int_distribution<> distP(20000.0, 100000.0);
    std::uniform_int_distribution<> distT(300.0, 20000.0);
    std::uniform_int_distribution<> distV(1.0, 2.0);


    GasType g = GasType::AIR11; // Set gas type
    ConstraintType constraint = ConstraintType::TV; // Set minimization procedure

    mix gas = common_air::create_mix(g);
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double T, T_base = 300.0, P;
    double U, V = 1.0;

    int N = 500000;
    auto start = NOW;
    for (int i = 0; i < N; ++i) {
        T = T_base + i * (20000.0 - 300.0) / N;
        CE.compute_equilibrium(T, V); // Solve minimization 
    }
    
    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    CE.print_properties();

    cout << endl << "OLD version total time: " << setprecision(8) << duration << "-- Time per call = " << fixed << setprecision(8) << duration / N << " s." << endl << endl;

    // ====================== NEW VERSION ========================

    auto start1 = NOW;
    for (int i = 0; i < N; ++i) {
        T = T_base + i * (20000.0 - 300.0) / N;
        CE.CFD_equilibrium(T, V); // Solve minimization 
    }
    
    auto end1 = NOW;
    auto duration1 = chrono::duration<double>(end1 - start1).count();

    CE.print_properties();

    cout << endl << "NEW version total time: " << setprecision(8) << duration1 << "-- Time per call = " << fixed << setprecision(8) << duration1 / N << " s." << endl << endl;

    return 0;
}
