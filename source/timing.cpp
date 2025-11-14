#include "../libraries/CESolver/CESolver.h"

using namespace std;

int main() {

    // Random Number Generators for testing time.
    std::random_device rd;                      // non-deterministic seed
    std::mt19937 gen(rd());                     // Mersenne Twister RNG
    std::uniform_int_distribution<> distP(20000.0, 100000.0);
    std::uniform_int_distribution<> distT(300.0, 20000.0);
    std::uniform_int_distribution<> distV(1.0, 100.0);
    std::uniform_int_distribution<> distU(5e5, 5e7);


    GasType g = GasType::AIR11;
    mix gas = common::air_mixture(g);
    ConstraintType constraint = ConstraintType::UV; // Set minimization procedure

    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double T, T_base = 300.0;
    double U, U_base = 5e5;
    double P = 101325.0, V = 1.0 / 1.225;

    int N = 500000;
    auto start = NOW;
    for (int i = 0; i < N; ++i) {
        U = distU(gen);
        V = distV(gen);
        CE.compute_equilibrium(U, V); // Solve minimization 
    }
    
    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    // CE.print_properties();

    cout << endl << "Total time: " << setprecision(8) << duration << "-- Time per call = " << fixed << setprecision(8) << duration / N << " s." << endl << endl;

    return 0;
}
