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


    GasType g = GasType::AIR11; // Set gas type
    ConstraintType constraint = ConstraintType::TP; // Set minimization procedure

    mix gas = common_air::create_air_mix(g); // Create gas mix
    CESolver CE(gas, constraint);   // Construct CESolver class for minimization.

    double U = 1e4;
    double P = 101325.0;
    double T = 300.0;
    double V = 1.0;

    int N = 500000;
    auto start = NOW;
    for (int i = 0; i < N; ++i) {
        T = distT(gen);
        P = distP(gen);
        CE.compute_equilibrium(T, P); // Solve minimization 
    }
    
    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    CE.print_properties();


    cout << endl << "Total time: " << setprecision(8) << duration << "-- Time per call = " << fixed << setprecision(8) << duration / N << " s." << endl << endl;

    return 0;
}
