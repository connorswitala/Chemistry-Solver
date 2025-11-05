
#include "../CESolver/CESolver.h"
#include<chrono>

using namespace std;

// int argc, char* argv[]

// if (argc < 4) {
//         cerr << "Usage: " << argv[0] << " <rho> <e> <gas_type>\n";
//         return 1;
//     }

//     // convert command-line args to doubles
//     double rho = atof(argv[1]);
//     double e   = atof(argv[2]);
//     string mix_type = argv[3];   


int main() {


    GasType g           = GasType::AIR5; 
    ConstraintType c    = ConstraintType::TP;

    CESolver solver(g, c);

    auto start = NOW;
    solver.compute_equilibrium(1e4, 1.0);

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}