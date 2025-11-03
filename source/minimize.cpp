
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


    GasType g = GasType::AIR5; 
    ConstraintType c = ConstraintType::UV;

    auto start = NOW;

    CESolver solver(g, c);
    solver.compute_equilibrium_T(1e5, 0.01);

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}