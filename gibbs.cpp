
#include "equilibrium/equilibrium.h"
#include<chrono>

using namespace std;

int main(int argc, char* argv[]) {

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <rho> <e>\n";
        return 1;
    }

    // convert command-line args to doubles
    double rho = atof(argv[1]);
    double e   = atof(argv[2]);

    string mix_type;
    mix_type = "air11";    

    auto start = NOW;

    equilibrium myeq(mix_type);
    myeq.compute_equilibrium(rho, e);

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    myeq.display_gas_properties();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}