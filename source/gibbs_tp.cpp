
#include "../equilibrium/equilibrium.h"
#include<chrono>

using namespace std;

int main(int argc, char* argv[]) {

    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <T> <p> <gas_type>\n";
        return 1;
    }

    // convert command-line args to doubles
    double T = atof(argv[1]);
    double p   = atof(argv[2]);
    string mix_type = argv[3];   

    auto start = NOW;

    equilibrium myeq(mix_type);
    
    myeq.compute_equilibrium_TP(T, p);

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    myeq.display_gas_properties();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(6) << duration << " s." << endl << endl;

    return 0;
}