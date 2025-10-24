
#include "../equilibrium/equilibrium.h"
#include<chrono>

using namespace std;

int main(int argc, char* argv[]) {

    if (argc < 1) {
        cerr << "Usage: " << argv[0] << " <gas_type>\n";
        return 1;
    }

    // convert command-line args to doubles
    string mix_type = argv[1]; 

    auto start = NOW;

    equilibrium myeq(mix_type);
    myeq.plot_concentrations_for_T_range();

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}