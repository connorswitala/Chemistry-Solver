
#include "solverLib/equilibrium.h"
#include<chrono>

using namespace std;

int main() {

    string mix_type;
    mix_type = "air11";    

    auto start = NOW;

    equilibrium myeq(mix_type, 0.1, 6e5);
    myeq.compute_equilibrium();

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();
    cout << "-- Time elaped: " << duration << " s." << endl << endl;

    myeq.diplay_gas_properties();

    return 0;
}