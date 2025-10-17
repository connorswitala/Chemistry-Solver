
#include "solverLib/equilibrium.h"
#include<chrono>

using namespace std;

int main() {

    string mix_type;
    mix_type = "air11";    

    auto start = NOW;

    equilibrium myeq(mix_type);
    myeq.compute_equilibrium(0.1, 1e7);

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();
    cout << endl << "-- Time elaped: " << duration << " s." << endl << endl;

    return 0;
}