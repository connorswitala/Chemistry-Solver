
#include "equilibrium/equilibrium.h"
#include<chrono>

using namespace std;

int main() {


    string mix_type;
    mix_type = "air11";    

    auto start = NOW;

    equilibrium myeq(mix_type);
    myeq.plot_concentrations_for_T_range();

    auto end = NOW;
    auto duration = chrono::duration<double>(end - start).count();

    cout << endl << endl << "-- Time elaped: " << fixed << setprecision(4) << duration << " s." << endl << endl;

    return 0;
}