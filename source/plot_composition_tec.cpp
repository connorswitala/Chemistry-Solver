#include "../libraries/CESolver/CESolver.h"

using namespace std;

int main() {
    GasType g = GasType::AIR11;                     // Set gas type
    ConstraintType constraint = ConstraintType::UV; // Set minimization procedure

    mix gas = common_mixture(g);
    CESolver CE(gas, constraint);                   // Construct CESolver

    // Output Tecplot file
    string filename = "../misc/files/equilibrium_" + gas.name + ".dat";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open " << filename << " for writing.\n";
        return 1;
    }

    // Sweep settings
    double V = 1.0/1.225;   // Pa
    double e_min = 5e5;  // K

    // Buffers for BLOCK output
    vector<double> Tvals;
    vector<vector<double>> Y(gas.NS, vector<double>(0.0));
    vector<double> dpdr;
    vector<double> dtdr;
    vector<double> cv;

    // Compute equilibrium across the sweep
    auto start = chrono::high_resolution_clock::now();

    int i = 0;
    while (gas.T < 19500.0) {
        double e = e_min + i * 500.0;
        CE.compute_equilibrium(e, V);

        Tvals.push_back(gas.T);
        for (int j = 0; j < gas.NS; ++j) {
            Y[j].push_back(gas.Y[j]);
        }
        dpdr.push_back(gas.dpdr);
        dtdr.push_back(gas.dtdr);
        cv.push_back(gas.cv);
        i++;
    }

    auto end = chrono::high_resolution_clock::now();
    const double duration = chrono::duration<double>(end - start).count();

    // ===== Tecplot ASCII Header (BLOCK format) =====
    // Title
    write << "TITLE = \"Equilibrium TP Sweep: " << gas.name << "\"\n";

    // Variables: T, P, and species mass fractions
    write << "VARIABLES = \"T [K]\"";
    for (int j = 0; j < gas.NS; ++j) {
        write << ", \"Y(" << gas.species[j].name << ")\"";
    }
    write << ", \"dpdr\", \"dtdr\", \"cv\"";
    write << "\n";

    // One 1D zone with N points
    write << "ZONE T=\"TP_Sweep\", I=" << Y[0].size() << ", F=BLOCK\n";

    // Helper lambda to dump a vector with consistent formatting
    auto dump_vec = [&](const vector<double>& v) {
        write << scientific << setprecision(8);
        int count = 0;
        for (double val : v) {
            write << val << " ";
            // keep lines a bit shorter for readability
            if (++count % 5 == 0) write << "\n";
        }
        if (count % 5 != 0) write << "\n";
    };

    // BLOCK order: all T, then all P, then each Y_j
    dump_vec(Tvals);
    for (int j = 0; j < gas.NS; ++j) dump_vec(Y[j]);
    dump_vec(dpdr);
    dump_vec(dtdr);
    dump_vec(cv);

    write.close();

    cout << "\nTotal time: " << fixed << setprecision(8) << duration
         << " s  -- Time per call = " << scientific << duration / i << " s.\n"
         << "-- Tecplot file saved to " << filename << "\n";

    return 0;
}
