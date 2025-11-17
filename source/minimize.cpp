#include "../libraries/CESolver/CESolver.h"

using namespace std;

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

enum class RunMode {
    Single, Sweep
};

// Configuration/options struct
struct Config {
    RunMode mode = RunMode::Single;
    string modestring = "single";
    ConstraintType constraint = ConstraintType::TP;
    string constraintstring = "TP";
    vector<string> species;
    vector<string> elements;
    vector<double> initial;
    mix gas;
    bool creation;
    double val1, val2;
    double maxval, minval;
};

// Simple helper to trim whitespace from both ends
string trimIO(const string& s) {
    size_t first = s.find_first_not_of(" \t\n\r");
    if (first == string::npos) return "";
    size_t last = s.find_last_not_of(" \t\n\r");
    return s.substr(first, last - first + 1);
}

// Split a comma-separated list into vector<string>
vector<string> splitCommaList(const string& s) {
    vector<string> result;
    string item;
    stringstream ss(s);

    while (getline(ss, item, ',')) {
        item = trimIO(item);
        if (!item.empty()) result.push_back(item);
    }
    return result;
}

// Split a comma-separated list into a vector<double>
vector<double> splitCommaDoubles(const string& s) {
    vector<double> result;
    string item;
    stringstream ss(s);

    while (getline(ss, item, ',')) {
        item = trim(item);
        if (item.empty()) continue;
        try {
            double val = stod(item);
            result.push_back(val);
        } catch (const exception& e) {
            cerr << "Could not parse '" << item 
                      << "' as a number in --Y list\n";
        }
    }
    return result;
}

// Handle a single line like "--mode=TPsweep" or "--species=N2,O2"
void handle_option_line(const string& line, Config& cfg, mix& gas) {
    if (line.rfind("--", 0) != 0) {  // doesn't start with "--"
        cerr << "Unrecognized option format: " << line << "\n";
        return;
    }

    // Find '='
    size_t eq = line.find('=');
    if (eq == string::npos) {
        cerr << "Expected --key=value format: " << line << "\n";
        return;
    }

    string key = line.substr(2, eq - 2);      // skip leading "--"
    string value = trimIO(line.substr(eq + 1)); // after '='

    if (key == "mode") {
        if (value == "sweep") {
            cfg.mode = RunMode::Sweep;
            cfg.modestring = value;
            cout << "-- Run mode: sweep\n\n";
        } else { // anything else (or "single")
            cfg.mode = RunMode::Single;
            cfg.modestring = value;
            cout << "-- Run mode: single\n\n";
        }
    }
    else if (key == "constraint"){
        if (value == "TP" || value == "UV" || value == "TV" || value == "CFD" ) {
            cfg.constraint = constraintFromString(value);
            cfg.constraintstring = value;
            cout << "-- Constraint set to: " << value;
            cout << "\n\n";
        }
        else cout << "== Invalid constraint" << endl << endl;
    }
    else if (key == "species") {
        // Optionally trim whitespace here if you have a trim() helper
        // value = trim(value);

        // Check if it's a single token (no comma) and one of the known presets
        bool is_preset =
            value.find(',') == std::string::npos &&
            (value == "air5" || value == "air7" || value == "air11" ||
            value == "air13" || value == "mars8");

        if (is_preset) {
            // Use predefined gas mix
            gas = mixFromString(value);
            cfg.creation = false;
            cfg.species.clear();
            cfg.elements = gas.elements;
            cfg.initial = gas.b;
            cfg.species = gas.speciesnames;

            cout << endl;

        } else {
            // Treat it as a user-specified list of species
            cfg.species = splitCommaList(value);
            cfg.creation = true;

            std::cout << "-- Species set to: ";
            for (const auto& s : cfg.species) std::cout << s << " ";
            std::cout << "\n\n";
        }
    }
    else if (key == "elements") {
        cfg.elements = splitCommaList(value);
        cout << "-- Elements set to: ";
        for (auto& e : cfg.elements) cout << e << " ";
        cout << "\n\n";
    }
    else if (key == "Y") {
        cfg.initial = splitCommaDoubles(value);
        cout << "-- Elemental mass fractions set to: ";
        for (auto& e : cfg.initial) cout << e << " ";
        cout << "\n\n";
    }
    else {
        cerr << "-- Unknown option key: " << key << "\n\n";
    }
}

// Just for debugging: show current configuration
void show_config(const Config& cfg) {
    cout << endl << endl << "=====: Current configuration :=====\n";
    cout << "-- mode: " << (cfg.modestring.empty() ? "(not set)" : cfg.modestring) << "\n";
    cout << "-- constraint: " << (cfg.constraintstring.empty() ? "(not set)" : cfg.constraintstring) << "\n";

    cout << "-- species: ";
    if (cfg.species.empty()) cout << "(none)";
    else {
        for (auto& s : cfg.species) cout << s << " ";
    }
    cout << "\n";

    cout << "-- elements: ";
    if (cfg.elements.empty())cout << "(none)";
    else {
        for (auto& e : cfg.elements) cout << e << " ";
    }
    cout << "\n";

    cout << "-- Initial mass fractions: ";
    if (cfg.initial.empty())cout << "(none)";
    else {
        for (auto& e : cfg.initial) cout << e << " ";
    }
    cout << endl << endl;
}

void run_sweep(double min, double max, double val2, Config cfg, mix gas, string filename_in, bool mass) {

    string filename = "../user-files/" + filename_in + ".dat";
    ofstream write(filename);
    if (!write) {
        cerr << "Failed to open " << filename << " for writing.\n";
        exit;
    }

    // Buffers for BLOCK output
    vector<double> Tvals;
    vector<double> Y;
    int N = 5000;

    CESolver CE(gas, cfg.constraint);

    for (int i = 0; i < N; ++i) {
        double v = min + (max - min) / (N - 1) * i;
        if (cfg.constraint == ConstraintType::CFD) {
            CE.CFD_equilibrium(v, val2);
        }
        else 
            CE.compute_equilibrium(v, val2);

        Tvals.push_back(gas.T);
        for (int j = 0; j < gas.NS; ++j) {
            if (mass == true)
                Y.push_back(gas.Y[j]);
            else 
                Y.push_back(gas.X[j]);
        }
        i++;
        if (gas.T > 19000) exit;
    }

    // ===== Tecplot ASCII Header (BLOCK format) =====
    // Title
    write << "TITLE = \"Equilibrium TP Sweep: " << gas.name << "\"\n";

    // Variables: T, P, and species mass fractions
    write << "VARIABLES = \"T [K]\"";
    for (int j = 0; j < gas.NS; ++j) {
        if (mass == true)
            write << ", \"Y(" << gas.species[j].name << ")\"";
        else    
            write << ", \"X(" << gas.species[j].name << ")\"";
    }
    write << "\n";

    // One 1D zone with N points
    write << "ZONE T=\"Sweep\", I=" << Tvals.size() << ", F=POINT\n";


    for (int i = 0; i < Tvals.size(); ++i) {
        write << Tvals[i];
        for (int j = 0; j < gas.NS; ++j) 
            write << ", " << Y[i * gas.NS + j];
        write << endl;
    }


    write.close();

    cout << "-- Tecplot file saved to " << filename << "\n";
}

// Stub for your actual minimization logic
void run_minimize(const Config& cfg, mix gas) {
    if (cfg.creation) 
        gas = create_mixture(cfg.species, cfg.elements, cfg.initial);
    
    CESolver CE(gas, cfg.constraint);
    double v1, v2, min, max;
    string filename;
    bool type;
    int type_in;

    switch (cfg.constraint) {
            case::ConstraintType::TP:
                switch (cfg.mode) {
                    case::RunMode::Single:
                        cout << "-- Enter temperature [K]: ";
                        cin >> v1;
                        cout << "-- Enter pressure [Pa]: ";
                        cin >> v2;
                        CE.compute_equilibrium(v1, v2);
                        print_properties(gas);
                    break;
                    case::RunMode::Sweep: 
                        cout << "-- Enter starting temperature [K]: ";
                        cin >> min;
                        cout << "-- Enter ending temperature [K]: ";
                        cin >> max;
                        cout << "-- Enter pressure [Pa]: ";
                        cin >> v2;
                        cout << "-- Enter 1 for plotting mass fraction, 0 for molar fractions ";
                        cin >> type_in;
                        if (type_in == 1)
                            type = true;
                        else 
                            type = false;
                        cout << "-- Enter filename (.dat is already appended, will be saved to the 'user-files' directory): ";
                        cin >> filename;
                        run_sweep(min, max, v2, cfg, gas, filename, type);
                    break;
                }
            break;
            case::ConstraintType::TV:
                switch (cfg.mode) {
                    case::RunMode::Single:
                        cout << "-- Enter temperature [K]: ";
                        cin >> v1;
                        cout << "-- Enter specific volume [m^3/kg]: ";
                        cin >> v2;
                        CE.compute_equilibrium(v1, v2);
                        print_properties(gas);
                    break;
                    case::RunMode::Sweep: 
                        cout << "-- Enter starting temperature [k]: ";
                        cin >> min;
                        cout << "-- Enter ending temperature [k]: ";
                        cin >> max;
                        cout << "-- Enter specific volume [m^3/kg]: ";
                        cin >> v2;
                        cout << "-- Enter 1 for plotting mass fraction, 0 for molar fractions ";
                        cin >> type_in;
                        if (type_in == 1)
                            type = true;
                        else 
                            type = false;
                        cout << "-- Enter filename (.dat is already appended, will be saved to the 'user-files' directory): ";
                        cin >> filename;
                        run_sweep(min, max, v2, cfg, gas, filename, type);
                    break;
                }
            break;
            case::ConstraintType::UV:
                switch (cfg.mode) {
                    case::RunMode::Single:
                        cout << "-- Enter internal energy [J/kg]: ";
                        cin >> v1;
                        cout << "-- Enter specific volume [m^3/kg]: ";
                        cin >> v2;
                        CE.compute_equilibrium(v1, v2);
                        print_properties(gas);
                    break;
                    case::RunMode::Sweep: 
                        cout << "-- Enter starting internal energy [J/kg]: ";
                        cin >> min;
                        cout << "-- Enter ending internal energy [J/kg]: ";
                        cin >> max;
                        cout << "-- Enter specific volume [m^3/kg]: ";
                        cin >> v2;
                        cout << "-- Enter 1 for plotting mass fraction, 0 for molar fractions ";
                        cin >> type_in;
                        if (type_in == 1)
                            type = true;
                        else 
                            type = false;
                        cout << "-- Enter filename (.dat is already appended, will be saved to the 'user-files' directory): ";
                        cin >> filename;
                        run_sweep(min, max, v2, cfg, gas, filename, type);
                    break;
                }
            break;
            case::ConstraintType::CFD:
                switch (cfg.mode) {
                    case::RunMode::Single:
                        cout << "-- Enter internal energy [J/kg]: ";
                        cin >> v1;
                        cout << "-- Enter density [kg/m^3]: ";
                        cin >> v2;
                        CE.CFD_equilibrium(v1, v2);
                        print_properties(gas);
                    break;
                    case::RunMode::Sweep: 
                        cout << "-- Enter starting internal energy [J/kg]: ";
                        cin >> min;
                        cout << "-- Enter ending internal energy [J/kg]: ";
                        cin >> max;
                        cout << "-- Enter density [kg/m^3]: ";
                        cin >> v2;
                        cout << "-- Enter 1 for plotting mass fraction, 0 for molar fractions ";
                        cin >> type_in;
                        if (type_in == 1)
                            type = true;
                        else 
                            type = false;
                        cout << "-- Enter filename (.dat is already appended, will be saved to the 'user-files' directory): ";
                        cin >> filename;
                        run_sweep(min, max, v2, cfg, gas, filename, type);
                    break;
                }
            break;
    }    
}

void show_options() {
    using std::cout;
    using std::endl;

    cout << "===============================================================================\n";
    cout << "  Chemical Equilibrium Solver - Command Line Options\n";
    cout << "===============================================================================\n";
    cout << "  General rules:\n";
    cout << "    * All options must be provided as '--option='\n";
    cout << "    * Lists use commas; spaces are optional.\n";
    cout << "      Example: --species=N2,O2,N  is the same as  --species= N2, O2, N\n";
    cout << "      Use species available in the NASA thermo.inp file. \n";
    cout << "      Species, elements, and Y lists are not checked, please make sure there are no errors.\n";
    cout << "-------------------------------------------------------------------------------\n\n";

    cout << "  Mixture definition\n";
    cout << "    --species=LIST       List of species to include\n";
    cout << "                         e.g. --species=N2,O2,NO,N,O\n";
    cout << "                         or a predefined mixture name:\n";
    cout << "                         [air5, air7, air11, air13, mars8]\n";
    cout << "                         e.g. --species= air11\n\n";

    cout << "    --elements=LIST      List of elements in the species list\n";
    cout << "                         e.g. --elements= N,O\n";
    cout << "                         (ignored when using a predefined mixture)\n\n";

    cout << "    --Y=LIST             Element mass fractions, in the same order\n";
    cout << "                         as given in --elements\n";
    cout << "                         e.g. --Y= 0.7572,0.2428\n";
    cout << "                         (ignored when using a predefined mixture)\n\n";

    cout << "    --constraint=TYPE    Energy constraint to minimize with respect to\n";
    cout << "                         Available: TP, TV, UV, CFD\n";
    cout << "                         e.g. --constraint= TP\n\n";

    cout << "    --mode=MODE          Run mode for the solver\n";
    cout << "                         Default (omit or anything else): single case\n";
    cout << "                         Use 'sweep' to run over a range for plotting\n";
    cout << "                         e.g. --mode= sweep\n";
    cout << "===============================================================================\n";
}



int main() {

    Config cfg;
    mix gas;

    cout << endl;
    cout << "=====: Interactive University of Minnesota chemical equilibrium solver :=====\n";
    cout << endl << endl << "Commands:\n";
    cout << "  help  - view options \n";
    cout << "  show  - show current configuration\n";
    cout << "  printNASA  - print NASA thermodynamic data for current mix\n";
    cout << "  run   - run minimization\n";
    cout << "  quit  - exit\n\n";

    string line;
    while (true) {
        cout << "> ";
        if (!getline(std::cin, line)) break;  // EOF

        line = trimIO(line);
        if (line.empty()) continue;

        if (line == "quit" || line == "exit") {
            break;
        }
        else if (line == "show") {
            show_config(cfg);
        }
        else if (line == "run") {
            run_minimize(cfg, gas);
        }
        else if (line == "help") {
            show_options();
        }
        else if (line == "printNASA") {
            print_NASA_mix(gas);
        }
        else if (line.rfind("--", 0) == 0) {
            // option line like --mode=...
            handle_option_line(line, cfg, gas);
        }
        else {
            std::cerr << "Unknown command: " << line << "\n";
        }
    }

    std::cout << "Exiting minimization program.\n";
    return 0;
}
