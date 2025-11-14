#include "mixes.h"


inline void print_NASA_mix(const mix& gas) {
    const int W_NAME   = 12;
    const int W_INT    = 8;
    const int W_DOUBLE = 20;

    // Header
    cout << left  << setw(W_NAME)   << "Species"
         << right << setw(W_INT)    << "Charge"
         << right << setw(W_DOUBLE) << "MW [g/mol]"
         << right << setw(W_DOUBLE) << "href [J/kmol]"
         << right << setw(W_DOUBLE) << "hf [J/kmol]\n";

    // Separator
    cout << std::string(W_NAME + W_INT + 3*W_DOUBLE, '-') << "\n";

    // Rows
    cout << fixed << setprecision(6);
    for (int j = 0; j < gas.NS; ++j) {
        const auto& sp = gas.species[j];

        cout << left  << setw(W_NAME)   << sp.name
             << right << setw(W_INT)    << sp.q
             << right << setw(W_DOUBLE) << sp.mw
             << right << setw(W_DOUBLE) << sp.href
             << right << setw(W_DOUBLE) << sp.hf
             << "\n";
    }
}

void print_NASA(vector<string> species_names) {

    int size = species_names.size();
    vector<SpeciesInfo> SI(size);

    for (int i = 0; i < size; ++i)
        SI[i] = read_species_info(species_names[i]);

    const int W_NAME   = 12;
    const int W_INT    = 8;
    const int W_DOUBLE = 20;

    // Header
    cout << left  << setw(W_NAME)   << "Species"
         << right << setw(W_INT)    << "Charge"
         << right << setw(W_DOUBLE) << "MW [g/mol]"
         << right << setw(W_DOUBLE) << "href [J/mol]"
         << right << setw(W_DOUBLE) << "hf [J/mol]\n";

    // Separator
    cout << string(W_NAME + W_INT + 3*W_DOUBLE, '-') << "\n";

    // Rows
    cout << fixed << setprecision(6);
    for (int j = 0; j < size; ++j) {
        cout << left  << setw(W_NAME)   << SI[j].name
             << right << setw(W_INT)    << SI[j].q
             << right << setw(W_DOUBLE) << SI[j].mw
             << right << setw(W_DOUBLE) << SI[j].href
             << right << setw(W_DOUBLE) << SI[j].hf
             << "\n";
    }
}


namespace common {

    mix air_mixture(GasType gastype) {
        mix gas;
        switch (gastype) {
            case::GasType::AIR5:
                gas = common::make_air5();
                return gas;
                break;

            case::GasType::AIR7:
                gas = common::make_air7();
                return gas;
                break;    

            case::GasType::AIR11:
                gas = common::make_air11();
                return gas;
                break;

            case::GasType::AIR11_AR:
                gas = common::make_air11_Ar();
                return gas;
                break;

            case::GasType::AIR13:
                gas = common::make_air13();
                return gas;
                break;
        }
    }

    mix make_air5() {

        vector<string> species = {"N2", "O2", "NO", "N", "O"};        
        vector<string> elements = {"N", "O"};
        vector<double> initial = {0.7572, 0.2428}; 

        mix gas = create_mixture(species, elements, initial);

        return gas;
    }

    mix make_air7() {
        vector<string> species = {"N2", "O2", "NO", "N", "O", "NO+", "e-"};        
        vector<string> elements = {"N", "O"};
        vector<double> initial = {0.7572, 0.2428}; 

        mix gas = create_mixture(species, elements, initial);

        return gas;
    }

    mix make_air11_Ar() {
        

        vector<string> species = {"N2", "O2", "NO", "N", "O", "Ar", "NO+", "N+", "O+", "Ar+", "e-"};        
        vector<string> elements = {"N", "O", "Ar"};
        vector<double> initial = { 0.755, 0.232, 0.0129}; 

        mix gas = create_mixture(species, elements, initial);

        return gas;
    }

    mix make_air11() {

        vector<string> species = {"N2", "O2", "NO", "N", "O", "N2+", "O2+", "NO+", "N+", "O+", "e-"};        
        vector<string> elements = {"N", "O"};
        vector<double> initial = {0.7572, 0.2428}; 

        mix gas = create_mixture(species, elements, initial);

        return gas;       
    }

    mix make_air13() {
        vector<string> species = {"N2", "O2", "NO", "N", "O", "Ar", "N2+", "O2+", "NO+", "N+", "O+", "Ar+", "e-"};        
        vector<string> elements = {"N", "O", "Ar"};
        vector<double> initial = { 0.755, 0.232, 0.0129}; 

        mix gas = create_mixture(species, elements, initial);

        return gas;
    }
}


mix create_mixture(vector<string>& speciesNames, vector<string>& elementNames, vector<double>& initial_Y) {
    mix gas;

    gas.name = "Personal Mix";

    gas.NS = speciesNames.size();
    gas.species = vector<SpeciesInfo>(gas.NS);

    // Read in species info (mw, charge, href, etc.)
    for (int j = 0; j < gas.NS; ++j) 
        gas.species[j] = read_species_info(speciesNames[j]);

    for (int j = 0; j < gas.NS; ++j)
        gas.species[j].name = speciesNames[j];

    // Find elements used in species and build 'a' matrix.
    build_elements_and_a(gas);
    gas.elemental_mw = vector<double>(gas.NE);

    for (int i = 0; i < gas.NE; ++i)  {
        for (int j = 0; j < gas.NS; ++j) {
            if (gas.elements[i] == gas.species[j].name) 
                gas.elemental_mw[i] = gas.species[j].mw;
        }
    }

    
    for (int i = 0; i < gas.NE; ++i) {
        for (int j = 0; j < gas.NS; ++j) {
            gas.a[i * gas.NS + j] *= gas.elemental_mw[i];
        }
    }

    gas.b = vector<double>(gas.NE, 0.0);   
    for (int i = 0; i < gas.NE; ++i) {
        if (gas.elements[i] == elementNames[i])
            gas.b[i] = initial_Y[i];            
    }
    
    gas.HAS_IONS = false;
    for (int j = 0; j < gas.NS; ++j) {
        if (gas.species[j].name == "e-")
            gas.HAS_IONS = true;

        gas.species[j].href *= 1000.0;
        gas.species[j].hf *= 1000.0;
    }        

    gas.H0_RT = vector<double>(gas.NS);
    gas.U0_RT = vector<double>(gas.NS);
    gas.S0_R = vector<double>(gas.NS);
    gas.mu0_RT = vector<double>(gas.NS);
    gas.CP0_R = vector<double>(gas.NS);
    gas.mu_RT = vector<double>(gas.NS);        
    gas.N = vector<double>(gas.NS);

    gas.Y = vector<double>(gas.NS);
    gas.X = vector<double>(gas.NS);
    gas.X0 = vector<double>(gas.NS);

    for (int j = 0; j < gas.NS; ++j) {
        if (gas.species[j].q > 0.1) 
            gas.X0[j] = 1e-12;
        else
            gas.X0[j] = 0.1 / gas.NS;
    }

    
    gas.T = 3800.0;         

    cout << endl << "-- " << gas.NS << " species gas mix created. Contains: ";
    for (int i = 0; i < gas.NS; ++i) cout << gas.species[i].name << " ";
    cout << endl;

    cout << endl;
    return gas;    
}


