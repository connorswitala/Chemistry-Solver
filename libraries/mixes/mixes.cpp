#include "mixes.h"


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

mix make_mars8() {

    vector<string> species = {"CO2", "N2", "O2", "CO", "O", "C", "NO", "N"};        
    vector<string> elements = {"C", "O", "N"};
    vector<double> initial = {0.264, 0.7186, 0.0174}; 

    mix gas = create_mixture(species, elements, initial);

    return gas;
}


// Creates common mixture from enum class type
mix common_mixture(GasType gastype) {
    mix gas;
    switch (gastype) {
        case::GasType::AIR5:
            gas = make_air5();
            gas.name = "AIR5";
            return gas;
            break;

        case::GasType::AIR7:
            gas = make_air7();
            gas.name = "AIR7";
            return gas;
            break;    

        case::GasType::AIR11:
            gas = make_air11();
            gas.name = "AIR11";
            return gas;
            break;

        case::GasType::AIR11_AR:
            gas = make_air11_Ar();
            gas.name = "AIR11_AR";
            return gas;
            break;

        case::GasType::AIR13:
            gas = make_air13();
            gas.name = "AIR13";
            return gas;
            break;

        case::GasType::MARS8:
            gas = make_mars8();
            gas.name = "MARS8";
            return gas;
            break;
    }
}

// Create mixture from user-specified inputs
mix create_mixture(vector<string> speciesNames, vector<string> elementNames, vector<double> initial_Y) {
    mix gas;

    gas.name = "personal_mix";

    gas.NS = speciesNames.size();
    gas.species = vector<SpeciesInfo>(gas.NS);
    gas.speciesnames = speciesNames;

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

    gas.N_tot = 0.0;
    for (int j = 0; j < gas.NS; ++j) {
        if (gas.species[j].q > 0.1) 
            gas.X0[j] = 1e-30;
        else
            gas.X0[j] = 0.1 / gas.NS;
        gas.N_tot += gas.X0[j];
    }
    
    gas.T = 3800.0;         

    cout << endl << "-- " << gas.NS << " species gas mix created. Contains: ";
    for (int i = 0; i < gas.NS; ++i) cout << gas.species[i].name << " ";
    cout << endl;
    return gas;    
}


// ===== For printing to terminal =====

void print_NASA_mix(const mix& gas) {
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

    for (int j = 0; j < gas.NS; ++j) {
        cout << "NASA polynomials for " << gas.species[j].name << endl << endl;
        for (int i = 0; i < 3; ++i) {
            cout << "Temperature range " << i << endl;
            for (int k = 0; k < NCOEF; ++k) 
                cout << scientific << setprecision(10) << gas.species[j].poly[i * NCOEF + k] << "\t";
            cout << endl;
        }            
        cout << endl << endl;
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

void print_properties(mix& gas) {

    struct Prop {
        string label;
        double value;
        string unit;
    };

    // Left column: normal decimal (fixed)
    vector<Prop> leftProps = {
        {"Mixture Temperature:",     gas.T,          "[K]"},
        {"Mixture Density:",         gas.rho,        "[kg/m^3]"},
        {"Mixture Specific Volume:", gas.V,          "[m^3/kg]"},
        {"Mixture MW:",              gas.MW,         "[g/mol]"},
        {"Mixture R:",               gas.R,          "[J/kg-K]"},
        {"Mixture cp:",              gas.cp,         "[J/kg-K]"},
        {"Mixture cv:",              gas.cv,         "[J/kg-K]"},
        {"Mixture a:",               gas.c,          "[m/s]"},
        {"Mixture gamma:",           gas.gamma,      ""}
    };

    // Right column: scientific
    vector<Prop> rightProps = {
        {"Mixture Pressure:",        gas.p, "[Pa]"},
        {"Mixture Internal Energy:", gas.up,"[J/kg]"},
        {"Mixture dpdr:",            gas.dpdr,       ""},
        {"Mixture dtdr:",            gas.dtdr,       ""}
    };

    const int labelW = 24;
    const int valueW = 14;
    const int unitW  = 10;

    cout << "\n";

    size_t nRows = max(leftProps.size(), rightProps.size());

    for (size_t i = 0; i < nRows; ++i) {

        // ----- Left column: fixed -----
        if (i < leftProps.size()) {
            const auto& a = leftProps[i];

            cout << left  << setw(labelW) << a.label;

            cout << fixed << setprecision(5);
            cout << right << setw(valueW) << a.value << " ";

            cout << left  << setw(unitW)  << a.unit;
        } else {
            // no left entry, just pad
            cout << setw(labelW + valueW + unitW + 1) << " ";
        }

        cout << "   "; // spacing between columns

        // ----- Right column: scientific -----
        if (i < rightProps.size()) {
            const auto& b = rightProps[i];

            cout << left  << setw(labelW) << b.label;

            cout << scientific << setprecision(5);
            cout << right << setw(valueW) << b.value << " ";

            cout << left  << setw(unitW)  << b.unit;
        }

        cout << "\n";
    }

    // Reset for species table
    cout.unsetf(ios::floatfield);
    cout << fixed << setprecision(6);

    // Species table: name, mass fraction, molar fraction
    cout << "\n"
         << "================= Species Concentrations =================\n";

    const int nameW  = 12;
    const int fracW  = 14;

    // Header
    cout << left  << setw(nameW) << "Species"
         << right << setw(fracW) << "Y (mass)"
         << right << setw(fracW) << "X (molar)"
         << "\n";

    // Underline
    cout << left  << setw(nameW) << string(nameW - 1, '-')
         << right << setw(fracW) << string(fracW - 1, '-')
         << right << setw(fracW) << string(fracW - 1, '-')
         << "\n";

    // Rows
    for (int i = 0; i < gas.NS; ++i) {
        cout << left  << setw(nameW) << gas.species[i].name
             << right << setw(fracW) << gas.Y[i]
             << right << setw(fracW) << gas.X[i]
             << "\n";
    }

    cout << endl;
}

mix mixFromString(string& s) {

    mix g;
    GasType gt;

    if (s == "air5") {
        gt = GasType::AIR5;
        g = common_mixture(gt);
    }
    else if (s == "air7") {
        gt = GasType::AIR7;
        g = common_mixture(gt);
    }      
    else if (s == "air11") {
        gt = GasType::AIR11;
        g = common_mixture(gt);
    }       
    else if (s == "air13") {
        gt = GasType::AIR13;
        g = common_mixture(gt);
    }        
    else if (s == "mars8") {
        gt = GasType::MARS8;
        g = common_mixture(gt);
    }

    return g;       
}