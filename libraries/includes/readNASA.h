#ifndef READNASA_H
#define READNASA_H

#include "common.h"
#include "thermoObjs.h"
#include <cctype>
#include <regex>
#include <sstream>
#include <map>
#include <unordered_map>

using namespace std;


inline string trim(const string& s) {
    size_t first = 0;
    while (first < s.size() &&
           isspace(static_cast<unsigned char>(s[first]))) ++first;

    ssize_t last = s.size();
    while (last > first &&
           isspace(static_cast<unsigned char>(s[last - 1]))) --last;

    return s.substr(first, last - first);
}

inline double fortran_to_double(string s) {
    std::replace(s.begin(), s.end(), 'D', 'E');
    std::replace(s.begin(), s.end(), 'd', 'e');
    try {
        return std::stod(s);
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("Cannot parse numeric value from field '" + s + "'");
    }
}

inline int compute_charge_from_name(const string& name) {
    int q = 0;
    for (char c : name) {
        if (c == '+') q++;
        else if (c == '-') q--;
    }
    return q;
}

inline vector<double> parse_poly_block(const string& line1, const string& line2, int max_coeffs) {
    string block = line1 + " " + line2;

    // Match Fortran-style reals like -7.453750000D+02 or 1.234567890E-01
    std::regex num_re(R"(([+-]?\d+\.\d+(?:[DdEe][+-]\d+)))");

    vector<double> coeffs;
    coeffs.reserve(max_coeffs);

    for (std::sregex_iterator it(block.begin(), block.end(), num_re),
                              end;
         it != end; ++it)
    {
        string s = (*it)[1].str();
        coeffs.push_back(fortran_to_double(s));
        if ((int)coeffs.size() >= max_coeffs)
            break;
    }

    if ((int)coeffs.size() != max_coeffs) {
        throw std::runtime_error(
            "Expected " + std::to_string(max_coeffs) +
            " polynomial coefficients, got " + std::to_string(coeffs.size()));
    }

    return coeffs;
}

inline SpeciesInfo read_species_info(const std::string& target_name) {

    string filename = "../misc/files/thermo.inp";

    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("Could not open thermo file: " + filename);
    }

    SpeciesInfo sp{};
    sp.name = target_name;
    sp.q    = compute_charge_from_name(target_name);
    sp.poly.clear();

    std::string line;

    // ===== 1. Find species header line (e.g. "Ag-" or "e-") =====
    bool found = false;
    while (std::getline(fin, line)) {
        if (line.size() >= target_name.size() &&
            line.compare(0, target_name.size(), target_name) == 0 &&
            (line.size() == target_name.size() ||
             std::isspace(static_cast<unsigned char>(line[target_name.size()]))))
        {
            found = true;
            break;
        }
    }

    if (!found) {
        throw std::runtime_error("Species '" + target_name + "' not found in file.");
    }

    // ===== 2. Next line: MW and Hf(298) at the end =====
    if (!std::getline(fin, line)) {
        throw std::runtime_error("Unexpected EOF after species line for " + target_name);
    }

    {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string tok;
        while (iss >> tok) tokens.push_back(tok);

        if (tokens.size() < 2) {
            throw std::runtime_error("Not enough tokens on MW/Hf line for " + target_name);
        }

        sp.hf = fortran_to_double(tokens.back());                 // Hf(298)
        sp.mw = fortran_to_double(tokens[tokens.size() - 2]);     // MW
    }

    // ===== 3. NASA-9 polynomials: 3 ranges × 9 coeffs =====
    const int NRANGES    = 3;
    const int NCOEFF_PER = 9;
    sp.poly.reserve(NRANGES * NCOEFF_PER);

    for (int r = 0; r < NRANGES; ++r) {
        // --- Line A: T_low, T_high, ..., href (last field) ---
        if (!std::getline(fin, line)) {
            throw std::runtime_error("Unexpected EOF reading T-range line for " + target_name);
        }

        if (r == 0) {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::string tok;
            while (iss >> tok) tokens.push_back(tok);
            if (!tokens.empty()) {
                sp.href = fortran_to_double(tokens.back());   // e.g. 6197.428
            } else {
                throw std::runtime_error("Could not parse href for " + target_name);
            }
        }

        // --- Lines B & C: 9 coeffs total (a1..a9) ---
        std::string lineB, lineC;
        if (!std::getline(fin, lineB)) {
            throw std::runtime_error("Unexpected EOF reading coeff line 1 for " + target_name);
        }
        if (!std::getline(fin, lineC)) {
            throw std::runtime_error("Unexpected EOF reading coeff line 2 for " + target_name);
        }

        std::vector<double> coeffs = parse_poly_block(lineB, lineC, NCOEFF_PER);

        sp.poly.insert(sp.poly.end(), coeffs.begin(), coeffs.end());
    }

    return sp;
}

inline static map<string, int> parse_species_formula(const string& name) {

    map<string, int> elem_counts;

    size_t i = 0;
    const size_t n = name.size();

    while (i < n) {
        // Look for element start: uppercase letter
        if (!isupper(static_cast<unsigned char>(name[i]))) {
            ++i;
            continue;
        }

        // Element symbol: capital + optional lowercase
        string elem;
        elem.push_back(name[i++]); // add uppercase

        if (i < n && islower(static_cast<unsigned char>(name[i]))) {
            elem.push_back(name[i++]); // add lowercase, e.g. 'Ag'
        }

        // Now parse digits for count (if any)
        int count = 0;
        while (i < n && isdigit(static_cast<unsigned char>(name[i]))) {
            count = 10 * count + (name[i] - '0');
            ++i;
        }
        if (count == 0) count = 1;

        elem_counts[elem] += count;
    }

    return elem_counts;
}

inline void build_elements_and_a(mix& gas) {
    const int NS = gas.NS;
    if (NS <= 0 || static_cast<int>(gas.species.size()) < NS) {
        throw runtime_error("mix.NS or mix.species not consistent");
    }

    // 1) Parse each species and gather unique elements
    vector<map<string, int>> species_elem_counts(NS);
    vector<string> element_list;
    unordered_map<string, int> elem_to_index;

    for (int j = 0; j < NS; ++j) {
        const string& spname = gas.species[j].name;
        auto counts = parse_species_formula(spname);

        species_elem_counts[j] = move(counts);

        for (const auto& kv : species_elem_counts[j]) {
            const string& elem = kv.first;
            if (elem_to_index.find(elem) == elem_to_index.end()) {
                int idx = static_cast<int>(element_list.size());
                element_list.push_back(elem);
                elem_to_index[elem] = idx;
            }
        }
    }
    

    // 2) Set NE and resize mix.elements and mix.a
    const int NE = static_cast<int>(element_list.size());
    gas.NE = NE;

    gas.elements = vector<string>(NE);
    gas.elements = element_list;              // copy element names
    gas.a = vector<double>(NE * NS, 0.0);

    // 3) Fill a[i*NS + j] = number of atoms of element i in species j
    for (int j = 0; j < NS; ++j) {
        for (const auto& kv : species_elem_counts[j]) {
            const std::string& elem = kv.first;
            int count               = kv.second;
            int i                   = elem_to_index[elem];      // row index for this element
            gas.a[static_cast<std::size_t>(i) * NS + j] = static_cast<double>(count);
        }
    }
}


#endif