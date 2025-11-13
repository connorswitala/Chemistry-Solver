#ifndef COMMON_H
#define COMMON_H

#include<vector>
#include<string>
#include<cmath>
#include<iostream>
#include<array>
#include<stdexcept>
#include<cassert>
#include<iomanip>
#include<chrono>
#include<fstream>
#include <algorithm>

constexpr double gcon = 8314.462;       // Universal gas constant: Joules / mol-Kelvin
constexpr double bcon = 1.380649e-23;   // Boltzmann constant: Joules / Kelvin
constexpr double acon  = 6.022140e23;   // Avogradro's constant: 1/mol 

typedef std::vector<double> Vector;
#define NOW std::chrono::high_resolution_clock::now();


#endif