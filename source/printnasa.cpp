#include "../libraries/CESolver/CESolver.h"

using namespace std;

int main(int argc, char* argv[]) {

    // vector to hold species names
    vector<string> species;

    // skip argv[0] (program name), start at 1
    for (int i = 1; i < argc; ++i) {
        species.emplace_back(argv[i]);
    }


    print_NASA(species);

    return 0;
}
