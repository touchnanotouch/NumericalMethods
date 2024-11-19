#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>

#include "../headers/NLE.h"
//#include "../headers/SoNLE.h"


int main() {
    try {
        NLE<int> nle(4);
        nle.set_func_vector("../data/vector.txt");

        std::cout << nle.solve_iter(-2, -2, 1, "chord") << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
