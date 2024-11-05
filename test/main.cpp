#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>

#include "../headers/NLE.h"


double func(double x) {
    return std::pow(x, 3) - x + 1;
}


int main() {
    try {
        NLE<double> eq(func);

        std::cout << eq.solve_iter(0, -2, 1, "si") << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
