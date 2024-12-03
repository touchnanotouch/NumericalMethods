#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>

#include "../headers/Vector.h"
#include "../headers/NLE.h"
#include "../headers/SoNLE.h"


double func1(double* x) {
    return 2 * std::pow(x[0], 2) - x[0] * x[1] - 5 * x[0] + 1;
}

double func2(double* x) {
    return x[0] + 3 * std::log(x[0]) - std::pow(x[1], 2);
}

double func1_si(double* x) {
    return std::sqrt((x[0] * (x[1] + 5) - 1) / 2);
}

double func2_si(double* x) {
    return std::sqrt(x[0] + 3 * std::log(x[0]));
}

double func3(double* x) {
    return std::pow(x[0], 2) + std::pow(x[1], 2) - 1;
}

double func4(double* x) {
    return std::pow(x[0], 3) - x[1];
}

double func5(double* x) {
    return x[0] + x[1] - 3;
}

double func6(double* x) {
    return std::pow(x[0], 2) + std::pow(x[1], 2) - 9;
}


int main() {
    try {
        const int n = 2;

        std::function<double(double*)>* funcs = new std::function<double(double*)>[n];
        
        funcs[0] = &func1_si;
        funcs[1] = &func2_si;

        Vector<double> x0(n);
        x0[0] = 1;
        x0[1] = 5;

        SoNLE<double> sonle(n, funcs);

        std::cout << "si: " << sonle.solve_iter(x0, "si", 1000, 1e-10) << std::endl; 
        std::cout << "newton: " << sonle.solve_iter(x0, "newton", 1000, 1e-10) << std::endl;
        std::cout << "sign: " << sonle.solve_iter(x0, "sign", 1000, 1e-10) << std::endl;

        delete[] funcs;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
