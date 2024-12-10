#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>

#include "../headers/Vector.h"
#include "../headers/NLE.h"
#include "../headers/SoNLE.h"


double func1_1(double* x) {
    return 2 * x[0] * x[0] - x[0] * x[1] - 5 * x[0] + 1;
}

double func1_2(double* x) {
    return x[0] + 3 * std::log10(x[0]) - x[1] * x[1];
}

//

double func2_1(double* x) {
    return x[0] + x[1] - 3;
}

double func2_2(double* x) {
    return x[0] * x[0] + x[1] * x[1] - 9;
}

//

double func3_1(double* x) {
    return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;
}

double func3_2(double* x) {
    return 2 * x[0] * x[0] + x[1] * x[1] - 4 * x[2];
}

double func3_3(double* x) {
    return 3 * x[0] * x[0] - 4 * x[1] + x[2] * x[2];
}

//

double func4_1(double* x) {
    return x[0] + x[0] * x[0] - 2 * x[1] * x[2] - 0.1;
}

double func4_2(double* x) {
    return x[1] - x[1] * x[1] + 3 * x[0] * x[2] + 0.2;
}

double func4_3(double* x) {
    return x[2] + x[2] * x[2] + 2 * x[0] * x[1] - 0.3;
}


int main() {
    try {
        const int n = 3;

        std::function<double(double*)>* funcs = new std::function<double(double*)>[n];
        
        // funcs[0] = &func1_1;
        // funcs[1] = &func1_2;
        
        // funcs[0] = &func2_1;
        // funcs[1] = &func2_2;
        
        funcs[0] = &func3_1;
        funcs[1] = &func3_2;
        funcs[2] = &func3_3;

        // funcs[0] = &func4_1;
        // funcs[1] = &func4_2;
        // funcs[2] = &func4_3;

        Vector<double> x0(n);
        // x0[0] = 3.5;
        // x0[1] = 2.2;    
        
        // x0[0] = 1;
        // x0[0] = 5; 
        
        x0[0] = 0.5;
        x0[1] = 0.5;
        x0[2] = 0.5;

        // x0[0] = 0;
        // x0[1] = 0;
        // x0[2] = 0;

        SoNLE<double> sonle(n, funcs);

        std::cout << "si: " << sonle.solve_iter(x0, "si", 1000, 1e-10) << std::endl; 
        std::cout << "newton: " << sonle.solve_iter(x0, "newton", 1000, 1e-10) << std::endl;
        std::cout << "sign: " << sonle.solve_iter(x0, "sign", 1000, 1e-10) << std::endl;
        std::cout << "grad: " << sonle.solve_iter(x0, "grad", 1000, 1e-10) << std::endl;
        std::cout << "test1: " << sonle.solve_iter(x0, "test1", 1000, 1e-10) << std::endl;
        std::cout << "test2: " << sonle.solve_iter(x0, "test2", 1000, 1e-10) << std::endl;
        std::cout << "test3: " << sonle.solve_iter(x0, "test3", 1000, 1e-10) << std::endl;

        delete[] funcs;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
