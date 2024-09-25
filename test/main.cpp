#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "../headers/Matrix.h"
#include "../headers/Vector.h"

#include "../headers/SoLE.h"


int main() {
    try {
        SoLE<double> sole(4, 4);
        sole.set_matrix("../data/matrix_n4.txt");
        sole.set_vector("../data/vector_n4.txt");

        if (!sole.is_compatible()) {
            throw std::runtime_error("Matrix isn't compatible");
        } else {
            if (!sole.is_diag_d()) {
                sole.to_diag_d();
            }

            if (!sole.is_solvable()) {
                throw std::runtime_error("Matrix isn't solvable");
            }
        }

        std::cout << sole.solve_iter("seidel") << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
