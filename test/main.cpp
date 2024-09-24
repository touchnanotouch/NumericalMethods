#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "../headers/Matrix.h"
#include "../headers/Vector.h"

#include "../headers/SoLE.h"


int main() {
    try {
        SoLE<double> sole(3, 3);
        sole.set_matrix("../data/matrix.txt");
        sole.set_vector("../data/vector.txt");

        std::cout << sole << std::endl;

        std::cout << sole.solve_iter("seidel") << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
