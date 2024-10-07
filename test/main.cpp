#include <iostream>
#include <string>

#include "../headers/SoLE.h"


int main() {
    try {
        const int n = 4;

        std::string m_path = "../data/matrix_n" + std::to_string(n) + ".txt";
        std::string v_path = "../data/vector_n" + std::to_string(n) + ".txt";

        SoLE<double> sole(n, n);
        sole.set_matrix(m_path);
        sole.set_vector(v_path);

        std::cout << sole << std::endl;

        Vector<double> solution = sole.solve("lu");

        // if (!sole.is_compatible()) {
        //    throw std::runtime_error("Matrix isn't compatible");
        // } else {
        //    if (!sole.is_diag_d()) {
        //        sole.to_diag_d();
        //    }

        //    if (!sole.is_diag_d()) {
        //        throw std::runtime_error("Matrix isn't transformable to d/d form");
        //    }
        // }
        // Vector<double> solution = sole.solve_iter("seidel");

        std::cout << solution << std::endl;
        std::cout << "is_solution: " << sole.is_solution(solution) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
