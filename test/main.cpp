#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "../headers/SoLE.h"
#include "../headers/Vector.h"
#include "../headers/Tester.h"


int main() {
    try {
        const int n_count = 9;

        int* n_list = new int[n_count]{10, 25, 50, 75, 100, 250, 500, 750, 1000};

        for (size_t i = 0; i < n_count; i++) {
            int n = n_list[i];

            std::cout << "n = " << n << " . . . ";

            std::string m_path = "../data/gens/matrix_n" + std::to_string(n) + ".txt";
            std::string v_path = "../data/gens/vector_n" + std::to_string(n) + ".txt";
            std::string s_path = "../data/gens/solution_n" + std::to_string(n) + ".txt";

            std::string i_path = "../data/info_all_n" + std::to_string(n) + ".txt";

            Tester<double> tester(n, m_path, v_path);

            Vector<double> solution(n);
            solution.set_vector(s_path);

            // SoLE<double> sole(n, n);
            // sole.set_matrix(m_path);
            // sole.set_vector(v_path);

            // std::cout << "\n\n\n\n" << sole.solve_iter("si") << std::endl;
            // std::cout << sole.solve_iter("seidel") << std::endl;

            std::stringstream ss;

            ss << "n = " << n << std::endl;
            ss << tester.measure_all(solution);

            std::ofstream i_file(i_path);
            i_file << ss.str();
            i_file.close();

            std::cout << "done" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
