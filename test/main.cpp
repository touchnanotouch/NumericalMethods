#include <iostream>
#include <string>

#include "../headers/SoLE.h"


int main() {
    try {
        const int n = 3;

        std::string m_path = "../data/matrix_n" + std::to_string(n) + ".txt";
        std::string v_path = "../data/vector_n" + std::to_string(n) + ".txt";

        SoLE<double> sole(n, n);
        sole.set_matrix(m_path);
        sole.set_vector(v_path);

        std::cout << sole << std::endl;

        Vector<double> thomas_sol = sole.solve("thomas");
        Vector<double> lu_sol = sole.solve("lu");
        Vector<double> qroots_sol = sole.solve("qroots");
        
        Vector<double> si_sol = sole.solve_iter("si");
        Vector<double> seidel_sol = sole.solve_iter("seidel");
        Vector<double> sor_sol = sole.solve_iter("sor");
        Vector<double> res_sol = sole.solve_iter("res");
        Vector<double> grad_sol = sole.solve_iter("grad");

        std::cout << "Exact methods\n" << std::endl;
        
        std::cout << "thomas: " << thomas_sol;
        std::cout << "lu: " << lu_sol;
        std::cout << "qroots: " << qroots_sol << std::endl;

        std::cout << "Iterative methods\n" << std::endl;

        std::cout << "si: " << si_sol;
        std::cout << "seidel: " << seidel_sol;
        std::cout << "sor: " << sor_sol;
        std::cout << "res: " << res_sol;
        std::cout << "grad: " << grad_sol;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
