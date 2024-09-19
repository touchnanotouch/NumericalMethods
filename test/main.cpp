#include <iostream>

#include "../src/headers/Matrix.h"
#include "../src/headers/Vector.h"

#include "../src/headers/SoLE.h"


int main() {
    try {
        SoLE<int> sole(2, 2);
        
        Matrix<int> mat = sole.get_matrix();
        Vector<int> vec = sole.get_vector();

        mat[0][0] = 2;
        mat[0][1] = -1;
        mat[1][0] = 5;
        mat[1][1] = -3;

        vec[0] = 20;
        vec[1] = 51;

        Vector<int> ans(2);
        ans[0] = 9;
        ans[1] = -2;

        sole.set_matrix(mat);
        sole.set_vector(vec);

        std::cout << sole.is_solution(ans) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
