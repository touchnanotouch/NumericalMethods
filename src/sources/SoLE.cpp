#include <string>
#include <sstream>
#include <fstream>

#include <omp.h>

#include "../headers/SoLE.h"


template class SoLE<int>;
template class SoLE<float>;
template class SoLE<double>;


template<typename T>
void SoLE<T>::elimination(

) {
    size_t n = row_count();

    Matrix<T> mat_copy(n, n);
    mat_copy.set_matrix(this->matrix().mat());

    Vector<T> vec_copy(n);
    vec_copy.set_vector(this->vector().vec());

    #pragma omp parallel for
    for (int k = 0; k < n; k++) {
        T diag_val = mat_copy[k][k];

        for (int i = k + 1; i < n; i++) {
            double mu = mat_copy[i][k] / diag_val;
            for (int j = 0; j < n; j++)
                mat_copy[i][j] -= mat_copy[k][j] * mu;

            vec_copy[i] -= vec_copy[k] * mu;
        }

        for (int i = 0; i < k; i++) {
            double mu = mat_copy[i][k] / diag_val;
            for (int j = 0; j < n; j++) {
                mat_copy[i][j] -= mat_copy[k][j] * mu;
            }

            vec_copy[i] -= vec_copy[k] * mu;
        }
    }

    set_matrix(mat_copy);
    set_vector(vec_copy);
}

template<typename T>
Vector<T> SoLE<T>::calc_next_x(
    Vector<T> x,
    std::string method
) const {
    int n = (int)row_count();

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        T sum = 0;
        if (method == "si") {
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += matrix()[i][j] * x[j];
                }
            }
        } else if (method == "seidel") {
            for (int j = 0; j < i; j++) {
                sum += matrix()[i][j] * x[j];
            }
            for (int j = i + 1; j < n; j++) {
                sum += matrix()[i][j] * x[j];
            }
        }

        x[i] = (vector()[i] - sum) / matrix()[i][i];
    }

    return x;
}

template<typename T>
size_t SoLE<T>::row_count(
    
) const {
    return _row_cnt;
}

template<typename T>
size_t SoLE<T>::col_count(

) const {
    return _col_cnt;
}

template<typename T>
Matrix<T> SoLE<T>::matrix(

) const {
    return _coeffs;
}

template<typename T>
Vector<T> SoLE<T>::vector(
    
) const {
    return _terms;
}

template<typename T>
void SoLE<T>::set_row_count(
    size_t row_count
) {
    _row_cnt = row_count;
}

template<typename T>
void SoLE<T>::set_col_count(
    size_t col_count
) {
    _col_cnt = col_count;
}

template<typename T>
void SoLE<T>::set_matrix(
    Matrix<T>& matrix
) {
    _coeffs = matrix;
}

template<typename T>
void SoLE<T>::set_matrix(
    std::string file_path,
    const char delimiter
) {
    Matrix<T> result(row_count(), col_count());

    std::ifstream file(file_path);
    std::string line;

    for (int i = 0; i < row_count(); i++) {
        std::getline(file, line);
        std::istringstream stream(line);
        T value;

        Vector<T> row(col_count());

        for (int j = 0; j < col_count(); j++) {
            stream >> value;
            row[j] = value;
            if (stream.peek() == delimiter) {
                stream.ignore();
            }
        }

        result.set_row(row.vec(), i);
    }

    set_matrix(result);

    file.close();
}

template <typename T>
void SoLE<T>::set_vector(
    Vector<T>& vector
) {
    _terms = vector;
}

template<typename T>
void SoLE<T>::set_vector(
    std::string file_path,
    const char delimiter
) {
    Vector<T> result(row_count());

    std::ifstream file(file_path);
    std::string line;

    for (int i = 0; i < row_count(); i++) {
        std::getline(file, line);
        std::istringstream stream(line);
        T value;
    
        for (int j = 0; j < col_count(); j++) {
            stream >> value;
            result[j] = value;
            if (stream.peek() == delimiter) {
                stream.ignore();
            }
        }
    }

    set_vector(result);

    file.close();
}

template<typename T>
Matrix<T> SoLE<T>::extended(
    
) {
    size_t n = row_count();
    size_t m = col_count();

    // TODO : rework for loops so it doesn't look like unlogical shit

    Matrix<T> result(n, m + 1);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            result[i][j] = matrix()[i][j];
        }
        result[i][m] = vector()[i];
    }

    return result;
}

template<typename T>
double SoLE<T>::det(

) {
    return matrix().det();
}

template<typename T>
int SoLE<T>::rank(

) {
    return matrix().rank();
}

template<typename T>
bool SoLE<T>::is_compatible(

) {
    return extended().rank() == matrix().rank();
}

template<typename T>
bool SoLE<T>::is_solution(
    Vector<T>& solution
) const {
    return matrix() * solution == vector();
}

template<typename T>
bool SoLE<T>::is_diag_d(

) {
    return matrix().is_diag_d();
}

template<typename T>
bool SoLE<T>::is_solvable(
    
) {
    return is_compatible() && is_diag_d();
}

template<typename T>
void SoLE<T>::to_diag_d(

) {
    elimination();

    // int n = (int)row_count();

    // Matrix<T> E(n, n);
    // for (int i = 0; i < n; i++) {
    //     E[i][i] = (T)1;
    // }

    // Matrix<T> mat_new = matrix().inversed() * E;
    
    // set_matrix(mat_new);
}

template<typename T>
Vector<T> SoLE<T>::solve_iter(
    std::string method,
    int iter_max,
    double epsilon
) const {
    Vector<T> x(row_count());

    for (int i = 0; i < iter_max; i++) {
        Vector<T> x_next = calc_next_x(x, method);

        double x_norm = (x_next - x).norm();

        std::cout << "iter " << i + 1 << ", norm: " << x_norm << std::endl;

        if (x_norm < epsilon) {
            break;
        } else {

        }

        x = x_next;
    }

    return x;
}
