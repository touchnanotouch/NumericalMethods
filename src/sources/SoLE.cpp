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
    for (size_t k = 0; k < n; k++) {
        T diag_val = mat_copy[k][k];

        for (size_t i = k + 1; i < n; i++) {
            double mu = mat_copy[i][k] / diag_val;
            for (size_t j = 0; j < n; j++)
                mat_copy[i][j] -= mat_copy[k][j] * mu;

            vec_copy[i] -= vec_copy[k] * mu;
        }

        for (size_t i = 0; i < k; i++) {
            double mu = mat_copy[i][k] / diag_val;
            for (size_t j = 0; j < n; j++) {
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
    size_t n = row_count();

    Matrix<T> matr = matrix();
    Vector<T> vect = vector();

    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        T sum = 0;
        if (method == "si") {
            for (size_t j = 0; j < n; j++) {
                if (j != i) {
                    sum += matr[i][j] * x[j];
                }
            }
        } else if (method == "seidel") {
            for (size_t j = 0; j < i; j++) {
                sum += matr[i][j] * x[j];
            }
            for (size_t j = i + 1; j < n; j++) {
                sum += matr[i][j] * x[j];
            }
        }

        x[i] = (vect[i] - sum) / matr[i][i];
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
    size_t n = row_count();
    size_t m = col_count();

    Matrix<T> result(n, m);

    std::ifstream file(file_path);
    std::string line;

    for (size_t i = 0; i < n; i++) {
        std::getline(file, line);
        std::istringstream stream(line);
        T value;

        Vector<T> row(m);

        for (size_t j = 0; j < m; j++) {
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
    size_t n = row_count();
    size_t m = col_count();

    Vector<T> result(n);

    std::ifstream file(file_path);
    std::string line;

    for (size_t i = 0; i < n; i++) {
        std::getline(file, line);
        std::istringstream stream(line);
        T value;
    
        for (size_t j = 0; j < m; j++) {
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

    Matrix<T> matr = matrix();
    Vector<T> vect = vector();

    Matrix<T> result(n, m + 1);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            result[i][j] = matr[i][j];
        }
        result[i][m] = vect[i];
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
    // TODO : review and remake it (done?)

    return is_compatible() && is_diag_d();
}

template<typename T>
void SoLE<T>::to_diag_d(

) {
    elimination();

    // size_t n = row_count();

    // Matrix<T> E(n, n);
    // for (size_t i = 0; i < n; i++) {
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
        }

        x = x_next;
    }

    return x;
}

template<typename T>
Vector<T> SoLE<T>::solve(
    std::string method
) const {
    size_t n = row_count();

    Matrix<T> matr = matrix();
    Vector<T> vect = vector();

    Vector<T> x(n);

    if (method == "thomas") {
        Vector<T> P(n);
        Vector<T> Q(n);

        P[0] = -matr[0][1] / matr[0][0];
        Q[0] = vect[0] / matr[0][0];

        for (size_t i = 1; i < n; i++) {
            P[i] = -matr[i][i + 1] / (matr[i][i] + matr[i][i - 1] * P[i - 1]);
            Q[i] = (vect[i] - matr[i][i - 1] * Q[i - 1]) / (matr[i][i] + matr[i][i - 1] * P[i - 1]);
        }

        x[n - 1] = Q[n - 1];

        for (size_t i = n - 1; i > 0; i--) {
            x[i - 1] = P[i - 1] * x[i] + Q[i - 1];
        }
    } else if (method == "lu") {
        Matrix<T> L(n, n);
        Matrix<T> U(n, n);

        Vector<T> y(n);

        for (int k = 0; k < n; k++) {
            U[k][k] = matr[k][k];

            for (int i = k + 1; i < n; i++) {
                L[i][k] = matr[i][k] / U[k][k];
                U[k][i] = matr[k][i];
            }
            
            for (int i = k + 1; i < n; i++) {
                for(int j = k + 1; j < n; j++) {
                    matr[i][j] = matr[i][j] - L[i][k] * U[k][j];
                }
            }
        }

        for (int i = 0; i < n; i++) {
            y[i] = vect[i];

            for (int j = 0; j < i; j++) {
                y[i] -= L[i][j] * y[j];
            }
        }

        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];

            for (int j = i + 1; j < n; j++) {
                x[i] -= U[i][j] * x[j];
            }

            x[i] /= U[i][i];
        }
    }

    return x;
}
