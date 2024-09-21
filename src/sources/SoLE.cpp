#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>

#include "../headers/SoLE.h"


template class SoLE<int>;
template class SoLE<float>;
template class SoLE<double>;

template<typename T>
Vector<T> SoLE<T>::calc_next_x(
    Vector<T> x,
    std::string method
) const {
    for (int i = 0; i < row_count(); i++) {
        T sum = 0;
        if (method == "si") {
            for (int j = 0; j < row_count(); j++) {
                if (j != i) {
                    sum += matrix()[i][j] * x[j];
                }
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
double SoLE<T>::determinant(

) {
    return matrix().determinant();
}

template<typename T>
int SoLE<T>::rank(

) {
    return matrix().rank();
}

template<typename T>
bool SoLE<T>::is_solution(
    Vector<T>& solution
) const {
    return matrix() * solution == vector();
}

template<typename T>
Vector<T> SoLE<T>::solve_si_method(
    int iter_max,
    double epsilon
) const {
    Vector<T> x(row_count());

    for (int i = 0; i < iter_max; i++) {
        Vector<T> x_next = calc_next_x(x, "si");

        if ((x_next - x).norm() < epsilon) {
            break;
        }

        x = x_next;
    }

    return x;
}

template<typename T>
bool SoLE<T>::to_diag_d(

) {
    int min_s = -10;
    int max_s = 10;
    int size = max_s - min_s + 1;

    size_t n = row_count();

    int* indices = new int[n];
    for (int i = 0; i < n; i++) {
        indices[i] = i;
    }

    do {
        Matrix<T> mat_new(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                mat_new[i][j] = matrix()[indices[i]][j];
            }
        }

        if (mat_new.is_diag_d()) {
            Vector<T> vector_new(n);
            for (int l = 0; l < n; l++) {
                vector_new[l] = vector()[indices[l]];
            }

            set_matrix(mat_new);
            set_vector(vector_new);

            delete[] indices;

            return true;
        }
    } while (std::next_permutation(indices, indices + n));



    delete[] indices;

    return false;
}

