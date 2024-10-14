#include <cmath>

#include <omp.h>

#include "../headers/Matrix.h"


template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;


template<typename T>
void Matrix<T>::elimination(

) {
    size_t n = row_count();

    T** mat = this->mat();

    #pragma omp parallel for
    for (size_t k = 0; k < n; k++) {
        T diag_val = mat[k][k];

        for (size_t i = k + 1; i < n; i++) {
            double mu = mat[i][k] / diag_val;

            for (size_t j = 0; j < n; j++)
                mat[i][j] -= mat[k][j] * mu;
        }

        for (size_t i = 0; i < k; i++) {
            double mu = mat[i][k] / diag_val;

            for (size_t j = 0; j < n; j++) {
                mat[i][j] -= mat[k][j] * mu;
            }
        }
    }
}

template<typename T>
size_t Matrix<T>::row_count(

) const { 
    return _row_cnt;
}

template<typename T>
size_t Matrix<T>::col_count(

) const {
    return _col_cnt;
}

template<typename T>
T** Matrix<T>::mat(

) const {
    return _mat;
}

template<typename T>
void Matrix<T>::set_row_count(
    size_t row_count
) {
    _row_cnt = row_count;
}

template<typename T>
void Matrix<T>::set_col_count(
    size_t col_count
) {
    _col_cnt = col_count;
}

template<typename T>
void Matrix<T>::set_row(
    T* row,
    int index
) {
    for (size_t i = 0; i < _row_cnt; i++) {
        _mat[index][i] = row[i];
    }
}

template<typename T>
void Matrix<T>::set_col(
    T* col,
    int index
) {
    for (size_t i = 0; i < _col_cnt; i++) {
        _mat[i][index] = col[i];
    }
}

template<typename T>
void Matrix<T>::set_matrix(
    T** matrix
) {
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            _mat[i][j] = matrix[i][j];
        }
    }
}

template<typename T>
void Matrix<T>::swap(
    T* a,
    T* b
) {
    for (size_t i = 0; i < _row_cnt; i++) {
        T temp = a[i];
        a[i] = b[i];
        b[i] = temp;
    }
}

template<typename T>
void Matrix<T>::swap(
    T& a,
    T& b
) {
    T temp = a;
    a = b;
    b = temp;
}

template<typename T>
Matrix<T>::Matrix(
    const Matrix<T>& other
) {
    _row_cnt = other._row_cnt;
    _col_cnt = other._col_cnt;

    _mat = new T*[_row_cnt];
    for (size_t i = 0; i < _row_cnt; i++) {
        _mat[i] = new T[_col_cnt];
        for (size_t j = 0; j < _col_cnt; j++) {
            _mat[i][j] = other[i][j];
        }
    }
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(
    const Matrix<T>& other
) {
    if (this != &other) {
        for (size_t i = 0; i < _row_cnt; i++) {
            delete[] _mat[i];
        }
        delete[] _mat;

        _row_cnt = other._row_cnt;
        _col_cnt = other._col_cnt;
        _mat = new T*[_row_cnt];
        for (size_t i = 0; i < _row_cnt; i++) {
            _mat[i] = new T[_col_cnt];
            for (size_t j = 0; j < _col_cnt; j++) {
                _mat[i][j] = other[i][j];
            }
        }
    }

    return *this;
}

template<typename T>
Matrix<T>::Matrix(
    Matrix<T>&& other
) : _row_cnt(other._row_cnt), _col_cnt(other._col_cnt), _mat(other._mat) {
    other._mat = nullptr;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(
    Matrix<T>&& other
) {
    if (this != &other) {
        _row_cnt = other._row_cnt;
        _col_cnt = other._col_cnt;
        _mat = other._mat;
        other._mat = nullptr;
    }

    return *this;
}

template<typename T>
T*& Matrix<T>::operator[](
    int index
) const {
    if (index < 0 || index >= _row_cnt && index >= _col_cnt) {
        throw std::out_of_range("Index out of range");
    }

    return _mat[index];
}

template<typename T>
bool Matrix<T>::operator==(
    const Matrix<T>& other
) const {
    if (_row_cnt != other._row_cnt || _col_cnt != other._col_cnt) {
        return false;
    }

    const double eps = 1e-5;
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            if (_mat[i][j] != other[i][j] && std::abs(_mat[i][j] - other[i][j]) > eps) {
                return false;
            }
        }
    }

    return true;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(
    const Matrix<T>& other
) const {
    if (_row_cnt != other._row_cnt || _col_cnt != other._col_cnt) {
        throw std::invalid_argument("Unmatching matrix sizes");
    }

    Matrix result(_row_cnt, _col_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result[i][j] = _mat[i][j] + other[i][j];
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(
    const T& scalar
) const {
    Matrix result(_row_cnt, _col_cnt);

    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result[i][j] = _mat[i][j] + scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(
    const Matrix& other
) const {
    if (_row_cnt != other._row_cnt || _col_cnt != other._col_cnt) {
        throw std::invalid_argument("Unmatching matrix sizes");
    }

    Matrix result(_row_cnt, _col_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result[i][j] = _mat[i][j] - other[i][j];
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(
    const T& scalar
) const {
    Matrix result(_row_cnt, _col_cnt);

    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result[i][j] = _mat[i][j] - scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(
    const Matrix& other
) const {
    if (_col_cnt != other._row_cnt) {
        throw std::invalid_argument("Unmatching matrix sizes");
    }

    Matrix result(_row_cnt, other._col_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < other._col_cnt; j++) {
            for (size_t k = 0; k < _col_cnt; k++) {
                result[i][j] += _mat[i][k] * other[k][j];
            }
        }
    }

    return result;
}

template<typename T>
Vector<T> Matrix<T>::operator*(
    const Vector<T>& vector
) const {
    if (_col_cnt != vector.row_count()) {
        throw std::invalid_argument("Unmatching matrix and vector sizes");
    }

    Vector<T> result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result[i] += _mat[i][j] * vector[j];
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(
    const T& scalar
) const {
    Matrix result(_row_cnt, _col_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result._mat[i][j] = _mat[i][j] * scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(
    const T& scalar
) const {
    Matrix result(_row_cnt, _col_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        for (size_t j = 0; j < _col_cnt; j++) {
            result._mat[i][j] = _mat[i][j] / scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::transposed(
    
) const {
    size_t n = row_count();
    size_t m = col_count();

    T** mat = this->mat();

    Matrix result(n, m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            result[j][i] = mat[i][j];
        }
    }
    
    return result;
}

template<typename T>
void Matrix<T>::transpose(

) {
    *this = transposed();
}

template<typename T>
Matrix<T> Matrix<T>::cofactors(

) {
    size_t n = row_count();

    T** mat = this->mat();

    Matrix<T> cof(n, n);

    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            Matrix<T> minor(n - 1, n - 1);

            for (size_t k = 0; k < n; k++) {
                for (size_t l = 0; l < n; l++) {
                    if (k != i && l != j) {
                        minor[k < i ? k : k - 1][l < j ? l : l - 1] = mat[k][l];
                    }
                }
            }

            cof[i][j] = std::pow(-1, i + j) * minor.det();
        }
    }

    return cof.transposed();
}

template<typename T>
Matrix<T> Matrix<T>::inversed(

) {
    double determ = det();

    if (determ == 0) {
        throw std::runtime_error("Singular matrix");
    }

    Matrix<T> cof = cofactors() / determ;

    return cof;
}

template<typename T>
void Matrix<T>::inverse(

) {
    *this = inversed();
}

template<typename T>
double Matrix<T>::det(

) const {
    size_t n = row_count();

    Matrix<T> copy(n, n);
    copy.set_matrix(this->mat());

    copy.elimination();

    double det = 1;
    for (size_t i = 0; i < n; i++) {
        det *= copy[i][i];
    }

    return det;
}

template<typename T>
int Matrix<T>::rank(
  
) const {
    size_t n = row_count();
    size_t m = col_count();

    Matrix<T> copy(n, m);
    copy.set_matrix(this->mat());

    copy.elimination();

    int rank = 0;
    for (size_t i = 0; i < n; i++) {
        bool is_zero_row = true;

        for (size_t j = 0; j < m; j++) {
            if (std::abs(copy[i][j]) > 1e-10) {
                is_zero_row = false;
                break;
            }
        }

        if (!is_zero_row) {
            rank++;
        }
    }

    return rank;
}

template<typename T>
bool Matrix<T>::is_diag_d(

) const {
    size_t n = row_count();
    size_t m = col_count();

    T** mat = this->mat();

    for (size_t i = 0; i < n; i++) {
        T diag = std::abs(mat[i][i]);
        
        T sum = 0;
        for (size_t j = 0; j < m; j++) {
            if (i != j) {
                sum += std::abs(mat[i][j]);
            }
        }

        if (sum > diag) {
            return false;
        }
    }

    return true;
}

template<typename T>
void Matrix<T>::to_diag_d(
    
) {
    elimination();
    
    //size_t n = row_count();

    //Matrix<T> E(n, n);
    //for (size_t i = 0; i < n; i++) {
    //    E[i][i] = 1;
    //}

    //Matrix<T> mat_new = inversed() * E;

    //set_matrix(mat_new.mat());
}
