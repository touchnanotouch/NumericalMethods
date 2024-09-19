#include "../headers/Matrix.h"
#include "../headers/Vector.h"


template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;

template<typename T>
Matrix<T> Matrix<T>::gaussian_elimination(
    Matrix<T>& matrix
) {
    Matrix result(matrix._row_cnt, matrix._col_cnt);
    for (size_t i = 0; i < matrix._row_cnt; i++) {
        for (size_t j = 0; j < matrix._col_cnt; j++) {
            result._mat[i][j] = matrix._mat[i][j];
        }
    }

    for (size_t r = 0; r < result._row_cnt; r++) {
        T diag_val = result._mat[r][r];

        for (size_t i = r + 1; i < result._row_cnt; i++) {
            T factor = result._mat[i][r] / diag_val;
            for (size_t j = r; j < result._col_cnt; j++) {
                result._mat[i][j] -= result._mat[r][j] * factor;
            }
        }

        for (size_t j = r + 1; j < result._col_cnt; j++) {
            T factor = result._mat[j][r] / diag_val;
            for (size_t k = r; k < result._col_cnt; k++) {
                result._mat[j][k] -= result._mat[r][k] * factor;
            }
        }
    }

    return result;
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
            _mat[i][j] = other._mat[i][j];
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
                _mat[i][j] = other._mat[i][j];
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
size_t Matrix<T>::get_row_count(

) const { 
    return _row_cnt;
}

template<typename T>
size_t Matrix<T>::get_col_count(

) const {
    return _col_cnt;
}

template<typename T>
T** Matrix<T>::get_matrix(

) const {
    return _mat;
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
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            if (_mat[i][j] != other._mat[i][j] ) {
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
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            result._mat[i][j] = _mat[i][j] + other._mat[i][j];
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
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            result._mat[i][j] = _mat[i][j] - other._mat[i][j];
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
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < other._col_cnt; j++) {
            for (int k = 0; k < _col_cnt; k++) {
                result._mat[i][j] += _mat[i][k] * other._mat[k][j];
            }
        }
    }

    return result;
}

template<typename T>
Vector<T> Matrix<T>::operator*(
    const Vector<T>& vector
) const {
    if (_col_cnt != vector.get_row_count()) {
        throw std::invalid_argument("Unmatching matrix and vector sizes");
    }

    Vector<T> result(_row_cnt);
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            T* vec = vector.get_vector();
            result[i] += _mat[i][j] * vec[j];
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(
    const T& scalar
) const {
    Matrix result(_row_cnt, _col_cnt);
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
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
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            result._mat[i][j] = _mat[i][j] / scalar;
        }
    }

    return result;
}

template<typename T>
void Matrix<T>::transpose(

) {
    Matrix result(_col_cnt, _row_cnt);
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            result._mat[j][i] = _mat[i][j];
        }
    }
    
    *this = result;
}

template<typename T>
Matrix<T> Matrix<T>::transposed(
    
) const {
    Matrix result(_col_cnt, _row_cnt);
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            result._mat[j][i] = _mat[i][j];
        }
    }
    
    return result;
}

template<typename T>
double Matrix<T>::determinant(

) {
    if (_row_cnt != _col_cnt) {
        throw std::invalid_argument("Unmatching matrix size");
    }

    Matrix temp = gaussian_elimination(*this);

    double det = 1;
    for (int i = 0; i < temp._row_cnt; i++) {
        det *= temp[i][i];
    }

    return det;
}

template<typename T>
int Matrix<T>::rank(

) {
    Matrix temp = gaussian_elimination(*this);

    int rank = 0;
    for (size_t i = 0; i < temp._row_cnt; i++) {
        bool is_zero_row = true;
        for (size_t j = 0; j < temp._col_cnt; j++) {
            if (std::abs(temp._mat[i][j]) > 1e-10) {
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
