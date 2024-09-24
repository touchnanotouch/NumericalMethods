#include "../headers/Matrix.h"


template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;


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
    for (int i = 0; i < _row_cnt; i++) {
        _mat[index][i] = row[i];
    }
}

template<typename T>
void Matrix<T>::set_col(
    T* col,
    int index
) {
    for (int i = 0; i < _col_cnt; i++) {
        _mat[i][index] = col[i];
    }
}

template<typename T>
void Matrix<T>::set_matrix(
    T** matrix
) {
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
            _mat[i][j] = matrix[i][j];
        }
    }
}

template<typename T>
void Matrix<T>::swap(
    T* a,
    T* b
) {
    for (int i = 0; i < row_count(); i++) {
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
    if (_col_cnt != vector.row_count()) {
        throw std::invalid_argument("Unmatching matrix and vector sizes");
    }

    Vector<T> result(_row_cnt);
    for (int i = 0; i < _row_cnt; i++) {
        for (int j = 0; j < _col_cnt; j++) {
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
Matrix<T> Matrix<T>::transposed(
    
) const {
    int n = (int)row_count();
    int m = (int)col_count();

    Matrix result(n, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[j][i] = mat()[i][j];
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
    int n = (int)row_count();

    Matrix<T> cof(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix<T> minor(n - 1, n - 1);
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    if (k != i && l != j) {
                        minor[k < i ? k : k - 1][l < j ? l : l - 1] = mat()[k][l];
                    }
                }
            }

            cof[i][j] = pow(-1, i + j) * minor.det();
        }
    }

    return cof.transposed();
}

template<typename T>
Matrix<T> Matrix<T>::inversed(

) {
    int n = (int)row_count();

    double determ = det();

    if (determ == 0) {
        throw std::runtime_error("Матрица не обратима");
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
    if (row_count() == 1) {
        return mat()[0][0];
    } else if (row_count() == 2) {
        return mat()[0][0] * mat()[1][1] - mat()[0][1] * mat()[1][0];
    } else {
        double det = 0;
        for (int i = 0; i < row_count(); i++) {
            Matrix<T> minor(row_count() - 1, col_count() - 1);
            for (int j = 1; j < row_count(); j++) {
                for (int k = 0; k < col_count(); k++) {
                    if (k < i) {
                        minor[j - 1][k] = mat()[j][k];
                    } else if (k > i) {
                        minor[j - 1][k - 1] = mat()[j][k];
                    }
                }
            }
            det += (i % 2 == 0 ? 1 : -1) * mat()[0][i] * minor.det();
        }

        return det;
    }
}

template<typename T>
int Matrix<T>::rank(

) const {
    if (row_count() != col_count()) {
        throw std::invalid_argument("Unmatching matrix size");
    }

    int n = (int)row_count();

    Matrix<T> copy(n, n);
    copy.set_matrix(this->mat());

    int rank = n - 1;

    for (int row = 0; row < rank; row++) {
        if (copy[row][row]) {
            for (int col = 0; col < n; col++) {
                if (col != row) {
                    T mult = copy[col][row] / copy[row][row];
                    for (int i = 0; i < rank; i++)
                    copy[col][i] -= mult * copy[row][i];
                }
            }
        } else {
            bool reduce = true;

            for (int i = row + 1; i < n; i++) {
                if (copy[i][row]) {
                    copy.swap(copy[i], copy[rank]);
                    reduce = false;
                    break ;
                }
            }

            if (reduce) {
                rank--;

                for (int i = 0; i < n; i ++)
                    copy[i][row] = copy[i][rank];
            }

            row--;
        }
    }

    return rank;
}

template<typename T>
bool Matrix<T>::is_diag_d(

) const {
    for (int i = 0; i < (int)row_count(); i++) {
        T diag = std::abs(mat()[i][i]);
        
        T sum = 0;
        for (int j = 0; j < (int)col_count(); j++) {
            if (i != j) {
                sum += std::abs(mat()[i][j]);
            }
        }

        if (sum > diag) {
            return false;
        }
    }

    return true;
}

template<typename T>
bool Matrix<T>::to_diag_d(
    
) {
    return false;
}
