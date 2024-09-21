#include "../headers/Matrix.h"


template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;


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
    T* temp = a;
    a = b;
    b = temp;
}

template<typename T>
void Matrix<T>::swap(
    T a,
    T b
) {
    T temp = a;
    a = b;
    b = temp;
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
            T* vec = vector.vec();
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
    Matrix result(col_count(), row_count());
    for (int i = 0; i < row_count(); i++) {
        for (int j = 0; j < col_count(); j++) {
            result[j][i] = mat()[i][j];
        }
    }
    
    *this = result;
}

template<typename T>
Matrix<T> Matrix<T>::transposed(
    
) const {
    Matrix result(col_count(), row_count());
    for (int i = 0; i < row_count(); i++) {
        for (int j = 0; j < col_count(); j++) {
            result[j][i] = mat()[i][j];
        }
    }
    
    return result;
}

//template<typename T>
//void Matrix<T>::inverse(
//
//) {
//    ;
//}
//
//template<typename T>
//Matrix<T> Matrix<T>::inversed(
//
//) {
//    return *this;
//}

template<typename T>
double Matrix<T>::determinant(

) const {
    if (row_count() != col_count()) {
        throw std::invalid_argument("Unmatching matrix size");
    }

    // TODO

    int n = row_count();

    double det = 1;

    return det;
}

template<typename T>
int Matrix<T>::rank(

) const {
    if (row_count() != col_count()) {
        throw std::invalid_argument("Unmatching matrix size");
    }

    // TODO

    int n = (int)row_count();

    int rank = n - 1;

    return rank;
}

template<typename T>
void Matrix<T>::to_up_diag(

) {
    int n = (int)row_count();

    for (int r = 0; r < n; r++) {
        T diag = mat()[r][r];
        for (int c = r + 1; c < n; c++) {
            
        }
    }
}

template<typename T>
bool Matrix<T>::is_diag_d(

) const {
    for (int i = 0; i < row_count(); i++) {
        T diag = std::abs(mat()[i][i]);
        
        T sum = 0;
        for (size_t j = 0; j < col_count(); j++) {
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
