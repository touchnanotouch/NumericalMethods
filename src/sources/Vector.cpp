#include <sstream>
#include <fstream>

#include <omp.h>

#include "../headers/Vector.h"


template class Vector<int>;
template class Vector<double>;


template<typename T>
size_t Vector<T>::row_count(
    
) const {
    return _row_cnt;
}

template<typename T>
T* Vector<T>::vec(

) const {
    return _vec;
}

template<typename T>
void Vector<T>::set_row_count(
    size_t row_count
) {
    _row_cnt = row_count;
}

template<typename T>
void Vector<T>::set_vector(
    T* vector
) {
    for (size_t i = 0; i < _row_cnt; i++) {
        _vec[i] = vector[i];
    }
}

template<typename T>
void Vector<T>::set_vector(
    std::string file_path,
    const char delimiter
) {
    size_t n = row_count();

    Vector<T> result(n);

    std::ifstream file(file_path);
    std::string line;

    size_t i = 0;
    while (std::getline(file, line) && i < n) {
        std::istringstream stream(line);

        T value;
        while (stream >> value) {
            result[i++] = value;

            if (stream.peek() == delimiter) {
                stream.ignore();
            }
        }
    }

    set_vector(result.vec());

    file.close();
}

template<typename T>
void Vector<T>::swap(
    T& a,
    T& b
) {
    T temp = a;
    a = b;
    b = temp;
}

template<typename T>
Vector<T>::Vector(
    const Vector<T>& other
) {
    _row_cnt = other._row_cnt;
    _vec = new T[_row_cnt];
    for (size_t i = 0; i < _row_cnt; i++) {
        _vec[i] = other[i];
    }
}

template<typename T>
Vector<T>& Vector<T>::operator=(
    const Vector<T>& other
) {
    if (this != &other) {
        delete[] _vec;

        _row_cnt = other._row_cnt;
        _vec = new T[_row_cnt];
        for (size_t i = 0; i < _row_cnt; i++) {
            _vec[i] = other[i];
        }
    }

    return *this;
}

template<typename T>
Vector<T>::Vector(
    Vector<T>&& other
) : _row_cnt(other._row_cnt), _vec(other._vec) {
    other._vec = nullptr;
}

template<typename T>
Vector<T>& Vector<T>::operator=(
    Vector<T>&& other
) {
    if (this != &other) {
        _row_cnt = other._row_cnt;
        _vec = other._vec;
        other._vec = nullptr;
    }

    return *this;
}

template<typename T>
T& Vector<T>::operator[](
    size_t index
) const {
    if (index < 0 || index >= _row_cnt) {
        throw std::out_of_range("Index out of range");
    }

    return _vec[index];
}

template<typename T>
bool Vector<T>::operator==(
    const Vector& other
) const {
    if (_row_cnt != other._row_cnt) {
        throw std::invalid_argument("Unmatching vector sizes");
    }

    const double eps = 1e-2;
    for (size_t i = 0; i < _row_cnt; i++) {
        if (_vec[i] != other[i] && std::abs(_vec[i] - other[i]) > eps) {
            return false;
        }
    }

    return true;
}

template<typename T>
Vector<T> Vector<T>::operator+(
    const Vector& other
) const {
    if (_row_cnt != other._row_cnt) {
        throw std::invalid_argument("Unmatching vector sizes");
    }

    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] + other[i];
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator+(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] + scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator+=(
    const Vector& other
) const {
    if (_row_cnt != other._row_cnt) {
        throw std::invalid_argument("Unmatching vector sizes");
    }

    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] + other[i];
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator+=(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] + scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator-(
    const Vector& other
) const {
    if (_row_cnt != other._row_cnt) {
        throw std::invalid_argument("Unmatching vector sizes");
    }

    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] - other[i];
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator-(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] - scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator-=(
    const Vector& other
) const {
    if (_row_cnt != other._row_cnt) {
        throw std::invalid_argument("Unmatching vector sizes");
    }

    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] - other[i];
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator-=(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] - scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator*(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] * scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator*=(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] * scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator/(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] / scalar;
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator/=(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result[i] = _vec[i] / scalar;
    }

    return result;
}

template<typename T>
T Vector<T>::dot(
    const Vector<T>& other
) const {
    size_t n = row_count();

    T result = 0;
    #pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < n; i++) {
        result += _vec[i] * other._vec[i];
    }

    return result;
}

template<typename T>
double Vector<T>::norm(

) const {
    size_t n = row_count();

    T* vec = this->vec();

    double result = 0;
    #pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < n; i++) {
        result += std::abs(vec[i]);
    }

    return result;
}
