#include "../headers/Vector.h"


template class Vector<int>;
template class Vector<float>;
template class Vector<double>;


template<typename T>
Vector<T>::Vector(
    const Vector<T>& other
) {
    _row_cnt = other._row_cnt;
    _vec = new T[_row_cnt];
    for (size_t i = 0; i < _row_cnt; i++) {
        _vec[i] = other._vec[i];
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
            _vec[i] = other._vec[i];
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
    for (int i = 0; i < _row_cnt; i++) {
        _vec[i] = vector[i];
    }
}

template<typename T>
void Vector<T>::swap(
    T a,
    T b
) {
    T temp = a;
    a = b;
    b = temp;
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
        return false;
    }

    const double eps = 1e-5;
    for (size_t i = 0; i < _row_cnt; i++) {
        if (_vec[i] != other._vec[i]) {
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
        result._vec[i] = _vec[i] + other._vec[i];
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
        result._vec[i] = _vec[i] - other._vec[i];
    }

    return result;
}

template<typename T>
Vector<T> Vector<T>::operator*(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result._vec[i] = _vec[i] * scalar;
    }

    return result;
}


template<typename T>
Vector<T> Vector<T>::operator/(
    const T& scalar
) const {
    Vector result(_row_cnt);
    for (size_t i = 0; i < _row_cnt; i++) {
        result._vec[i] = _vec[i] / scalar;
    }

    return result;
}

template<typename T>
T Vector<T>::norm(

) const {
    T result = 0;
    for (size_t i = 0; i < row_count(); i++) {
        result += std::abs(vec()[i]);
    }

    return result;
}
