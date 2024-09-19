#include "../headers/Matrix.h"
#include "../headers/Vector.h"
#include "../headers/SoLE.h"


template class SoLE<int>;
template class SoLE<float>;
template class SoLE<double>;

template<typename T>
Matrix<T> SoLE<T>::get_matrix(

) const {
    return _coeffs;
}

template<typename T>
Vector<T> SoLE<T>::get_vector(
    
) const {
    return _terms;
}

template<typename T>
void SoLE<T>::set_matrix(
    Matrix<T>& matrix
) {
    _coeffs = matrix;
}

template <typename T>
void SoLE<T>::set_vector(
    Vector<T>& vector
) {
    _terms = vector;
}

template<typename T>
bool SoLE<T>::is_solution(
    Vector<T>& solution
) const {
    return _coeffs * solution == _terms;
}
