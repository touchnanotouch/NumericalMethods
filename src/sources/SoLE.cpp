#include <string>
#include <cmath>

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

        #pragma omp parallel for
        for (size_t i = k + 1; i < n; i++) {
            double mu = mat_copy[i][k] / diag_val;

            for (size_t j = 0; j < n; j++)
                mat_copy[i][j] -= mat_copy[k][j] * static_cast<T>(mu);

            vec_copy[i] -= vec_copy[k] * static_cast<T>(mu);
        }

        #pragma omp parallel for
        for (size_t i = 0; i < k; i++) {
            double mu = mat_copy[i][k] / diag_val;

            for (size_t j = 0; j < n; j++) {
                mat_copy[i][j] -= mat_copy[k][j] * static_cast<T>(mu);
            }

            vec_copy[i] -= vec_copy[k] * static_cast<T>(mu);
        }
    }

    #pragma omp critical
    {
        set_matrix(mat_copy);
        set_vector(vec_copy);
    }
}

template<typename T>
Vector<T> SoLE<T>::calc_next_x(
    Vector<T> x,
    std::string method,
    double res_omega,
    int grad_iter
) const {
    size_t n = row_count();

    Matrix<T> matr = matrix();
    Vector<T> vect = vector();

    int mtd = 0;
    
    if (method == "si") {
        mtd = 1;
    } else if (method == "seidel") {
        mtd = 2;
    } else if (method == "sor") {
        mtd = 3;
    } else if (method == "res") {
        mtd = 4;
    } else if (method == "grad") {
        mtd = 5;
    }

    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        T sum = 0;
 
        switch (mtd) {
            case 1: {
                #pragma omp parallel for
                for (size_t j = 0; j < n; j++) {
                    if (j != i) {
                        sum += matr[i][j] * x[j];
                    }
                }

                #pragma omp critical
                {
                    x[i] = (vect[i] - sum) / matr[i][i];
                }

                break;
            }
            case 2: {
                #pragma omp parallel for
                for (size_t j = 0; j < n; j++) {
                    if (j != i) {
                        sum += matr[i][j] * x[j];
                    }
                }

                #pragma omp critical
                {
                    x[i] = (vect[i] - sum) / matr[i][i];
                }

                break;
            }
            case 3: {
                #pragma omp parallel for
                for (size_t j = 0; j < n; j++) {
                    if (j != i) {
                        sum += matr[i][j] * x[j];
                    }
                }

                #pragma omp critical
                {
                    x[i] = (1 - res_omega) * x[i] + res_omega * (vect[i] - sum) / matr[i][i];   
                }

                break;
            }
            case 4: {
                Vector<T> r = vect - matr * x;
                Vector<T> matr_r = matr * r;
        
                x[i] += (r.dot(r)) / matr_r.dot(matr_r) * r[i];
                break;
            }
            case 5: {
                double cf = 0.9;

                Vector<T> grad = vect - matr * x;
                
                double alpha = 1 / (1 + grad.norm());
                double new_alpha = alpha;

                #pragma omp parallel for
                for (size_t k = 0; k < n; k++) {
                    Vector<T> new_x(n);
                    new_x[k] = x[k] + new_alpha * grad[k];

                    Vector<T> r = vect - matr * new_x;
                    Vector<T> matr_r = matr * r;

                    double r_dot_r = r.dot(r);
                    double matr_r_dot_matr_r = matr_r.dot(matr_r);

                    if (r_dot_r <= matr_r_dot_matr_r * new_alpha) {
                        break;
                    }

                    new_alpha *= cf;
                }

                #pragma omp parallel for
                for (size_t k = 0; k < n; k++) {
                    x[k] = x[k] + new_alpha * grad[k];
                }

                break;
            }
        }
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
    _coeffs.set_matrix(file_path, delimiter);
}

template <typename T>
void SoLE<T>::set_vector(
    Vector<T>& vector
) {
    _terms = vector;
}

template<typename T>
void SoLE<T>::set_vector(
    std::string file_path
) {
    _terms.set_vector(file_path);
}

template<typename T>
Matrix<T> SoLE<T>::extended(
    
) const {
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

) const {
    return matrix().det();
}

template<typename T>
int SoLE<T>::rank(

) const {
    return matrix().rank();
}

template<typename T>
bool SoLE<T>::is_compatible(

) const {
    return extended().rank() == matrix().rank();
}

template<typename T>
bool SoLE<T>::is_diag(

) const {
    return matrix().is_diag();
}

template<typename T>
bool SoLE<T>::is_diag_d(

) const {
    return matrix().is_diag_d();
}

template<typename T>
bool SoLE<T>::is_solvable(
    
) const {
    // TODO : review and remake it (need new concept)

    return is_compatible() && is_diag_d();
}

template<typename T>
bool SoLE<T>::is_solution(
    Vector<T>& solution
) const {
    return matrix() * solution == vector();
}

template<typename T>
void SoLE<T>::to_diag(
    
) {
    elimination();
}

template<typename T>
void SoLE<T>::to_diag_d(

) {
    // to_diag();

    size_t n = row_count();

    Matrix<T> E(n, n);
    for (size_t i = 0; i < n; i++) {
        E[i][i] = (T)1;
    }

    Matrix<T> mat_new = matrix().inversed() * E;
    
    set_matrix(mat_new);
}

template<typename T>
Vector<T> SoLE<T>::solve_iter(
    std::string method,
    int iter_max,
    double epsilon,
    double res_omega,
    int grad_iter
) const {
    Vector<T> x(row_count());

    for (int i = 0; i < iter_max; i++) {
        Vector<T> x_next = calc_next_x(x, method, res_omega, grad_iter);

        double x_norm = (x_next - x).norm();

        // std::cout << "iter " << i + 1 << ", norm: " << x_norm << std::endl;

        if (x_norm < epsilon) {
            // std::cout << std::endl;
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

    if (method == "gauss") {
        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            size_t max_el = i;

            #pragma omp parallel for reduction(max:max_el)
            for (size_t k = i + 1; k < n; k++) {
                if (abs(matr[k][i]) > abs(matr[max_el][i])) {
                    max_el = k;
                }
            }

            if (max_el != i) {
                matr.swap(matr[i], matr[max_el]);
                vect.swap(vect[i], vect[max_el]);
            }

            #pragma omp parallel for
            for (size_t k = i + 1; k < n; k++) {
                T c = matr[k][i] / matr[i][i];
                T local_vect = vect[k];

                #pragma omp parallel for
                for (size_t j = i; j < n; j++) {
                    local_vect -= c * matr[i][j];
                }

                vect[k] = local_vect;
            }
        }

        #pragma omp parallel for
        for (size_t i = n - 1; i > 0; i--) {
            x[i] = vect[i];

            #pragma omp parallel for reduction(-:x[i])
            for (size_t k = i + 1; k < n - 1; k++) {
                x[i] -= matr[i][k] * x[k];
            }

            x[i] /= matr[i][i];
        }
    } else if (method == "thomas") {
        Vector<T> P(n);
        Vector<T> Q(n);

        P[0] = -matr[0][1] / matr[0][0];
        Q[0] = vect[0] / matr[0][0];

        #pragma omp parallel for
        for (size_t i = 1; i < n; i++) {
            #pragma omp critical
            {
                P[i] = -matr[i][i + 1] / (matr[i][i] + matr[i][i - 1] * P[i - 1]);
                Q[i] = (vect[i] - matr[i][i - 1] * Q[i - 1]) / (matr[i][i] + matr[i][i - 1] * P[i - 1]);
            }
        }

        x[n - 1] = Q[n - 1];

        #pragma omp parallel for
        for (size_t i = n - 1; i > 0; i--) {
            #pragma omp critical
            {
                x[i - 1] = P[i - 1] * x[i] + Q[i - 1];
            }
        }
    } else if (method == "lu") {
        Matrix<T> L(n, n);
        Matrix<T> U(n, n);

        Vector<T> y(n);

        #pragma omp parallel for
        for (size_t k = 0; k < n; k++) {
            U[k][k] = matr[k][k];

            #pragma omp parallel for
            for (size_t i = k + 1; i < n; i++) {
                #pragma omp critical
                {
                    L[i][k] = matr[i][k] / U[k][k];
                    U[k][i] = matr[k][i];
                }
            }
            
            #pragma omp parallel for
            for (size_t i = k + 1; i < n; i++) {
                #pragma omp parallel for
                for (size_t j = k + 1; j < n; j++) {
                    #pragma omp critical
                    {
                        matr[i][j] = matr[i][j] - L[i][k] * U[k][j];
                    }
                }
            }
        }

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            y[i] = vect[i];

            #pragma omp parallel for
            for (size_t j = 0; j < i; j++) {
                #pragma omp critical
                {
                    y[i] -= L[i][j] * y[j];
                }
            }
        }

        #pragma omp parallel for
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];

            #pragma omp parallel for
            for (size_t j = i + 1; j < n; j++) {
                #pragma omp critical
                {
                    x[i] -= U[i][j] * x[j];
                }
            }

            x[i] /= U[i][i];
        }
    } else if (method == "qroots") {
        Matrix<T> L(n, n);

        Vector<T> y(n);

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            #pragma omp parallel for
            for (size_t j = 0; j <= i; j++) {
                T sum = 0;

                #pragma omp parallel for
                for (size_t k = 0; k < j; k++) {
                    #pragma omp critical
                    {
                        sum += L[i][k] * L[j][k];
                    }
                }

                if (i == j) {
                    #pragma omp critical
                    {
                        L[i][j] = std::sqrt(matr[i][i] - sum);
                    }
                } else {
                    #pragma omp critical
                    {
                        L[i][j] = (matr[i][j] - sum) / L[j][j];
                    }
                }
            }
        }

        Matrix<T> U = L.transposed();

        #pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            y[i] = vect[i];

            #pragma omp parallel for
            for (size_t j = 0; j < i; j++) {
                y[i] -= L[i][j] * y[j];
            }

            y[i] /= L[i][i];
        }

        #pragma omp parallel for
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];

            #pragma omp parallel for
            for (size_t j = i + 1; j < n; j++) {
                x[i] -= U[i][j] * x[j];
            }

            x[i] /= U[i][i];
        }
    }

    return x;
}
