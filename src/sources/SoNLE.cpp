#include <string>
#include <iostream>
#include <functional>
#include <cmath>

#include <omp.h>

#include "../headers/Vector.h"
#include "../headers/Matrix.h"
#include "../headers/SoLE.h"
#include "../headers/SoNLE.h"


//template class SoNLE<int>;
template class SoNLE<double>;


template<typename T>
Vector<double> SoNLE<T>::calc_next_x(
    Vector<double> x,
    Vector<double> x_prev,
    std::string method
) const {
    size_t n = row_count();

    int mtd = 0;
    
    if (method == "si") {
        mtd = 1;
    } else if (method == "newton") {
        mtd = 2;
    } else if (method == "sign") {
        mtd = 3;
    }

    const double epsilon = 1e-10;

    switch (mtd) {
        case 1: {
            x.set_vector(
                vals(x).vec()
            );

            break;
        }
        case 2: {
            x.set_vector(
                (x - jacobian(x).inversed() * vals(x)).vec()
            );

            break;
        }
        case 3: {
            Vector<double> residual = vals(x);

            Vector<double> sign(n);
            for (size_t i = 0; i < n; i++) {
                sign[i] = (residual[i] >= 0) ? 1.0 : -1.0;
            }

            Vector<double> delta_x = jacobian(x).inversed() * sign;

            double alpha = 1.0;
            double beta = 0.5;

            while (alpha > epsilon) {
                Vector<double> x_new = x - delta_x * alpha;
                Vector<double> residual_new = vals(x_new);

                if (residual_new.norm() < residual.norm()) {
                    x.set_vector(x_new.vec());
                    break;
                } else {
                    alpha *= beta;
                }
            }

            break;
        }
    }

    return x;
}

template<typename T>
size_t SoNLE<T>::row_count(

) const {
    return _row_cnt;
}

template<typename T>
std::function<double(T*)>* SoNLE<T>::funcs(

) const {
    return _funcs;
}

template<typename T>
void SoNLE<T>::set_row_count(
    size_t row_count
) {
    _row_cnt = row_count;
}

template<typename T>
void SoNLE<T>::set_funcs(
    std::function<double(T*)>* funcs
) {
    for (size_t i = 0; i < _row_cnt; i++) {
        _funcs[i] = funcs[i];
    }
}

template<typename T>
Vector<double> SoNLE<T>::vals(
    Vector<T> x
) const {
    size_t n = row_count();

    Vector<double> result(n);
    for (size_t i = 0; i < n; i++) {
        result[i] = _funcs[i](x.vec());
    }

    return result;
}

int factorial(int n) {
    return n <= 1 ? 1 : n * factorial(n - 1);
}

template<typename T>
double SoNLE<T>::diff(
    Vector<double> x,
    int order,
    size_t func_index,
    size_t var_index,
    double h
) const {
    size_t n = row_count();

    double result = 0.0;

    if (order % 2 == 0) {
        #pragma omp parallel for
        for (int k = -order / 2; k <= order / 2; ++k) {
            Vector<double> x_diff(n);

            x_diff.set_vector(x.vec());
            x_diff[var_index] += k * h;

            double coeff = (k == 0) ? 1.0 : (k % 2 == 0 ? 1.0 : -1.0) / (factorial(2 * std::abs(k)));

            result += coeff * vals(x_diff)[func_index];
        }
    } else {
        #pragma omp parallel for
        for (int k = -order; k <= order; ++k) {
            if (k == 0) {
                continue;
            }

            Vector<double> x_diff(n);

            x_diff.set_vector(x.vec());
            x_diff[var_index] += k * h;

            result += (k % 2 == 0 ? -1.0 : 1.0) * vals(x_diff)[func_index] / (k * factorial(order));
        }
    }

    return result / (h * std::pow(2, order));
}


template<typename T>
Matrix<double> SoNLE<T>::jacobian(
    Vector<double> x,
    double h
) const {
    size_t n = row_count();

    Matrix<double> result(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            result[i][j] = diff(x, 1, i, j, h);
        }
    }

    return result;
}

template<typename T>
Vector<double> SoNLE<T>::solve_iter(
    Vector<double> x0,
    std::string method,
    int iter_max,
    double epsilon
) const {
    size_t n = row_count();

    Vector<double> x(n);
    Vector<double> x_prev(n);

    x.set_vector(x0.vec());
    x_prev = x - (x0 / 100);

    for (int i = 0; i < iter_max; i++) {
        Vector<double> x_next = calc_next_x(x, x_prev, method);

        double x_norm = 0.0;
        for (size_t j = 0; j < n; j++) {
            double diff_abs = std::abs(x[j] - x_next[j]);

            if (diff_abs > x_norm) {
                x_norm = diff_abs;
            }
        }

        // std::cout << "i: " << i + 1 << ", norm: " << x_norm << std::endl;

        if (x_norm < epsilon) {
            break;
        }

        x_prev = x;
        x = x_next;
    }

    return x;
}
