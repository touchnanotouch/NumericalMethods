#include <string>
#include <functional>
#include <cmath>

#include <omp.h>

#include "../headers/NLE.h"


template class NLE<int>;
template class NLE<double>;


template<typename T>
double NLE<T>::calc_next_x(
    T x,
    T x_prev,
    double a,
    double b,
    std::string method
) const {
    int mtd = 0;

    if (method == "si") {
        mtd = 1;
    } else if (method == "newton") {
        mtd = 2;
    } else if (method == "secant") {
        mtd = 3;
    } else if (method == "chord") {
        mtd = 4;
    } else if (method == "bisec") {
        mtd = 5;
    } else if (method == "muller") {
        mtd = 6;
    } else if (method == "brent") {
        mtd = 7;
    }

    switch (mtd) {
        case 1: {
            double h = 1e-1;

            double f_x = func()(x);
            double f_x_h = func()(x + h);
            double df_x = (f_x_h - f_x) / h;

            double h_min = 1e-6;
            double h_max = 1e-2;
            double epsilon = 1e-5;

            if (std::abs(f_x_h - f_x) < epsilon) {
                h = std::min(h_max, std::max(h_min, std::abs(f_x) / std::abs(df_x)));
            } else {
                h /= 2;
            }

            x = x - f_x / df_x;

            break;
        }
        case 2: {
            x = x - func()(x) / diff(x, 1);

            break;
        }
        case 3: {
            x = x - func()(x) / ((func()(x) - func()(x_prev)) / (x - x_prev));

            break;
        }
        case 5: {
            double fa = func()(a);
            double fb = func()(b);

            if (fa * fb > 0) {
                return x;
            }

            double c = (a + b) / 2;
            double fc = func()(c);

            if (fc == 0) {
                return c;
            }

            if (fa * fc < 0) {
                b = c;
            } else {
                a = c;
            }

            x = (a + b) / 2;

            break;
        }
    }

    return x;
}

template<typename T>
std::function<T(T)> NLE<T>::func(

) const {
    return _func;
}

template<typename T>
void NLE<T>::set_func(
    std::function<T(T)> func
) {
    _func = func;
}

template<typename T>
double NLE<T>::diff(
    T x,
    int n,
    double h
) const {
    if (n == 0) {
        return func()(x);
    } else {
        return (diff(x + h, n - 1, h) - diff(x, n - 1, h)) / h;
    }
}

template<typename T>
double NLE<T>::solve_iter(
    T x0,
    double a,
    double b,
    std::string method,
    int iter_max,
    double epsilon
) {
    T x = x0;
    T x_prev = x + ((b - a) / 2 * 1e-10);

    std::cout << "\n\n" << method << std::endl;

    for (int i = 0; i < iter_max; i++) {
        T x_next = calc_next_x(x, x_prev, a, b, method);

         std::cout << "i: " << i + 1 << ", norm: " << std::abs(x - x_next) << std::endl;

        if (std::abs(x - x_next) < epsilon) {
            break;
        }

        x_prev = x;
        x = x_next;
    }

    return x;
}
