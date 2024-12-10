#include <string>
#include <iostream>
#include <functional>
#include <cmath>

#include <omp.h>

#include "../headers/NLE.h"


// template class NLE<int>;
template class NLE<double>;


template<typename T>
double NLE<T>::calc_next_x(
    double x,
    double x_prev,
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
    } else if (method == "bisec") {
        mtd = 4;
    } else if (method == "chord") {
        mtd = 5;
    }

    const double epsilon = 1e-8;

    switch (mtd) {
        case 1: {
            double h = 1e-1;

            double f_x = val(x);
            double f_x_h = val(x + h);
            double df_x = (f_x_h - f_x) / h;

            double h_min = 1e-6;
            double h_max = 1e-2;

            if (std::abs(f_x_h - f_x) < epsilon) {
                h = std::min(h_max, std::max(h_min, std::abs(f_x) / std::abs(df_x)));
            } else {
                h /= 2;
            }

            x -= f_x / df_x;

            break;
        }
        case 2: {
            x -= val(x) / diff(x, 1);

            break;
        }
        case 3: {
            x -= (val(x) / (val(x) - val(x_prev))) * (x - x_prev);

            break;
        }
        case 4: {
            double c = (a + b) / 2;

            while (std::abs(val(c)) > epsilon) {
                if (val(a) * val(c) < 0) {
                    b = c;
                } else {
                    a = c;
                }

                c = (a + b) / 2;
            }

            x = c;

            break;
        }
        case 5: {
            if (val(a) * val(x) < 0) {
                x -= (val(x) / (val(x) - val(a))) * (x - a);
            } else {
                x -= (val(x) / (val(b) - val(x))) * (b - x);
            }

            break;
        }
    }

    return x;
}

template<typename T>
std::function<double(T)> NLE<T>::func(
    
) const {
    return _func;
}

template<typename T>
void NLE<T>::set_func(
    std::function<double(T)>& func
) {
    _func = func;
}

template<typename T>
double NLE<T>::val(
    T x
) const {
    return _func(x);
}

template<typename T>
double NLE<T>::diff(
    double x,
    int order,
    double h
) const {
    if (order == 0) {
        return val(x);
    } else {
        return (diff(x + h, order - 1, h) - diff(x, order - 1, h)) / h;
    }
}

template<typename T>
double NLE<T>::solve_iter(
    double x0,
    double a,
    double b,
    std::string method,
    int iter_max,
    double epsilon
) const {
    double x = x0;
    double x_prev = x + ((b - a) / 2 * 1e-10);

    // std::cout << "\n\n" << method << std::endl;

    for (int i = 0; i < iter_max; i++) {
        double x_next = calc_next_x(x, x_prev, a, b, method);

        // std::cout << "i: " << i + 1 << ", norm: " << std::abs(x - x_next) << std::endl;

        if (std::abs(x - x_next) < epsilon) {
            std::cout << i + 1 << " iters of " << iter_max << std::endl;
            break;
        }

        x_prev = x;
        x = x_next;
    }

    return x;
}
