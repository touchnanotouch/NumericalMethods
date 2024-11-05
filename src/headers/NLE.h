#pragma once

#include <iostream>
#include <functional>
#include <string>

#include "../headers/Vector.h"


template<typename T>
class NLE {
    private:
        std::function<T(T)> _func;

        double calc_next_x(
            T x,
            T x_prev,
            double a,
            double b,
            std::string method
        ) const;
    public:
        NLE(
            std::function<T(T)> func
        ) {
            _func = func;
        };

        ~NLE() {

        };

        std::function<T(T)> func(

        ) const;

        void set_func(
            std::function<T(T)> func
        );

        double diff(
            T x,
            int n,
            double h = 1e-5
        ) const;

        double solve_iter(
            T x0,
            double a,
            double b,
            std::string method = "si",
            int iter_max = 100,
            double epsilon = 1e-5
        );
};
