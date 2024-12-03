#pragma once

#include <string>
#include <functional>

#include "../headers/Vector.h"


template<typename T>
class NLE {
    private:
        std::function<double(T)> _func;

        double calc_next_x(
            double x,
            double x_prev,
            double a,
            double b,
            std::string method
        ) const;
    public:
        NLE(
            std::function<double(T)> func
        ) : _func(func) {
            ;
        };

        ~NLE(
            
        ) {
            ;
        };

        std::function<double(T)> func(

        ) const;

        void set_func(
            std::function<double(T)>& func
        );

        double val(
            T x
        ) const;

        double diff(
            double x,
            int order,
            double h = 1e-5
        ) const;

        double solve_iter(
            double x0,
            double a,
            double b,
            std::string method = "si",
            int iter_max = 100,
            double epsilon = 1e-5
        ) const;
};
