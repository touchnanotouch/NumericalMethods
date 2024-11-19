#pragma once

#include <string>

#include "../headers/Vector.h"


template<typename T>
class NLE {
    private:
        Vector<T> _func_vec;

        double calc_next_x(
            double x,
            double x_prev,
            double a,
            double b,
            std::string method
        ) const;
    public:
        NLE(
            size_t row_count
        ) : _func_vec(row_count) {
            ;
        };

        ~NLE(
            
        ) {
            ;
        };

        size_t row_count(

        ) const;

        Vector<T> func_vector(

        ) const;

        void set_row_count(
            size_t row_count
        );

        void set_func_vector(
            Vector<T>& vector
        );

        void set_func_vector(
            std::string file_path,
            const char delimiter = ' '
        );

        double val(
            double x
        ) const;

        double diff(
            double x,
            int n,
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
