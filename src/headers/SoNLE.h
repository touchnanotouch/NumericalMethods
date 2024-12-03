#pragma once

#include <string>
#include <functional>

#include "../headers/Vector.h"
#include "../headers/Matrix.h"


template<typename T>
class SoNLE {
    private:
        size_t _row_cnt;

        std::function<double(T*)>* _funcs;

        Vector<double> calc_next_x(
            Vector<double> x,
            Vector<double> x_prev,
            std::string method
        ) const;
    public:
        SoNLE(
            size_t row_count,
            std::function<double(T*)>* funcs
        ) : _row_cnt(row_count) {
            _funcs = new std::function<double(T*)>[row_count];
            for (size_t i = 0; i < row_count; i++) {
                _funcs[i] = funcs[i];
            }
        };

        ~SoNLE(

        ) {
            delete[] _funcs;
        };

        size_t row_count(

        ) const;

        std::function<double(T*)>* funcs(

        ) const;

        void set_row_count(
            size_t row_count
        );

        void set_funcs(
            std::function<double(T*)>* funcs
        );

        Vector<double> vals(
            Vector<T> x
        ) const;

        double diff(
            Vector<double> x,
            int order,
            size_t func_index,
            size_t var_index,
            double h = 1e-5
        ) const;

        Matrix<double> jacobian(
            Vector<double> x,
            double h = 1e-5
        ) const;

        Vector<double> solve_iter(
            Vector<double> x0,
            std::string method = "si",
            int iter_max = 100,
            double epsilon = 1e-5
        ) const;
};
