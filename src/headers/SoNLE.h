#pragma once

#include <string>

#include "../headers/Matrix.h"


template<typename T>
class SoNLE {
    private:
        Matrix<T> _funcs_mat;
    
    public:
        SoNLE(
            size_t row_count,
            size_t col_count
        ) : _funcs_mat(row_count, col_count) {

        };

        ~SoNLE(

        ) {
            ;
        };

        size_t row_count(

        ) const;

        Matrix<T> funcs_matrix(

        ) const;

        void set_row_count(
            size_t row_count
        );

        void set_funcs_matrix(
            Matrix<T>& matrix
        );

        void set_funcs_matrix(
            std::string file_path,
            const char delimiter = ' '
        );

        double val(
            int index,
            double x
        );

        double diff(
            int index,
            double x,
            int n,
            double h = 1e-5
        )
};
