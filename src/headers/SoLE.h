#pragma once

#include <iostream>
#include <string>

#include "../headers/Matrix.h"
#include "../headers/Vector.h"


template <typename T>
class SoLE {
    private:
        Matrix<T> _coeffs;
        Vector<T> _terms;

        void elimination(

        );

        Vector<T> calc_next_x(
            Vector<T> x,
            Vector<T> x_prev,
            std::string method,
            int grad_iter
        ) const;
    public:
        SoLE(
            size_t row_count,
            size_t col_count
        ) : _coeffs(row_count, col_count), _terms(row_count) {

        };

        ~SoLE(

        ) {
            ;
        };

        size_t row_count(
            
        ) const;

        size_t col_count(

        ) const;

        Matrix<T> matrix(

        ) const;
        Vector<T> vector(

        ) const;

        void set_row_count(
            size_t row_count
        );

        void set_col_count(
            size_t col_count
        );

        void set_matrix(
            Matrix<T>& matrix
        );

        void set_matrix(
            std::string file_path,
            const char delimiter = ' '
        );

        void set_vector(
            Vector<T>& vector
        );

        void set_vector(
            std::string file_path,
            const char delimiter = ' '
        );

        friend std::ostream& operator<<(
            std::ostream& out,
            const SoLE& sole
        ) {
            out << sole._coeffs << std::endl;
            out << sole._terms;

            return out;
        };

        Matrix<T> extended(
            
        ) const;

        double det(

        ) const;

        int rank(

        ) const;

        bool is_compatible(

        ) const;

        bool is_diag(

        ) const;

        bool is_diag_d(

        ) const;

        bool is_solvable(
    
        ) const;

        bool is_solution(
            Vector<T>& solution
        ) const;

        void to_diag(
        
        );

        void to_diag_d(
        
        );

        Vector<T> solve_iter(
            std::string method = "si",
            int iter_max = 100,
            double epsilon = 1e-5,
            int grad_iter = 100
        ) const;

        Vector<T> solve(
            std::string method
        ) const;
};
