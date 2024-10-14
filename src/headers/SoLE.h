#pragma once

#include <iostream>
#include <string>

#include "../headers/Matrix.h"
#include "../headers/Vector.h"


template <typename T>
class SoLE {
    private:
        size_t _row_cnt;
        size_t _col_cnt;
        Matrix<T> _coeffs;
        Vector<T> _terms;

        void elimination(

        );

        Vector<T> calc_next_x(
            Vector<T> x,
            std::string method,
            double res_omega,
            int grad_iter
        ) const;
    public:
        SoLE(
            size_t row_count,
            size_t col_count
        ) : _coeffs(row_count, col_count), _terms(row_count) {
            _row_cnt = row_count;
            _col_cnt = col_count;
        }

        ~SoLE() {

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

        );

        int rank(

        );

        bool is_compatible(

        );

        bool is_solution(
            Vector<T>& solution
        ) const;

        bool is_diag_d(

        );

        bool is_solvable(
    
        );

        void to_diag_d(
        
        );

        Vector<T> solve_iter(
            std::string method = "si",
            int iter_max = 100,
            double epsilon = 1e-5,
            double res_omega = 1.0,
            int grad_iter = 1000
        ) const;

        Vector<T> solve(
            std::string method
        ) const;
};
