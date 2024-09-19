#pragma once

#include <iostream>


template <typename T>
class SoLE {
    private:
        Matrix<T> _coeffs;
        Vector<T> _terms;
    public:
        SoLE(
            size_t row_count,
            size_t col_count
        ) : _coeffs(row_count, col_count), _terms(row_count) {

        };

        ~SoLE() {

        };

        Matrix<T> get_matrix(

        ) const;
        Vector<T> get_vector(

        ) const;

        void set_matrix(
            Matrix<T>& matrix
        );

        void set_vector(
            Vector<T>& vector
        );

        friend std::ostream& operator<<(
            std::ostream& out,
            const SoLE& sole
        ) {
            out << sole._coeffs << std::endl;
            out << sole._terms << std::endl;

            return out;
        };

        bool is_solution(
            Vector<T>& solution
        ) const;
};
