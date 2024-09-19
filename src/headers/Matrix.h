#pragma once

#include <iostream>

#include "../headers/Vector.h"


template<typename T>
class Matrix {
    private:
        size_t _row_cnt;
        size_t _col_cnt;
        T** _mat;

        Matrix gaussian_elimination(
            Matrix& matrix
        );
    public:
        Matrix(
            size_t row_count,
            size_t col_count
        ) : _row_cnt(row_count), _col_cnt(col_count) {
            _mat = new T*[_row_cnt];
            for (size_t i = 0; i < _row_cnt; i++) {
                _mat[i] = new T[_col_cnt];
                for (size_t j = 0; j < _col_cnt; j++) {
                    _mat[i][j] = 0;
                }
            }
        };

        ~Matrix(

        ) {
            for (size_t i = 0; i < _row_cnt; i++) {
                delete[] _mat[i];
            }

            delete[] _mat;
        };

        Matrix(
            const Matrix& other
        );

        Matrix& operator=(
            const Matrix& other
        );

        Matrix(
            Matrix&& other
        );

        Matrix& operator=(
            Matrix&& other
        );

        size_t get_row_count(

        ) const;

        size_t get_col_count(

        ) const;

        T** get_matrix(
            
        ) const;

        friend std::ostream& operator<<(
            std::ostream& out,
            const Matrix& matrix
        ) {
            for (size_t i = 0; i < matrix._row_cnt; i++) {
                for (size_t j = 0; j < matrix._col_cnt; j++) {
                    out << matrix._mat[i][j] << " ";
                }
                out << std::endl;
            }

            return out;
        };

        T*& operator[](
            int index
        ) const;

        bool operator==(
            const Matrix& other
        ) const;

        Matrix operator+(
            const Matrix& other
        ) const;

        Matrix operator-(
            const Matrix& other
        ) const;

        Matrix operator*(
            const Matrix& other
        ) const;

        Vector<T> operator*(
            const Vector<T>& vector
        ) const;

        Matrix operator*(
            const T& scalar
        ) const;

        Matrix operator/(
            const T& scalar
        ) const;

        void transpose(

        );

        Matrix transposed(
            
        ) const;

        double determinant(

        );

        int rank(

        );
};
