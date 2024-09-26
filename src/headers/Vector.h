#pragma once

#include <iostream>


template<typename T>
class Vector {
    private:
        size_t _row_cnt;
        T* _vec;

    public:
        Vector(
            size_t row_count
        ) : _row_cnt(row_count) {
            _vec = new T[_row_cnt];
            for (size_t i = 0; i < _row_cnt; i++) {
                _vec[i] = 0;
            }
        };

        ~Vector(

        ) {
            delete[] _vec;
        };

        size_t row_count(

        ) const;

        T* vec(

        ) const;

        void set_row_count(
            size_t row_count
        );

        void set_vector(
            T* vector
        );

        void swap(
            T& a,
            T& b
        );

        Vector(
            const Vector& other
        );

        Vector& operator=(
            const Vector& other
        );

        Vector(
            Vector&& other
        );

        Vector& operator=(
            Vector&& other
        );

        friend std::ostream& operator<<(
            std::ostream& out,
            const Vector& vector
        ) {
            for (size_t i = 0; i < vector._row_cnt; i++) {
                out << vector._vec[i] << " ";
            }
            out << std::endl;

            return out;
        }

        T& operator[](
            size_t index
        ) const;
        
        bool operator==(
            const Vector& other
        ) const;

        Vector operator+(
            const Vector& other
        ) const;

        Vector operator-(
            const Vector& other
        ) const;

        Vector operator*(
            const T& scalar
        ) const;

        Vector operator/(
            const T& scalar
        ) const;

        double norm(
            
        ) const;
};
