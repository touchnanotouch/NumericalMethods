#pragma once

#include <string>
#include <functional>

#include "../headers/Vector.h"
#include "../headers/SoLE.h"


template<typename T>
class Tester {
    private:
        size_t _n;

        std::string _m_path;
        std::string _v_path;

        SoLE<T> _sole;

        double execution_time(
            std::function<void()> func
        ) const;

        double solution_accuracy(
            Vector<T> sol,
            Vector<T> actual
        ) const;
    public:
        Tester(
            size_t n,
            std::string matrix_path,
            std::string vector_path
        ) : _n(n), _sole(n, n) {
            _m_path = matrix_path;
            _v_path = vector_path;

            _sole.set_matrix(_m_path);
            _sole.set_vector(_v_path);
        };

        ~Tester(
            
        ) {
            ;
        };

        std::string solution_iter(

        ) const;

        std::string solution_exact(

        ) const;

        std::string solution_all(

        ) const;

        std::string time_iter(

        ) const;

        std::string time_exact(

        ) const;

        std::string time_all(

        ) const;

        std::string accuracy_iter(
            Vector<T> actual
        ) const;

        std::string accuracy_exact(
            Vector<T> actual
        ) const;

        std::string accuracy_all(
            Vector<T> actual
        ) const;

        std::string measure_iter(
            Vector<T> actual
        ) const;

        std::string measure_exact(
            Vector<T> actual
        ) const;

        std::string measure_all(
            Vector<T> actual
        ) const;
};
