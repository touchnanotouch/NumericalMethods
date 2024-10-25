#include <functional>
#include <iostream>
#include <chrono>
#include <sstream>
#include <cmath>

#include "../headers/Vector.h"
#include "../headers/Tester.h"


template class Tester<int>;
template class Tester<double>;


template<typename T>
double Tester<T>::execution_time(
    std::function<void()> func
) const {
    auto start_time = std::chrono::high_resolution_clock::now();
    func();
    auto end_time = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
}

template<typename T>
double Tester<T>::solution_accuracy(
    Vector<T> sol,
    Vector<T> actual
) const {
    size_t n = sol.row_count();

    double rmse = 0.0;

    for (size_t i = 0; i < n; i++) {
        rmse += (sol[i] - actual[i]) * (sol[i] - actual[i]);
    }

    rmse = std::sqrt(rmse / n);

    double accuracy_percent = std::max(0.0, (1 - (rmse / actual.norm())) * 100);

    return accuracy_percent;
}

template<typename T>
std::string Tester<T>::solution_iter(
    
) const {
    Vector<T> si_sol = _sole.solve_iter("si");
    Vector<T> seidel_sol = _sole.solve_iter("seidel");
    Vector<T> rel_sol = _sole.solve_iter("rel");
    Vector<T> res_sol = _sole.solve_iter("res");
    Vector<T> grad_sol = _sole.solve_iter("grad");

    std::stringstream ss;

    ss << "\nIterative methods (solution)\n" << std::endl;
    ss << "si: " << si_sol;
    ss << "seidel: " << seidel_sol;
    ss << "rel: " << rel_sol;
    ss << "res: " << res_sol;
    ss << "grad: " << grad_sol;

    return ss.str();
}

template<typename T>
std::string Tester<T>::solution_exact(
    
) const {
    Vector<T> gauss_sol = _sole.solve("gauss");
    Vector<T> thomas_sol = _sole.solve("thomas");
    Vector<T> lu_sol = _sole.solve("lu");
    Vector<T> qroots_sol = _sole.solve("qroots");

    std::stringstream ss;

    ss << "\nExact methods (solution)\n" << std::endl;
    ss << "gauss: " << gauss_sol;
    ss << "thomas: " << thomas_sol;
    ss << "lu: " << lu_sol;
    ss << "qroots: " << qroots_sol;

    return ss.str();
}

template<typename T>
std::string Tester<T>::solution_all(

) const {
    std::stringstream ss;

    ss << solution_iter();
    ss << solution_exact();

    return ss.str();
}

template<typename T>
std::string Tester<T>::time_iter(
    
) const {
    double si_time = execution_time([&](){ _sole.solve_iter("si"); });
    double seidel_time = execution_time([&](){ _sole.solve_iter("seidel"); });
    double rel_time = execution_time([&](){ _sole.solve_iter("rel"); });
    double res_time = execution_time([&](){ _sole.solve_iter("res"); });
    double grad_time = execution_time([&](){ _sole.solve_iter("grad"); });

    std::stringstream ss;

    ss << "\nIterative methods (execution time)\n" << std::endl;
    ss << "si: " << si_time << "s" << std::endl;
    ss << "seidel: " << seidel_time << "s" << std::endl;
    ss << "rel: " << rel_time << "s" << std::endl;
    ss << "res: " << res_time << "s" << std::endl;
    ss << "grad: " << grad_time << "s" << std::endl;

    return ss.str();
}

template<typename T>
std::string Tester<T>::time_exact(
    
) const {
    double gauss_time = execution_time([&](){ _sole.solve("gauss"); });
    double thomas_time = execution_time([&](){ _sole.solve("thomas"); });
    double lu_time = execution_time([&](){ _sole.solve("lu"); });
    double qroots_time = execution_time([&](){ _sole.solve("qroots"); });

    std::stringstream ss;

    ss << "\nExact methods (execution time)\n" << std::endl;
    ss << "gauss: " << gauss_time << "s" << std::endl;
    ss << "thomas: " << thomas_time << "s" << std::endl;
    ss << "lu: " << lu_time << "s" << std::endl;
    ss << "qroots: " << qroots_time << "s" << std::endl;

    return ss.str();
}

template<typename T>
std::string Tester<T>::time_all(
    
) const {
    std::stringstream ss;

    ss << time_iter();
    ss << time_exact();

    return ss.str();
}

template<typename T>
std::string Tester<T>::accuracy_iter(
    Vector<T> actual
) const {
    Vector<T> si_sol = _sole.solve_iter("si");
    Vector<T> seidel_sol = _sole.solve_iter("seidel");
    Vector<T> rel_sol = _sole.solve_iter("rel");
    Vector<T> res_sol = _sole.solve_iter("res");
    Vector<T> grad_sol = _sole.solve_iter("grad");

    std::stringstream ss;

    ss << "\nIterative methods (accuracy)\n" << std::endl;
    ss << "si: " << solution_accuracy(si_sol, actual) << std::endl;
    ss << "seidel: " << solution_accuracy(seidel_sol, actual) << std::endl;
    ss << "rel: " << solution_accuracy(rel_sol, actual) << std::endl;
    ss << "res: " << solution_accuracy(res_sol, actual) << std::endl;
    ss << "grad: " << solution_accuracy(grad_sol, actual) << std::endl;

    return ss.str();
}

template<typename T>
std::string Tester<T>::accuracy_exact(
    Vector<T> actual
) const {
    Vector<T> gauss_sol = _sole.solve("gauss");
    Vector<T> thomas_sol = _sole.solve("thomas");
    Vector<T> lu_sol = _sole.solve("lu");
    Vector<T> qroots_sol = _sole.solve("qroots");

    std::stringstream ss;

    ss << "\nExact methods (accuracy)\n" << std::endl;
    ss << "gauss: " << solution_accuracy(gauss_sol, actual) << std::endl;
    ss << "thomas: " << solution_accuracy(thomas_sol, actual) << std::endl;
    ss << "lu: " << solution_accuracy(lu_sol, actual) << std::endl;
    ss << "qroots: " << solution_accuracy(qroots_sol, actual) << std::endl;

    return ss.str();
}

template<typename T>
std::string Tester<T>::accuracy_all(
    Vector<T> actual
) const {
    std::stringstream ss;

    ss << accuracy_iter(actual);
    ss << accuracy_exact(actual);

    return ss.str();
}

template<typename T>
std::string Tester<T>::measure_iter(
    Vector<T> actual
) const {
    auto start_time = std::chrono::high_resolution_clock::now();
    Vector<T> si_sol = _sole.solve_iter("si");
    auto end_time = std::chrono::high_resolution_clock::now();

    double si_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    auto start_time1 = std::chrono::high_resolution_clock::now();
    Vector<T> seidel_sol = _sole.solve_iter("seidel");
    auto end_time1 = std::chrono::high_resolution_clock::now();

    double seidel_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time1 - start_time1).count();

    auto start_time2 = std::chrono::high_resolution_clock::now();
    Vector<T> rel_sol = _sole.solve_iter("rel");
    auto end_time2 = std::chrono::high_resolution_clock::now();

    double rel_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time2 - start_time2).count();

    auto start_time3 = std::chrono::high_resolution_clock::now();
    Vector<T> res_sol = _sole.solve_iter("res");
    auto end_time3 = std::chrono::high_resolution_clock::now();

    double res_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time3 - start_time3).count();

    auto start_time4 = std::chrono::high_resolution_clock::now();
    Vector<T> grad_sol = _sole.solve_iter("grad");
    auto end_time4 = std::chrono::high_resolution_clock::now();

    double grad_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time4 - start_time4).count();

    std::stringstream ss;

    ss << "\nIterative methods (execution time)\n" << std::endl;

    ss << "si: " << si_time << "s" << std::endl;
    ss << "seidel: " << seidel_time << "s" << std::endl;
    ss << "rel: " << rel_time << "s" << std::endl;
    ss << "res: " << res_time << "s" << std::endl;
    ss << "grad: " << grad_time << "s" << std::endl;

    ss << "\nIterative methods (accuracy)\n" << std::endl;

    ss << "si: " << solution_accuracy(si_sol, actual) << "%" << std::endl;
    ss << "seidel: " << solution_accuracy(seidel_sol, actual) << "%" << std::endl;
    ss << "rel: " << solution_accuracy(rel_sol, actual) << "%" << std::endl;
    ss << "res: " << solution_accuracy(res_sol, actual) << "%" << std::endl;
    ss << "grad: " << solution_accuracy(grad_sol, actual) << "%" << std::endl;

    return ss.str();
}

template<typename T>
std::string Tester<T>::measure_exact(
    Vector<T> actual
) const {
    auto start_time = std::chrono::high_resolution_clock::now();
    Vector<T> gauss_sol = _sole.solve("gauss");
    auto end_time = std::chrono::high_resolution_clock::now();

    double gauss_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    auto start_time1 = std::chrono::high_resolution_clock::now();
    Vector<T> thomas_sol = _sole.solve("thomas");
    auto end_time1 = std::chrono::high_resolution_clock::now();

    double thomas_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time1 - start_time1).count();

    auto start_time2 = std::chrono::high_resolution_clock::now();
    Vector<T> lu_sol = _sole.solve("lu");
    auto end_time2 = std::chrono::high_resolution_clock::now();

    double lu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time2 - start_time2).count();

    auto start_time3 = std::chrono::high_resolution_clock::now();
    Vector<T> qroots_sol = _sole.solve("qroots");
    auto end_time3 = std::chrono::high_resolution_clock::now();

    double qroots_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time3 - start_time3).count();

    std::stringstream ss;

    ss << "\nExact methods (execution time)\n" << std::endl;

    ss << "gauss: " << gauss_time << "s" << std::endl;
    ss << "thomas: " << thomas_time << "s" << std::endl;
    ss << "lu: " << lu_time << "s" << std::endl;
    ss << "qroots: " << qroots_time << "s" << std::endl;

    ss << "\nExact methods (accuracy)\n" << std::endl;

    ss << "gauss: " << solution_accuracy(gauss_sol, actual) << "%" << std::endl;
    ss << "thomas: " << solution_accuracy(thomas_sol, actual) << "%" << std::endl;
    ss << "lu: " << solution_accuracy(lu_sol, actual) << "%" << std::endl;
    ss << "qroots: " << solution_accuracy(qroots_sol, actual) << "%" << std::endl;

    return ss.str();
}

template<typename T>
std::string Tester<T>::measure_all(
    Vector<T> actual
) const {
    std::stringstream ss;

    ss << measure_iter(actual);
    ss << measure_exact(actual);

    return ss.str();
}
