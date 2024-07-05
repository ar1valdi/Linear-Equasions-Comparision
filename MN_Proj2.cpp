﻿#include <iostream>
#include "matplotlib-cpp/matplotlibcpp.h"
#include "Matrix.h"
#include "MatrixCreator.h"
#include "Consts.h"
#include "Operators.h"
#include "ExReturnData.h"
#include <ctime>
#include <list>

namespace plt = matplotlibcpp;
namespace cns = Consts;
using namespace std;

Matrix calc_residuum(const Matrix& A, const Matrix& x, const Matrix& b) {
    return A * x - b;
}

ExReturnData solve(const Matrix& A, const Matrix& b, int max_iter, int method) {
    pair<int, int> size_b = b.getSize();
    pair<int, int> size_A = A.getSize();

    if (size_b.first != 1) {
        throw new runtime_error("b is not a column not vector");
    }
    if (size_A.second != size_b.second || size_A.first != size_A.second) {
        throw new runtime_error("incompatible sizes");
    }

    int x_m = size_b.second;
    Matrix x(1, x_m, 1);
    int iter = 0;
    double err = cns::RESIDUUM_ACCEPTANCE + 1;
    list<double> errors;

    clock_t start = clock();
    while (iter != max_iter && err > cns::RESIDUUM_ACCEPTANCE) {
        Matrix x_old = x;
        Matrix* x_to_use_sum1 = method == cns::JACOBI ? &x_old : &x;
        for (int i = 0; i < x_m; i++) {
            double sum1, sum2;
            sum1 = sum2 = 0;

            for (int j = 0; j < i; j++) {
                sum1 += A[j][i] * (*x_to_use_sum1)[0][j];
            }
            for (int j = i + 1; j < x_m; j++) {
                sum2 += A[j][i] * x_old[0][j];
            }

            x[0][i] = (b[0][i] - sum1 - sum2) / A[i][i];
        }

        err = calc_residuum(A, x, b).get_norm();
        errors.push_back(err);
        iter++;
    }
    clock_t stop = clock();
    double seconds = double(stop - start) / CLOCKS_PER_SEC;

    return { x, errors, iter, seconds };
}

void LU_factorization(const Matrix& A, Matrix& L, Matrix& U) {
    int n = A.getSize().first;
    U = A;
    L = MatrixCreator::diagonal_1(n, 1);
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            L[i][j] = U[i][j] / U[j][j];
            U[i] = U[i] - U[j] * L[i][j];
        }
    }
}

ExReturnData solve_LU(const Matrix& A, const Matrix& b) {
    int n = A.getSize().first;
    Matrix L, U;
    clock_t start = clock();
    LU_factorization(A, L, U);

    Matrix y(1, n);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[0][j];
        }
        y[0][i] = (b[0][i] - sum) / L[i][i];
    }

    Matrix x(1, n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[0][j];
        }
        x[0][i] = (y[0][i] - sum) / U[i][i];
    }
    clock_t stop = clock();
    double time = double(stop - start) / CLOCKS_PER_SEC;

    return { x, {calc_residuum(A, x, b).get_norm()}, 1, time};
}

void zad_ab() {
    Matrix A = MatrixCreator::diagonal_5(cns::N, cns::A1_1, cns::A2_1, cns::A3_1);
    Matrix B = MatrixCreator::vector_sin_n_f_1(cns::N);

    ExReturnData x1 = solve(A, B, cns::MAX_ITER, cns::JACOBI);
    ExReturnData x2 = solve(A, B, cns::MAX_ITER, cns::GAUSS_SEIDEL);
    ExReturnData x3 = solve_LU(A, B);
    vector<int> y1(x1.errors.size());
    vector<int> y2(x2.errors.size());

    vector<double> err_x1(x1.errors.size());
    vector<double> err_x2(x2.errors.size());
    int i = 0;
    for (double x : x1.errors) {
        err_x1[i++] = x;
    }
    i = 0;
    for (double x : x2.errors) {
        err_x2[i++] = x;
    }
    for (int i = 0; i < y1.size(); i++) {
        y1[i] = i + 1;
    }
    for (int i = 0; i < y2.size(); i++) {
        y2[i] = i + 1;
    }

    cout << "Zadanie A i B:\n";
    cout << "Jacobi       | iters=" << x1.iters << " | time=" << x1.time << "| error=" << *(--x1.errors.end()) << '\n';
    cout << "Gauss-Seidel | iters=" << x2.iters << " | time=" << x2.time << "| error=" << *(--x2.errors.end()) << '\n';
    cout << "xd: " << *x3.errors.begin();

    plt::named_semilogy("Jacobi", y1, err_x1, "");
    plt::named_semilogy("Gauss-Seidel", y2, err_x2, "");
    plt::title("Errors on subsequent iterations");
    plt::xlabel("Size");
    plt::ylabel("Time");
    plt::legend();
    plt::show();
}

void zad_c() {
    Matrix A = MatrixCreator::diagonal_5(cns::N, cns::A1_2, cns::A2_2, cns::A3_2);
    Matrix B = MatrixCreator::vector_sin_n_f_1(cns::N);

    ExReturnData x1 = solve(A, B, cns::MAX_ITER, cns::JACOBI);
    ExReturnData x2 = solve(A, B, cns::MAX_ITER, cns::GAUSS_SEIDEL);
    vector<int> y1(x1.errors.size());
    vector<int> y2(x2.errors.size());

    vector<double> err_x1(x1.errors.size());
    vector<double> err_x2(x2.errors.size());
    int i = 0;
    for (double x : x1.errors) {
        err_x1[i++] = x;
    }
    i = 0;
    for (double x : x2.errors) {
        err_x2[i++] = x;
    }
    for (int i = 0; i < y1.size(); i++) {
        y1[i] = i + 1;
    }
    for (int i = 0; i < y2.size(); i++) {
        y2[i] = i + 1;
    }

    cout << "ZADANIE C:\nNie\n";
    cout << "Jacobi: " << err_x1[99] << " " << x1.time << "\n";
    cout << "Gauss-Seidel: ERR=" << err_x2[99] << " TIME=" << x2.time << "\n";
    plt::named_semilogy("Jacobi", y1, err_x1, "");
    plt::named_semilogy("Gauss-Seidel", y2, err_x2, "");
    plt::title("Errors on subsequent iterations");
    plt::xlabel("Size");
    plt::ylabel("Time");
    plt::legend();
    plt::show();
}

void zad_d() {
    Matrix A = MatrixCreator::diagonal_5(cns::N, cns::A1_2, cns::A2_2, cns::A3_2);
    Matrix B = MatrixCreator::vector_sin_n_f_1(cns::N);

    Matrix x = solve_LU(A, B).m;
    cout << "ZADANIE D:\nNorma wektora bledu = " << calc_residuum(A, x, B).get_norm() << '\n';
}

void zad_e() {
    list<double> measured_times;
    for (int N : cns::LENGHTS_TO_COMPARE) {
        Matrix A = MatrixCreator::diagonal_5(N, cns::A1_1, cns::A2_1, cns::A3_1);
        Matrix B = MatrixCreator::vector_sin_n_f_1(N);

        ExReturnData x_jc = solve(A, B, cns::MAX_ITER, cns::JACOBI);
        ExReturnData x_gs = solve(A, B, cns::MAX_ITER, cns::GAUSS_SEIDEL);
        ExReturnData x_lu = solve_LU(A, B);

        measured_times.push_back(x_jc.time);
        measured_times.push_back(x_gs.time);
        measured_times.push_back(x_lu.time);

        cout << "solved for N = " << N << '\n';
    }

    int samples = measured_times.size() / 3;
    vector<double> jc_times(samples);
    vector<double> gs_times(samples);
    vector<double> lu_times(samples);
    vector<int> plot_x(samples);

    int i = 0;
    for (auto it = measured_times.begin(); it != measured_times.end(); it++) {
        jc_times[i] = *(it++);
        gs_times[i] = *(it++);
        lu_times[i++] = *(it);
    }
    i = 0;
    for (int N : cns::LENGHTS_TO_COMPARE) {
        plot_x[i++] = N;
    }

    cout << "Jacobi:\n";
    for (double t : jc_times) {
        cout << t << '\n';
    }
    cout << "\nGauss-Seidel:\n";
    for (double t : gs_times) {
        cout << t << '\n';
    }
    cout << "LU:\n";
    for (double t : lu_times) {
        cout << t << '\n';
    }

    plt::named_plot("Jacobi", plot_x, jc_times);
    plt::named_plot("Gauss-Seidel", plot_x, gs_times);
    plt::named_plot("LU", plot_x, lu_times);
    plt::title("Time required to solve A*x=b for different matricies sizes");
    plt::xlabel("Size");
    plt::ylabel("Time");
    plt::legend();
    plt::show();
}

int main() {
    plt::backend("WXAgg");
    //zad_ab();
    //zad_c();
    zad_d();
    zad_e();
    return 0;
}

