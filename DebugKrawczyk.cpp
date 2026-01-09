#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "include/interval_krawczyk/KaucherInterval.h"
#include "ik::IntervalVector<2>.h"
#include "ik::IntervalMatrix<2, 2>.h"

int main()
{
    std::cout << "========================================" << std::endl;
    std::cout << "Simple Interval Operation Test" << std::endl;
    std::cout << "========================================" << std::endl;

    KaucherInterval x(1.0, 2.0);
    PSGMDirectedInterval y(0.5, 1.5);

    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;

    auto f = x * x + y * y - PSGMDirectedInterval(4.0, 4.0);
    std::cout << "f(x,y) = x^2 + y^2 - 4 = " << f << std::endl;

    std::cout << "\nTest interval contains point:" << std::endl;
    PSGMDirectedInterval test(1.4, 1.5);
    std::cout << "test = " << test << std::endl;
    std::cout << "test.contains(1.4142) = " << (test.contains(1.4142) ? "true" : "false") << std::endl;
    std::cout << "test.contains(1.6) = " << (test.contains(1.6) ? "true" : "false") << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "Interval Vector Test" << std::endl;
    std::cout << "========================================" << std::endl;

    ik::IntervalVector<2> vec1(2);
    vec1[0] = PSGMDirectedInterval(1.0, 2.0);
    vec1[1] = PSGMDirectedInterval(1.0, 2.0);

    std::cout << "vec1 = " << vec1 << std::endl;
    std::cout << "vec1.maxWidth() = " << vec1.maxWidth() << std::endl;
    std::cout << "vec1.middle() = " << vec1.middle() << std::endl;

    ik::IntervalVector<2> mid = vec1.middle();
    std::cout << "vec1.midpoint() = " << mid << std::endl;

    auto eigen_vec = vec1.midpoint();
    std::cout << "Eigen vector size: " << eigen_vec.size() << std::endl;
    std::cout << "Eigen vector: " << eigen_vec.transpose() << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "Interval Matrix Test" << std::endl;
    std::cout << "========================================" << std::endl;

    ik::IntervalMatrix<2, 2> J(2, 2);
    J(0, 0) = x * PSGMDirectedInterval(2.0);
    J(0, 1) = y * PSGMDirectedInterval(2.0);
    J(1, 0) = PSGMDirectedInterval(1.0);
    J(1, 1) = PSGMDirectedInterval(-1.0);

    std::cout << "J = " << std::endl << J << std::endl;

    auto J_mid = J.midpoint();
    std::cout << "J.midpoint() = " << std::endl << J_mid << std::endl;

    double det = J_mid.determinant();
    std::cout << "det(J.midpoint()) = " << det << std::endl;

    if (std::abs(det) > 1e-12)
    {
        Eigen::MatrixXd J_inv = J_mid.inverse();
        std::cout << "J.midpoint()^-1 = " << std::endl << J_inv << std::endl;

        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
        Eigen::MatrixXd I_minus_YJ = I - J_inv * J.lower();
        std::cout << "I - Y * J.lower() = " << std::endl << I_minus_YJ << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Krawczyk Iteration Test" << std::endl;
    std::cout << "========================================" << std::endl;

    ik::IntervalVector<2> X(2);
    X[0] = PSGMDirectedInterval(1.0, 2.0);
    X[1] = PSGMDirectedInterval(1.0, 2.0);

    std::cout << "初始区间 X = " << X << std::endl;
    std::cout << "X.maxWidth() = " << X.maxWidth() << std::endl;

    int maxIter = 20;
    for (int iter = 0; iter < maxIter; ++iter)
    {
        std::cout << "\nIteration " << iter + 1 << ":" << std::endl;

        ik::IntervalVector<2> c = X.middle();
        std::cout << "  Midpoint c = " << c << std::endl;

        Eigen::VectorXd c_eigen = c.midpoint();
        std::cout << "  c_eigen = " << c_eigen.transpose() << std::endl;

        ik::IntervalVector<2> fX(2);
        fX[0] = X[0] * X[0] + X[1] * X[1] - PSGMDirectedInterval(4.0, 4.0);
        fX[1] = X[0] - X[1];
        std::cout << "  f(X) = " << fX << std::endl;

        ik::IntervalVector<2> fc(2);
        fc[0] = c[0] * c[0] + c[1] * c[1] - PSGMDirectedInterval(4.0, 4.0);
        fc[1] = c[0] - c[1];
        std::cout << "  f(c) = " << fc << std::endl;

        ik::IntervalMatrix<2, 2> JX(2, 2);
        JX[0][0] = X[0] * PSGMDirectedInterval(2.0);
        JX[0][1] = X[1] * PSGMDirectedInterval(2.0);
        JX[1][0] = PSGMDirectedInterval(1.0);
        JX[1][1] = PSGMDirectedInterval(-1.0);
        std::cout << "  J(X) = " << std::endl << JX << std::endl;

        Eigen::MatrixXd J_mid = JX.midpoint();
        double det = J_mid.determinant();
        std::cout << "  det(J(c)) = " << det << std::endl;

        if (std::abs(det) < 1e-12)
        {
            std::cout << "  Jacobian is singular, stopping iteration" << std::endl;
            break;
        }

        Eigen::MatrixXd Y = J_mid.inverse();

        ik::IntervalVector<2> term1(2);
        for (int i = 0; i < 2; ++i)
        {
            double val = c_eigen(i);
            for (int j = 0; j < 2; ++j)
            {
                val -= Y(i, j) * fc[j].middle();
            }
            term1[i] = PSGMDirectedInterval(val, val);
        }
        std::cout << "  term1 = c - Y*f(c) = " << term1 << std::endl;

        ik::IntervalVector<2> X_minus_c(2);
        for (int i = 0; i < 2; ++i)
        {
            double delta = (X[i].upper() - X[i].lower()) / 2.0;
            X_minus_c[i] = PSGMDirectedInterval(-delta, delta);
        }
        std::cout << "  X - c = " << X_minus_c << std::endl;

        Eigen::MatrixXd I_minus_YJ = Eigen::MatrixXd::Identity(2, 2) - Y * JX.lower();
        std::cout << "  I - Y*J.lower() = " << std::endl << I_minus_YJ << std::endl;

        ik::IntervalVector<2> term3(2);
        for (int i = 0; i < 2; ++i)
        {
            PSGMDirectedInterval sum;
            for (int j = 0; j < 2; ++j)
            {
                double coeff = I_minus_YJ(i, j);
                if (coeff >= 0)
                    sum = sum + PSGMDirectedInterval(coeff) * X_minus_c[j];
                else
                    sum = sum + PSGMDirectedInterval(coeff) * X_minus_c[j].dual();
            }
            term3[i] = sum;
        }
        std::cout << "  term3 = (I - Y*J(X)) * (X - c) = " << term3 << std::endl;

        ik::IntervalVector<2> K = term1 + term3;
        std::cout << "  K = term1 + term3 = " << K << std::endl;

        ik::IntervalVector<2> X_new = X.intersection(K);
        std::cout << "  X_new = X ∩ K = " << X_new << std::endl;

        if (X_new.empty())
        {
            std::cout << "  X_new is empty, stopping iteration" << std::endl;
            break;
        }

        bool contains = true;
        for (int i = 0; i < 2; ++i)
        {
            if (!X_new[i].contains(X[i]))
            {
                contains = false;
                break;
            }
        }

        X = X_new;
        std::cout << "  X 更新为 = " << X << std::endl;
        std::cout << "  X.maxWidth() = " << X.maxWidth() << std::endl;

        if (contains)
        {
            std::cout << "  X_new ⊆ X, converged!" << std::endl;
            break;
        }

        if (X.maxWidth() < 1e-8)
        {
            std::cout << "  Interval is small enough, converged!" << std::endl;
            break;
        }
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Final Result" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "解区间 X = " << X << std::endl;
    std::cout << "中点 x* = " << X.middle() << std::endl;
    std::cout << "宽度 = " << X.maxWidth() << std::endl;

    auto mid_vec = X.middle();
    double x_val = mid_vec[0].middle();
    double y_val = mid_vec[1].middle();
    std::cout << "验证: f(" << x_val << ", " << y_val << ") = ";
    std::cout << "(" << x_val*x_val + y_val*y_val - 4 << ", " << x_val - y_val << ")" << std::endl;

    return 0;
}
