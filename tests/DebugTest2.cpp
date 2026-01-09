#include <iostream>
#include <iomanip>
#include "FixedIntervalKrawczyk.h"

void debugTest()
{
    std::cout << "\n=== Debug Test 2 ===" << std::endl;
    std::cout << "Equation: x^2 + y = 3, x*y = 1" << std::endl;

    IntervalVec<2> x;
    x[0] = Interval(0.5, 1.5);
    x[1] = Interval(1.5, 2.5);

    std::cout << "\n1. Test function calculation:" << std::endl;
    auto f = [](const IntervalVec<2>& x) -> IntervalVec<2>
    {
        IntervalVec<2> r;
        r[0] = x[0] * x[0] + x[1] - Interval(3.0, 3.0);
        r[1] = x[0] * x[1] - Interval(1.0, 1.0);
        return r;
    };

    auto f_x = f(x);
    std::cout << "f(x) = " << f_x << std::endl;

    auto mid = x.midPoint();
    std::cout << "x.mid() = (" << mid(0) << ", " << mid(1) << ")" << std::endl;
    IntervalVec<2> x_mid = IntervalVec<2>::fromPoint(mid);
    auto f_mid = f(x_mid);
    std::cout << "f(mid) = " << f_mid << std::endl;

    std::cout << "\n2. Test Jacobian matrix:" << std::endl;
    auto jacobian = [](const IntervalVec<2>& x) -> IntervalMat<2, 2>
    {
        IntervalMat<2, 2> J;
        J(0, 0) = x[0] * Interval(2.0, 2.0);
        J(0, 1) = Interval(1.0, 1.0);
        J(1, 0) = x[1];
        J(1, 1) = x[0];
        return J;
    };

    auto J = jacobian(x);
    std::cout << "J(x) = " << J << std::endl;
    auto J_mid = jacobian(x_mid);
    std::cout << "J(mid) = " << std::endl << J_mid.midPoint() << std::endl;

    double det = J_mid.midPoint().determinant();
    std::cout << "det(J(mid)) = " << det << std::endl;

    if (std::abs(det) > 1e-12)
    {
        Eigen::MatrixXd J_inv = J_mid.midPoint().inverse();
        std::cout << "J(mid)^-1 = " << std::endl << J_inv << std::endl;

        Eigen::VectorXd c = x.midPoint();
        Eigen::VectorXd f_mid_eigen;
        f_mid_eigen << f_mid[0].mid(), f_mid[1].mid();

        Eigen::VectorXd term1 = c - J_inv * f_mid_eigen;
        std::cout << "term1 = c - Y*f(c) = (" << term1(0) << ", " << term1(1) << ")" << std::endl;
    }

    std::cout << "\n3. Test I - YJ:" << std::endl;
    Eigen::MatrixXd J_lower = J.lower();
    std::cout << "J.lower() = " << std::endl << J_lower << std::endl;

    Eigen::MatrixXd J_inv = J_mid.midPoint().inverse();
    Eigen::MatrixXd I_minus_YJ_lower = Eigen::MatrixXd::Identity(2, 2) - J_inv * J_lower;
    std::cout << "I - Y*J.lower() = " << std::endl << I_minus_YJ_lower << std::endl;

    std::cout << "\n4. Test interval multiplied by coefficient matrix:" << std::endl;
    IntervalMat<2, 2> I_minus_YJ;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            I_minus_YJ(i, j) = Interval(I_minus_YJ_lower(i, j));

    IntervalVec<2> box_minus_c;
    for (int i = 0; i < 2; ++i)
    {
        double delta = x[i].width() / 2.0;
        box_minus_c[i] = Interval(-delta, delta);
    }
    std::cout << "box - c = " << box_minus_c << std::endl;

    IntervalVec<2> term3 = I_minus_YJ * box_minus_c;
    std::cout << "term3 = " << term3 << std::endl;

    std::cout << "\n=== Debug Completed ===" << std::endl;
}

int main()
{
    debugTest();
    return 0;
}
