#include <iostream>
#include <iomanip>
#include "FixedIntervalKrawczyk.h"

template<size_t N>
void printResult(const std::string& title, bool success, const std::string& msg,
                 const IntervalVec<N>& sol, double width)
{
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << title << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    if (success)
    {
        std::cout << "✓ Success!" << std::endl;
        std::cout << "Solution Interval: " << sol << std::endl;
        std::cout << "Width: " << std::scientific << width << std::endl;
        std::cout << "Midpoint: (";
        for (size_t i = 0; i < N; ++i)
        {
            std::cout << sol[i].mid();
            if (i < N - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;
    }
    else
    {
        std::cout << "✗ Failed: " << msg << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
}

void testSimpleSystem()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Test 1: x^2 + y^2 = 4, x = y" << std::endl;
    std::cout << "Expected Solution: (1.414, 1.414)" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    auto f = [](const IntervalVec<2>& x) -> IntervalVec<2>
    {
        IntervalVec<2> r;
        r[0] = x[0] * x[0] + x[1] * x[1] - Interval(4.0, 4.0);
        r[1] = x[0] - x[1];
        return r;
    };

    auto jacobian = [](const IntervalVec<2>& x) -> IntervalMat<2, 2>
    {
        IntervalMat<2, 2> J;
        J(0, 0) = x[0] * Interval(2.0, 2.0);
        J(0, 1) = x[1] * Interval(2.0, 2.0);
        J(1, 0) = Interval(1.0, 1.0);
        J(1, 1) = Interval(-1.0, -1.0);
        return J;
    };

    FixedKrawczykSolver<2> solver(f, jacobian, 1e-10, 100, 100);

    IntervalVec<2> initial;
    initial[0] = Interval(0.5, 2.5);
    initial[1] = Interval(0.5, 2.5);

    std::cout << "Initial Interval: " << initial << std::endl;
    auto result = solver.solve(initial);
    printResult<2>("Result", result.success, result.message, result.solution, result.finalWidth);

    if (result.success)
    {
        double x = result.solution[0].mid();
        double y = result.solution[1].mid();
        std::cout << "验证: f(" << x << ", " << y << ") = ("
                  << (x*x + y*y - 4) << ", " << (x - y) << ")" << std::endl;
    }
}

void testCubic()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Test 2: x^3 - x - 1 = 0 (Single Variable)" << std::endl;
    std::cout << "Expected Solution: x ≈ 1.3247" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    auto f = [](const IntervalVec<1>& x) -> IntervalVec<1>
    {
        IntervalVec<1> r;
        r[0] = x[0] * x[0] * x[0] - x[0] - Interval(1.0, 1.0);
        return r;
    };

    auto jacobian = [](const IntervalVec<1>& x) -> IntervalMat<1, 1>
    {
        IntervalMat<1, 1> J;
        J(0, 0) = x[0] * x[0] * Interval(3.0, 3.0) - Interval(1.0, 1.0);
        return J;
    };

    FixedKrawczykSolver<1> solver(f, jacobian, 1e-12, 100, 100);

    IntervalVec<1> initial;
    initial[0] = Interval(1.0, 2.0);

    std::cout << "Initial Interval: " << initial << std::endl;
    auto result = solver.solve(initial);
    printResult<1>("Result", result.success, result.message, result.solution, result.finalWidth);

    if (result.success)
    {
        double x = result.solution[0].mid();
        std::cout << "验证: f(" << x << ") = " << (x*x*x - x - 1) << std::endl;
    }
}

void testExponential()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Test 3: e^x + y = 2, x = e^y" << std::endl;
    std::cout << "Expected Solution: (0.5, 0.5)" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    auto f = [](const IntervalVec<2>& x) -> IntervalVec<2>
    {
        IntervalVec<2> r;
        r[0] = exp(x[0]) + x[1] - Interval(2.0, 2.0);
        r[1] = x[0] - exp(x[1]);
        return r;
    };

    auto jacobian = [](const IntervalVec<2>& x) -> IntervalMat<2, 2>
    {
        IntervalMat<2, 2> J;
        J(0, 0) = exp(x[0]);
        J(0, 1) = Interval(1.0, 1.0);
        J(1, 0) = Interval(1.0, 1.0);
        J(1, 1) = -exp(x[1]);
        return J;
    };

    FixedKrawczykSolver<2> solver(f, jacobian, 1e-12, 100, 100);

    IntervalVec<2> initial;
    initial[0] = Interval(0.0, 1.0);
    initial[1] = Interval(0.0, 1.0);

    std::cout << "初始区间: " << initial << std::endl;
    auto result = solver.solve(initial);
    printResult<2>("结果", result.success, result.message, result.solution, result.finalWidth);

    if (result.success)
    {
        double x = result.solution[0].mid();
        double y = result.solution[1].mid();
        std::cout << "验证: f(" << x << ", " << y << ") = ("
                  << (std::exp(x) + y - 2) << ", " << (x - std::exp(y)) << ")" << std::endl;
    }
}

void testTrigonometric()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Test 4: sin(x) = 0.5" << std::endl;
    std::cout << "Expected Solution: x = π/6 ≈ 0.5236" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    auto f = [](const IntervalVec<1>& x) -> IntervalVec<1>
    {
        IntervalVec<1> r;
        r[0] = sin(x[0]) - Interval(0.5, 0.5);
        return r;
    };

    auto jacobian = [](const IntervalVec<1>& x) -> IntervalMat<1, 1>
    {
        IntervalMat<1, 1> J;
        J(0, 0) = cos(x[0]);
        return J;
    };

    FixedKrawczykSolver<1> solver(f, jacobian, 1e-12, 100, 100);

    IntervalVec<1> initial;
    initial[0] = Interval(0.0, 1.0);

    std::cout << "初始区间: " << initial << std::endl;
    auto result = solver.solve(initial);
    printResult<1>("结果", result.success, result.message, result.solution, result.finalWidth);

    if (result.success)
    {
        double x = result.solution[0].mid();
        std::cout << "验证: f(" << x << ") = " << (std::sin(x) - 0.5) << std::endl;
    }
}

void testCircleIntersection()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Test 5: x^2 + y = 2, x - y = 0" << std::endl;
    std::cout << "Expected Solution: (1, 1)" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    auto f = [](const IntervalVec<2>& x) -> IntervalVec<2>
    {
        IntervalVec<2> r;
        r[0] = x[0] * x[0] + x[1] - Interval(2.0, 2.0);
        r[1] = x[0] - x[1];
        return r;
    };

    auto jacobian = [](const IntervalVec<2>& x) -> IntervalMat<2, 2>
    {
        IntervalMat<2, 2> J;
        J(0, 0) = x[0] * Interval(2.0, 2.0);
        J(0, 1) = Interval(1.0, 1.0);
        J(1, 0) = Interval(1.0, 1.0);
        J(1, 1) = Interval(-1.0, -1.0);
        return J;
    };

    FixedKrawczykSolver<2> solver(f, jacobian, 1e-12, 100, 100);

    IntervalVec<2> initial;
    initial[0] = Interval(0.5, 1.5);
    initial[1] = Interval(0.5, 1.5);

    std::cout << "初始区间: " << initial << std::endl;
    auto result = solver.solve(initial);
    printResult<2>("结果", result.success, result.message, result.solution, result.finalWidth);

    if (result.success)
    {
        double x = result.solution[0].mid();
        double y = result.solution[1].mid();
        std::cout << "验证: f(" << x << ", " << y << ") = ("
                  << (x*x + y - 2) << ", " << (x - y) << ")" << std::endl;
    }
}

int main()
{
    std::cout << std::string(70, '#') << std::endl;
    std::cout << "#" << std::setw(68) << std::left << " Fixed Dimension Interval-Krawczyk Solver" << "#" << std::endl;
    std::cout << std::string(70, '#') << std::endl;

    testSimpleSystem();
    testCircleIntersection();
    testCubic();
    testTrigonometric();
    testExponential();

    std::cout << "\n" << std::string(70, '#') << std::endl;
    std::cout << "#" << std::setw(68) << std::left << " Testing Complete" << "#" << std::endl;
    std::cout << std::string(70, '#') << std::endl;

    return 0;
}
