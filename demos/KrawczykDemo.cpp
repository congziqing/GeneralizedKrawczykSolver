#include <iostream>
#include <iomanip>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"
#include "KrawczykSolver.h"

void printResult(const KrawczykResult& result)
{
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Solution Result:" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    if (result.success)
    {
        std::cout << "✓ Successfully found solution!" << std::endl;
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Final interval width: " << std::scientific << std::setprecision(6)
                  << result.finalWidth << std::endl;
        std::cout << "\nSolution Interval:" << std::endl;

        for (size_t i = 0; i < result.solution.size(); ++i)
        {
            std::cout << "  x[" << i << "] = [" << std::fixed << std::setprecision(10)
                      << result.solution[i].lower() << ", "
                      << result.solution[i].upper() << "]" << std::endl;
            std::cout << "      Midpoint: " << result.solution[i].middle() << std::endl;
        }

        double mid = result.solution[0].middle();
        if (result.solution.size() > 1)
            mid = result.solution[1].middle();

        std::cout << "\nMessage: " << result.message << std::endl;
    }
    else
    {
        std::cout << "✗ Solution failed" << std::endl;
        std::cout << "Status code: " << static_cast<int>(result.status) << std::endl;
        std::cout << "Message: " << result.message << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
}

void example1_simpleSystem()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Example 1: Simple Nonlinear System" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "Equations: " << std::endl;
    std::cout << "  f1(x,y) = x^2 + y^2 - 4 = 0" << std::endl;
    std::cout << "  f2(x,y) = x - y = 0" << std::endl;
    std::cout << "Expected solutions: x = ±√2, y = ±√2" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = x[0] * x[0] + x[1] * x[1] - ik::KaucherInterval(4.0, 4.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = x[0] * ik::KaucherInterval(2.0);
        J[0][1] = x[1] * ik::KaucherInterval(2.0);
        J[1][0] = ik::KaucherInterval(1.0);
        J[1][1] = ik::KaucherInterval(-1.0);
        return J;
    };

    std::cout << "\nTest 1.1: Initial interval [-3, 3] × [-3, 3]" << std::endl;
    IntervalVector initial1(2);
    initial1[0] = ik::KaucherInterval(-3.0, 3.0);
    initial1[1] = ik::KaucherInterval(-3.0, 3.0);

    KrawczykSolver solver(2, f, jacobian, 1e-8, 50, 30);
    KrawczykResult result1 = solver.solve(initial1);
    printResult(result1);

    std::cout << "\nTest 1.2: Initial interval [0, 3] × [0, 3]" << std::endl;
    IntervalVector initial2(2);
    initial2[0] = ik::KaucherInterval(0.0, 3.0);
    initial2[1] = ik::KaucherInterval(0.0, 3.0);

    KrawczykResult result2 = solver.solve(initial2);
    printResult(result2);
}

void example2_polynomial()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Example 2: Polynomial System" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "Equations: " << std::endl;
    std::cout << "  f1(x,y) = x^2 + y - 3 = 0" << std::endl;
    std::cout << "  f2(x,y) = x*y - 1 = 0" << std::endl;
    std::cout << "Expected solutions: (1,2) and (-1,-2)" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = x[0] * x[0] + x[1] - ik::KaucherInterval(3.0, 3.0);
        result[1] = x[0] * x[1] - ik::KaucherInterval(1.0, 1.0);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = x[0] * ik::KaucherInterval(2.0);
        J[0][1] = ik::KaucherInterval(1.0);
        J[1][0] = x[1];
        J[1][1] = x[0];
        return J;
    };

    std::cout << "\nTest 2.1: Initial interval [-3, 3] × [-3, 3]" << std::endl;
    IntervalVector initial(2);
    initial[0] = ik::KaucherInterval(-3.0, 3.0);
    initial[1] = ik::KaucherInterval(-3.0, 3.0);

    KrawczykSolver solver(2, f, jacobian, 1e-8, 50, 30);
    KrawczykResult result = solver.solve(initial);
    printResult(result);

    std::cout << "\nTest 2.2: Search all solutions" << std::endl;
    auto allResults = solver.solveAll(initial, 5);
    std::cout << "Found " << allResults.size() << " solutions" << std::endl;
    for (size_t i = 0; i < allResults.size(); ++i)
    {
        std::cout << "\nSolution #" << (i + 1) << ":" << std::endl;
        for (size_t j = 0; j < allResults[i].solution.size(); ++j)
        {
            std::cout << "  x[" << j << "] = [" << std::fixed << std::setprecision(8)
                      << allResults[i].solution[j].lower() << ", "
                      << allResults[i].solution[j].upper() << "]" << std::endl;
        }
    }
}

void example3_cubic()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Example 3: Cubic Equation" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "Equation: x^3 - x - 1 = 0" << std::endl;
    std::cout << "Expected solution: x ≈ 1.3247" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(1);
        result[0] = x[0] * x[0] * x[0] - x[0] - ik::KaucherInterval(1.0, 1.0);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(1, 1);
        J[0][0] = x[0] * x[0] * ik::KaucherInterval(3.0) - ik::KaucherInterval(1.0);
        return J;
    };

    std::cout << "\nTest 3.1: Initial interval [0, 2]" << std::endl;
    IntervalVector initial(1);
    initial[0] = ik::KaucherInterval(0.0, 2.0);

    KrawczykSolver solver(1, f, jacobian, 1e-10, 100, 20);
    KrawczykResult result = solver.solve(initial);
    printResult(result);
}

void example4_trigonometric()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Example 4: Trigonometric Equation" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "Equation: sin(x) - x/2 = 0" << std::endl;
    std::cout << "Expected solution: x ≈ 0" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(1);
        result[0] = interval::sin(x[0]) - x[0] * ik::KaucherInterval(0.5);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(1, 1);
        J[0][0] = interval::cos(x[0]) - ik::KaucherInterval(0.5);
        return J;
    };

    std::cout << "\nTest 4.1: Initial interval [-π, π]" << std::endl;
    IntervalVector initial(1);
    initial[0] = ik::KaucherInterval(-3.14159, 3.14159);

    KrawczykSolver solver(1, f, jacobian, 1e-10, 50, 20);
    KrawczykResult result = solver.solve(initial);
    printResult(result);
}

void example5_exponential()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Example 5: Exponential System" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "Equations: " << std::endl;
    std::cout << "  f1(x,y) = e^x + y - 2 = 0" << std::endl;
    std::cout << "  f2(x,y) = x - e^y = 0" << std::endl;
    std::cout << "Expected solution: x ≈ 0.5, y ≈ 0.5" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = interval::exp(x[0]) + x[1] - ik::KaucherInterval(2.0, 2.0);
        result[1] = x[0] - interval::exp(x[1]);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = interval::exp(x[0]);
        J[0][1] = ik::KaucherInterval(1.0);
        J[1][0] = ik::KaucherInterval(1.0);
        J[1][1] = -interval::exp(x[1]);
        return J;
    };

    std::cout << "\nTest 5.1: Initial interval [0, 1] × [0, 1]" << std::endl;
    IntervalVector initial(2);
    initial[0] = ik::KaucherInterval(0.0, 1.0);
    initial[1] = ik::KaucherInterval(0.0, 1.0);

    KrawczykSolver solver(2, f, jacobian, 1e-10, 50, 20);
    KrawczykResult result = solver.solve(initial);
    printResult(result);
}

void demonstrateIntervalArithmetic()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "Interval Arithmetic Demonstration" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    ik::KaucherInterval a(1.0, 3.0);
    ik::KaucherInterval b(2.0, 4.0);

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;

    std::cout << "\nBasic Operations:" << std::endl;
    std::cout << "a + b = " << (a + b) << std::endl;
    std::cout << "a - b = " << (a - b) << std::endl;
    std::cout << "a * b = " << (a * b) << std::endl;
    std::cout << "a / b = " << (a / b) << std::endl;

    std::cout << "\nCommon Functions:" << std::endl;
    ik::KaucherInterval x(0.5, 2.0);
    std::cout << "x = " << x << std::endl;
    std::cout << "x² = " << (x * x) << std::endl;
    std::cout << "√x = " << interval::sqrt(x) << std::endl;
    std::cout << "e^x = " << interval::exp(x) << std::endl;
    std::cout << "sin(x) = " << interval::sin(x) << std::endl;
    std::cout << "cos(x) = " << interval::cos(x) << std::endl;
}

int main()
{
    std::cout << std::string(70, '#') << std::endl;
    std::cout << "#" << std::setw(68) << std::left << " Interval-Krawczyk Method for Nonlinear Equations" << "#" << std::endl;
    std::cout << std::string(70, '#') << std::endl;

    demonstrateIntervalArithmetic();

    try
    {
        example1_simpleSystem();
        example2_polynomial();
        example3_cubic();
        example4_trigonometric();
        example5_exponential();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n" << std::string(70, '#') << std::endl;
    std::cout << "#" << std::setw(68) << std::left << " All Examples Completed" << "#" << std::endl;
    std::cout << std::string(70, '#') << std::endl;

    return 0;
}
