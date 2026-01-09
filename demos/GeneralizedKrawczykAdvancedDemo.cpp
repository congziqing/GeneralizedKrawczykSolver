#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "../include/interval_krawczyk/PSGMDirectedInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
#include "../include/interval_krawczyk/GeneralizedKrawczykSolver.h"

// Auxiliary function: verify solution
template<size_t N>
void verifySolution(const std::string& desc, 
                   const ik::IntervalVector<N>& solution, 
                   const typename ik::GeneralizedKrawczykSolver<N>::Function& f)
{
    std::cout << "\nVerification: " << desc << std::endl;
    
    // Verify at midpoint
    ik::IntervalVector<N> c_vec;
    for (size_t i = 0; i < solution.Size; ++i)
    {
        double mid = solution[i].middle();
        c_vec[i] = PSGMDirectedInterval(mid, mid);
        std::cout << "  x[" << i << "] = " << mid << std::endl;
    }
    
    ik::IntervalVector<N> f_val = f(c_vec);
    std::cout << "  f(x) = " << f_val << std::endl;
    
    bool valid = true;
    for (size_t i = 0; i < solution.Size; ++i)
    {
        if (!f_val[i].contains(0.0))
        {
            valid = false;
            break;
        }
    }
    
    std::cout << "  Validity: " << (valid ? "✓ Valid" : "✗ Invalid") << std::endl;
}

// Test 1: Simple nonlinear system (proper intervals)
void testSimpleSystem()
{
    std::cout << "=== Test 1: Simple System (x^2 + y^2 = 4, x = y) ===\n";
    
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2> {
        ik::IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1] - PSGMDirectedInterval(4.0, 4.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2> {
        ik::IntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * PSGMDirectedInterval(2.0);
        J(0, 1) = x[1] * PSGMDirectedInterval(2.0);
        J(1, 0) = PSGMDirectedInterval(1.0);
        J(1, 1) = PSGMDirectedInterval(-1.0);
        return J;
    };
    
    ik::GeneralizedKrawczykSolver<2> solver(f, jacobian);
    
    // Use proper interval
    ik::IntervalVector<2> initialBox;
    initialBox[0] = PSGMDirectedInterval(0.5, 2.5);
    initialBox[1] = PSGMDirectedInterval(0.5, 2.5);
    
    std::cout << "Initial interval: " << initialBox << std::endl;
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution interval: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
        
        verifySolution("Solution validity", result.solution, f);
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

// Test 2: Cubic equation
void testCubicEquation()
{
    std::cout << "=== Test 2: Cubic Equation (x^3 - x - 1 = 0) ===\n";
    
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1> {
        ik::IntervalVector<1> result;
        result[0] = x[0] * x[0] * x[0] - x[0] - PSGMDirectedInterval(1.0, 1.0);
        return result;
    };

    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1> {
        ik::IntervalMatrix<1, 1> J;
        J(0, 0) = x[0] * x[0] * PSGMDirectedInterval(3.0) - PSGMDirectedInterval(1.0);
        return J;
    };
    
    ik::GeneralizedKrawczykSolver<1> solver(f, jacobian);
    
    ik::IntervalVector<1> initialBox;
    initialBox[0] = PSGMDirectedInterval(1.0, 2.0);
    
    std::cout << "Initial interval: " << initialBox << std::endl;
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution interval: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
        
        verifySolution("Solution validity", result.solution, f);
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

// Test 3: Quadratic system (demonstrating AE solution sets)
void testQuadraticSystem()
{
    std::cout << "=== Test 3: Quadratic System (x^2 + y = 2, x - y = 0) ===\n";
    
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2> {
        ik::IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] - PSGMDirectedInterval(2.0, 2.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2> {
        ik::IntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * PSGMDirectedInterval(2.0);
        J(0, 1) = PSGMDirectedInterval(1.0);
        J(1, 0) = PSGMDirectedInterval(1.0);
        J(1, 1) = PSGMDirectedInterval(-1.0);
        return J;
    };
    
    ik::GeneralizedKrawczykSolver<2> solver(f, jacobian);
    
    ik::IntervalVector<2> initialBox;
    initialBox[0] = PSGMDirectedInterval(0.5, 1.5);
    initialBox[1] = PSGMDirectedInterval(0.5, 1.5);
    
    std::cout << "Initial interval: " << initialBox << std::endl;
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution interval: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
        
        verifySolution("Solution validity", result.solution, f);
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

// Test 4: Using improper intervals (demonstrating Kaucher arithmetic advantages)
void testImproperInterval()
{
    std::cout << "=== Test 4: Using Improper Intervals ===\n";
    
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1> {
        ik::IntervalVector<1> result;
        result[0] = x[0] * x[0] - PSGMDirectedInterval(2.0, 2.0);
        return result;
    };

    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1> {
        ik::IntervalMatrix<1, 1> J;
        J(0, 0) = x[0] * PSGMDirectedInterval(2.0);
        return J;
    };
    
    ik::GeneralizedKrawczykSolver<1> solver(f, jacobian);
    
    // Use improper interval (demonstrating Kaucher arithmetic)
    ik::IntervalVector<1> initialBox;
    initialBox[0] = PSGMDirectedInterval(2.0, 1.0); // Improper interval
    
    std::cout << "Initial interval: " << initialBox << std::endl;
    std::cout << "Initial interval type: " << (initialBox[0].isProper() ? "Proper" : "Improper") << std::endl;
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution interval: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
        
        verifySolution("Solution validity", result.solution, f);
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

// Test 5: Trigonometric equation
void testTrigonometricEquation()
{
    std::cout << "=== Test 5: Trigonometric Equation (sin(x) - 0.5 = 0) ===\n";
    
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1> {
        ik::IntervalVector<1> result;
        result[0] = ::interval::sin(x[0]) - PSGMDirectedInterval(0.5, 0.5);
        return result;
    };

    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1> {
        ik::IntervalMatrix<1, 1> J;
        J(0, 0) = ::interval::cos(x[0]);
        return J;
    };
    
    ik::GeneralizedKrawczykSolver<1> solver(f, jacobian);
    
    ik::IntervalVector<1> initialBox;
    initialBox[0] = PSGMDirectedInterval(0.0, 1.0);
    
    std::cout << "Initial interval: " << initialBox << std::endl;
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution interval: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
        
        verifySolution("Solution validity", result.solution, f);
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

// Test 6: Find all solutions
void testAllSolutions()
{
    std::cout << "=== Test 6: Finding All Solutions ===\n";
    
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2> {
        ik::IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1] - PSGMDirectedInterval(4.0, 4.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2> {
        ik::IntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * PSGMDirectedInterval(2.0);
        J(0, 1) = x[1] * PSGMDirectedInterval(2.0);
        J(1, 0) = PSGMDirectedInterval(1.0);
        J(1, 1) = PSGMDirectedInterval(-1.0);
        return J;
    };
    
    ik::GeneralizedKrawczykSolver<2> solver(f, jacobian);
    
    // Use larger initial interval to find all solutions
    ik::IntervalVector<2> initialBox;
    initialBox[0] = PSGMDirectedInterval(-3.0, 3.0);
    initialBox[1] = PSGMDirectedInterval(-3.0, 3.0);
    
    std::cout << "Initial interval: " << initialBox << std::endl;
    
    auto solutions = solver.solveAll(initialBox, 10);
    
    std::cout << "Found " << solutions.size() << " solutions:\n";
    
    for (size_t i = 0; i < solutions.size(); ++i)
    {
        if (solutions[i].success)
        {
            std::cout << "\nSolution " << i + 1 << ": " << std::endl;
            std::cout << "  Interval: " << solutions[i].solution << std::endl;
            std::cout << "  Iterations: " << solutions[i].iterations << std::endl;
            std::cout << "  Width: " << solutions[i].finalWidth << std::endl;
        }
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Generalized Krawczyk Method Demo ===\n\n";
    
    // Show basic properties of Kaucher interval arithmetic
    std::cout << "=== Kaucher Interval Arithmetic Properties ===\n";
    PSGMDirectedInterval proper(1.0, 2.0);
    PSGMDirectedInterval improper(2.0, 1.0);
    
    std::cout << "Proper interval a = " << proper << std::endl;
    std::cout << "Improper interval b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "Inverse of a = " << (-proper) << std::endl;
    std::cout << "Dual of b = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // Run all tests
    testSimpleSystem();
    testQuadraticSystem();
    testCubicEquation();
    testImproperInterval();
    testTrigonometricEquation();
    testAllSolutions();
    
    std::cout << "=== Demo Completed ===\n";
    
    return 0;
}