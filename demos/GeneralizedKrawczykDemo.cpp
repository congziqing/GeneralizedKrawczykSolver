#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"
#include "GeneralizedKrawczykSolver.h"

void testSimpleSystem()
{
    std::cout << "=== Test: Simple System (Circle and Line Intersection) ===\n";
    
    using namespace generalized_krawczyk;
    
    GeneralizedKrawczykSolver solver(2, simpleFunction, simpleJacobian, 1e-8, 100, 100);
    
    IntervalVector initialBox(2);
    initialBox[0] = ik::KaucherInterval(0.5, 2.5);
    initialBox[1] = ik::KaucherInterval(0.5, 2.5);
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

void testCubicEquation()
{
    std::cout << "=== Test: Cubic Equation (x^3 - x - 1 = 0) ===\n";
    
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
    
    GeneralizedKrawczykSolver solver(1, f, jacobian, 1e-8, 100, 100);
    
    IntervalVector initialBox(1);
    initialBox[0] = ik::KaucherInterval(1.0, 2.0);
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

void testQuadraticEquation()
{
    std::cout << "=== Test: Quadratic Equation (x^2 + y = 2, x - y = 0) ===\n";
    
    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = x[0] * x[0] + x[1] - ik::KaucherInterval(2.0, 2.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = x[0] * ik::KaucherInterval(2.0);
        J[0][1] = ik::KaucherInterval(1.0);
        J[1][0] = ik::KaucherInterval(1.0);
        J[1][1] = ik::KaucherInterval(-1.0);
        return J;
    };
    
    GeneralizedKrawczykSolver solver(2, f, jacobian, 1e-8, 100, 100);
    
    IntervalVector initialBox(2);
    initialBox[0] = ik::KaucherInterval(0.5, 1.5);
    initialBox[1] = ik::KaucherInterval(0.5, 1.5);
    
    auto result = solver.solve(initialBox);
    
    if (result.success)
    {
        std::cout << "✓ Success!\n";
        std::cout << "Iterations: " << result.iterations << std::endl;
        std::cout << "Solution: " << result.solution << std::endl;
        std::cout << "Final width: " << result.finalWidth << std::endl;
    }
    else
    {
        std::cout << "✗ Failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

void testAllSolutions()
{
    std::cout << "=== Test: Finding All Solutions ===\n";
    
    using namespace generalized_krawczyk;
    
    GeneralizedKrawczykSolver solver(2, simpleFunction, simpleJacobian);
    
    IntervalVector initialBox(2);
    initialBox[0] = ik::KaucherInterval(0.5, 2.5);
    initialBox[1] = ik::KaucherInterval(0.5, 2.5);
    
    auto solutions = solver.solveAll(initialBox);
    
    std::cout << "Found " << solutions.size() << " solutions:\n";
    
    for (size_t i = 0; i < solutions.size(); ++i)
    {
        if (solutions[i].success)
        {
            std::cout << "Solution " << i + 1 << ": " << solutions[i].solution << std::endl;
            std::cout << "  Iterations: " << solutions[i].iterations << std::endl;
            std::cout << "  Width: " << solutions[i].finalWidth << std::endl;
        }
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Generalized Krawczyk Method Demo ===\n\n";
    
    testCubicEquation();
    testQuadraticEquation();
    testSimpleSystem();
    testAllSolutions();
    
    std::cout << "=== Demo Complete ===\n";
    
    return 0;
}