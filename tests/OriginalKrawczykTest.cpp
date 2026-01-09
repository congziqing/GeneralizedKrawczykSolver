#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"
#include "KrawczykSolver.h"

void testSimpleSystem()
{
    std::cout << "=== Test: Simple System (Circle and Line Intersection) ===\n";
    
    using namespace krawczyk;
    
    KrawczykSolver solver(2, simpleFunction, simpleJacobian);
    
    IntervalVector initialBox(2);
    initialBox[0] = ik::KaucherInterval(-3.0, 3.0);
    initialBox[1] = ik::KaucherInterval(-3.0, 3.0);
    
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

int main()
{
    std::cout << "=== Original Krawczyk Method Demo ===\n\n";
    
    testSimpleSystem();
    
    std::cout << "=== Demo Complete ===\n";
    
    return 0;
}