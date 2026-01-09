#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"

// Newton method pre-iteration function to find high-precision center point
auto newton_method_2d(double x1_init, double x2_init)
{
    // Define function f(x1, x2) = (x1^2 + x2^2 - 4, x1 - x2)
    auto f = [](double x1, double x2) -> std::array<double, 2> {
        return {
            x1 * x1 + x2 * x2 - 4.0,
            x1 - x2
        };
    };
    
    // Define Jacobian matrix point function
    auto df = [](double x1, double x2) -> std::array<std::array<double, 2>, 2> {
        std::array<std::array<double, 2>, 2> result;
        result[0][0] = 2.0 * x1;
        result[0][1] = 2.0 * x2;
        result[1][0] = 1.0;
        result[1][1] = -1.0;
        return result;
    };
    
    // Newton iteration
    double x1 = x1_init;
    double x2 = x2_init;
    int max_iter = 10;
    double tol = 1e-12;
    
    for (int iter = 0; iter < max_iter; ++iter)
    {
        auto fx = f(x1, x2);
        double fx1 = fx[0];
        double fx2 = fx[1];
        
        if (std::abs(fx1) < tol && std::abs(fx2) < tol)
        {
            break;
        }
        
        auto J = df(x1, x2);
        double J11 = J[0][0];
        double J12 = J[0][1];
        double J21 = J[1][0];
        double J22 = J[1][1];
        
        double det = J11 * J22 - J12 * J21;
        if (std::abs(det) < 1e-12)
        {
            break;
        }
        
        double invJ11 = J22 / det;
        double invJ12 = -J12 / det;
        double invJ21 = -J21 / det;
        double invJ22 = J11 / det;
        
        double delta1 = -(invJ11 * fx1 + invJ12 * fx2);
        double delta2 = -(invJ21 * fx1 + invJ22 * fx2);
        
        x1 += delta1;
        x2 += delta2;
        
        if (std::abs(delta1) < tol && std::abs(delta2) < tol)
        {
            break;
        }
    }
    
    return std::make_pair(x1, x2);
}

int main()
{
    std::cout << "=== Tolerance Embedding Simple Implementation ===\n\n";
    
    // 1. Use Newton method pre-iteration to find high-precision center point x*
    std::cout << "1. Newton method pre-iteration to find center point...\n";
    double x1_init = 1.5;
    double x2_init = 1.5;
    auto [x1_star, x2_star] = newton_method_2d(x1_init, x2_init);
    
    std::cout << "   Initial guess: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   Newton method result: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // 2. Construct small initial interval
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.01; // Use smaller initial interval radius
    FixedIntervalVector<2> X;
    X[0] = ik::KaucherInterval(x1_star - delta, x1_star + delta);
    X[1] = ik::KaucherInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   Initial interval X: " << X << std::endl;
    
    // 3. Define function f
    auto f = [](const FixedIntervalVector<2>& x) -> FixedIntervalVector<2> {
        FixedIntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
    };
    
    // 4. Calculate f(X)
    std::cout << "\n3. Calculating f(X)...\n";
    FixedIntervalVector<2> fX = f(X);
    
    std::cout << "   f(X) = " << fX << std::endl;
    
    // 5. Define target Y (improper interval)
    std::cout << "\n4. Defining target Y (improper interval)...\n";
    FixedIntervalVector<2> Y;
    Y[0] = ik::KaucherInterval(3.9, 4.1).dual(); // [4.1, 3.9]
    Y[1] = ik::KaucherInterval(-0.1, 0.1).dual(); // [0.1, -0.1]
    
    std::cout << "   Y = " << Y << std::endl;
    
    // 6. Verify inner enclosure
    std::cout << "\n5. Verifying inner enclosure...\n";
    
    // Check if fX is strictly contained within the dual of Y
    bool inner_inclusion = true;
    for (size_t i = 0; i < 2; ++i)
    {
        ik::KaucherInterval Y_dual = Y[i].dual(); // Convert to proper interval
        std::cout << "   Y_dual[" << i << "] = " << Y_dual << std::endl;
        std::cout << "   fX[" << i << "] = " << fX[i] << std::endl;
        
        bool contains = Y_dual.contains(fX[i]);
        std::cout << "   fX[" << i << "] contained in Y_dual[" << i << "]: " << (contains ? "Yes" : "No") << std::endl;
        
        if (!contains)
        {
            inner_inclusion = false;
        }
    }
    
    // 7. Final result
    std::cout << "\n6. Final result: " << (inner_inclusion ? "✓ Inner enclosure verified" : "✗ Inner enclosure failed") << std::endl;
    
    return 0;
}