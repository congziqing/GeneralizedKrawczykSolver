#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
using namespace ik;

// Helper function: create midpoint interval
KaucherInterval midPointInterval(double x)
{
    return KaucherInterval(x, x);
}

// Helper function: create dual interval (Kaucher dual)
KaucherInterval dualInterval(const KaucherInterval& interval)
{
    return KaucherInterval(interval.upper(), interval.lower());
}

// Test 1: AE-Solution Set Test - x² - p = 0, p ∈ [4, 9]
void testAESolutionSet()
{
    std::cout << "=== Test 1: AE-Solution Set Test (x² - p = 0, p ∈ [4, 9]) ===\n\n";
    
    // Define equation f(x, p) = x² - p = 0
    // Here p is a parameter, using dual form to represent universal quantifier
    auto f = [](const ik::IntervalVector<1>& x, const KaucherInterval& p) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        result[0] = x[0] * x[0] - p;
        return result;
    };
    
    // Calculate Jacobian with respect to x
    auto jacobian_x = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = 2x
        J(0, 0) = midPointInterval(2.0) * x[0];
        return J;
    };
    
    // Normal interval for parameter p: [4, 9]
    KaucherInterval p_normal(4.0, 9.0);
    // Dual interval for parameter p (representing universal quantifier): [9, 4]
    KaucherInterval p_dual = dualInterval(p_normal);
    
    std::cout << "Parameter p proper interval: " << p_normal << std::endl;
    std::cout << "Parameter p dual interval (universal quantifier): " << p_dual << std::endl;
    
    // Initial interval X: [0, 10]
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(0.0, 10.0);
    
    int maxIterations = 20;
    double tolerance = 1e-8;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // 1. Calculate midpoint
        double c = xn[0].middle();
        
        // 2. Calculate f at midpoint
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec, p_dual);
        
        // 3. Calculate Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian_x(xn);
        Eigen::MatrixXd J_mid = J.midpoint();
        
        // Check if Jacobian is invertible
        double det = J_mid(0, 0);
        if (std::abs(det) < 1e-12)
        {
            std::cout << "  Jacobian is singular, using approximate inverse...\n";
            det = (det > 0) ? 1e-12 : -1e-12;
        }
        Eigen::MatrixXd Y(1, 1);
        Y(0, 0) = 1.0 / det;
        
        // 4. Calculate term1: c - Y * f(c)
        double term1_val = c - Y(0, 0) * f_c[0].middle();
        ik::IntervalVector<1> term1;
        term1[0] = midPointInterval(term1_val);
        
        // 5. Calculate I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(1, 1);
        Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
        ik::IntervalMatrix<1, 1> I_minus_YJ_interval;
        I_minus_YJ_interval(0, 0) = midPointInterval(I_minus_YJ(0, 0));
        
        // 6. Calculate box - c
        double delta = (xn[0].upper() - xn[0].lower()) / 2.0;
        ik::IntervalVector<1> box_minus_c;
        box_minus_c[0] = KaucherInterval(-delta, delta);
        
        // 7. Calculate term3
        ik::IntervalVector<1> term3 = I_minus_YJ_interval * box_minus_c;
        
        // 8. Calculate K = term1 + term3
        ik::IntervalVector<1> K = term1 + term3;
        
        // 9. Calculate X_new = X ∩ K
        ik::IntervalVector<1> x_new = xn.intersection(K);
        
        std::cout << "Iteration " << iter + 1 << ":\n";
        std::cout << "  xn = " << xn << std::endl;
        std::cout << "  c = " << c << std::endl;
        std::cout << "  J = " << J(0, 0) << std::endl;
        std::cout << "  p_dual = " << p_dual << std::endl;
        std::cout << "  f(c) = " << f_c[0] << std::endl;
        std::cout << "  K = " << K << std::endl;
        std::cout << "  x_new = " << x_new << std::endl;
        std::cout << "  Width = " << x_new[0].width() << std::endl;
        
        // Check convergence
        if (x_new[0].width() < tolerance)
        {
            std::cout << "\n✓ Converged to " << x_new << ", Iterations: " << iter + 1 << std::endl;
            break;
        }
        
        xn = x_new;
        
        // If interval is empty or width increases, try using Kaucher arithmetic
        if (x_new[0].isEmpty() || x_new[0].width() > xn[0].width())
        {
            std::cout << "  Proper interval method encountered issues, using Kaucher arithmetic...\n";
            // Use Krawczyk interval as new interval, allowing improper intervals
            xn = K;
        }
        std::cout << std::endl;
    }
}

// Test 2: Complex AE-Solution Set Test - x - p + x³ = 0, p ∈ [-0.1, 0.1]
void testComplexAESolutionSet()
{
    std::cout << "\n=== Test 2: Complex AE-Solution Set Test (x - p + x³ = 0, p ∈ [-0.1, 0.1]) ===\n\n";
    
    // Define equation f(x, p) = x - p + x³ = 0
    // Here p is a parameter, using dual form to represent universal quantifier
    auto f = [](const ik::IntervalVector<1>& x, const KaucherInterval& p) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        result[0] = x[0] - p + x[0] * x[0] * x[0];
        return result;
    };
    
    // Calculate Jacobian with respect to x
    auto jacobian_x = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = 1 + 3x²
        J(0, 0) = midPointInterval(1.0) + midPointInterval(3.0) * x[0] * x[0];
        return J;
    };
    
    // Normal interval for parameter p: [-0.1, 0.1]
    KaucherInterval p_normal(-0.1, 0.1);
    // Dual interval for parameter p (representing universal quantifier): [0.1, -0.1]
    KaucherInterval p_dual = dualInterval(p_normal);
    
    std::cout << "Parameter p normal interval: " << p_normal << std::endl;
    std::cout << "Parameter p dual interval (universal quantifier): " << p_dual << std::endl;
    
    // Initial interval X: [-1, 1]
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(-1.0, 1.0);
    
    int maxIterations = 20;
    double tolerance = 1e-8;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // 1. Calculate midpoint
        double c = xn[0].middle();
        
        // 2. Calculate f at midpoint
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec, p_dual);
        
        // 3. Calculate Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian_x(xn);
        Eigen::MatrixXd J_mid = J.midpoint();
        
        // Check if Jacobian is invertible
        double det = J_mid(0, 0);
        if (std::abs(det) < 1e-12)
        {
            std::cout << "  Jacobian is singular, using approximate inverse...\n";
            det = (det > 0) ? 1e-12 : -1e-12;
        }
        Eigen::MatrixXd Y(1, 1);
        Y(0, 0) = 1.0 / det;
        
        // 4. Calculate term1: c - Y * f(c)
        double term1_val = c - Y(0, 0) * f_c[0].middle();
        ik::IntervalVector<1> term1;
        term1[0] = midPointInterval(term1_val);
        
        // 5. Calculate I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(1, 1);
        Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
        ik::IntervalMatrix<1, 1> I_minus_YJ_interval;
        I_minus_YJ_interval(0, 0) = midPointInterval(I_minus_YJ(0, 0));
        
        // 6. Calculate box - c
        double delta = (xn[0].upper() - xn[0].lower()) / 2.0;
        ik::IntervalVector<1> box_minus_c;
        box_minus_c[0] = KaucherInterval(-delta, delta);
        
        // 7. Calculate term3
        ik::IntervalVector<1> term3 = I_minus_YJ_interval * box_minus_c;
        
        // 8. Calculate K = term1 + term3
        ik::IntervalVector<1> K = term1 + term3;
        
        // 9. Calculate X_new = X ∩ K
        ik::IntervalVector<1> x_new = xn.intersection(K);
        
        std::cout << "Iteration " << iter + 1 << ":\n";
        std::cout << "  xn = " << xn << std::endl;
        std::cout << "  c = " << c << std::endl;
        std::cout << "  J = " << J(0, 0) << std::endl;
        std::cout << "  p_dual = " << p_dual << std::endl;
        std::cout << "  f(c) = " << f_c[0] << std::endl;
        std::cout << "  K = " << K << std::endl;
        std::cout << "  x_new = " << x_new << std::endl;
        std::cout << "  Width = " << x_new[0].width() << std::endl;
        
        // Check convergence
        if (x_new[0].width() < tolerance)
        {
            std::cout << "\n✓ Converged to " << x_new << ", Iterations: " << iter + 1 << std::endl;
            break;
        }
        
        xn = x_new;
        
        // If interval is empty or width increases, try using Kaucher arithmetic
        if (x_new[0].isEmpty() || x_new[0].width() > xn[0].width())
        {
            std::cout << "  Proper interval method encountered issues, using Kaucher arithmetic...\n";
            // Use Krawczyk interval as new interval, allowing improper intervals
            xn = K;
        }
        std::cout << std::endl;
    }
}

int main()
{
    std::cout << "=== Kaucher AE-Solution Set Test ===\n\n";
    
    // Run tests
    testAESolutionSet();
    testComplexAESolutionSet();
    
    std::cout << "\n=== Tests Completed ===\n";
    return 0;
}
