#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"

using namespace ik;

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

// Fixed-dimension Generalized Krawczyk Solver
template<size_t N>
class FixedGeneralizedKrawczykSolver
{
public:
    using Function = std::function<IntervalVector<N>(const IntervalVector<N>&)>;
    using JacobianFunction = std::function<IntervalMatrix<N, N>(const IntervalVector<N>&)>;
    
    FixedGeneralizedKrawczykSolver(Function f, JacobianFunction jacobian, 
                                  double tolerance = 1e-8, int maxIterations = 50)
        : f_(f), jacobian_(jacobian), tolerance_(tolerance), maxIterations_(maxIterations) {}
    
    struct Result
    {
        bool success;
        IntervalVector<N> solution;
        int iterations;
        double finalWidth;
        std::string message;
    };
    
    Result solve(const IntervalVector<N>& initialBox) const
    {
        IntervalVector<N> current = initialBox;
        
        for (int iter = 0; iter < maxIterations_; ++iter)
        {
            double width = current.maxWidth();
            if (width < tolerance_)
            {
                return {true, current, iter, width, "Converged successfully"};
            }
            
            // Calculate midpoint c
            Eigen::VectorXd c = current.middle();
            
            // Create interval vector for midpoint
            IntervalVector<N> c_vec;
            for (size_t i = 0; i < N; ++i)
            {
                c_vec[i] = KaucherInterval(c(i), c(i));
            }
            
            // Calculate f(c)
            IntervalVector<N> f_c = f_(c_vec);
            
            // Calculate Jacobian matrix J(x)
            IntervalMatrix<N, N> J = jacobian_(current);
            
            // Calculate midpoint of Jacobian matrix
            Eigen::MatrixXd J_mid = J.middle();
            
            // Check if Jacobian is singular
            double det = J_mid.determinant();
            if (std::abs(det) < 1e-12)
            {
                return {false, current, iter, width, "Singular Jacobian matrix"};
            }
            
            // Calculate inverse of Jacobian matrix
            Eigen::MatrixXd Y = J_mid.inverse();
            
            // Calculate term1: c - Y * f(c)
            IntervalVector<N> term1;
            for (size_t i = 0; i < N; ++i)
            {
                double val = c(i);
                for (size_t j = 0; j < N; ++j)
                {
                    val -= Y(static_cast<int>(i), static_cast<int>(j)) * f_c[j].middle();
                }
                term1[i] = KaucherInterval(val, val);
            }
            
            // Calculate I - Y * J(x)
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(static_cast<int>(N), static_cast<int>(N));
            Eigen::MatrixXd I_minus_YJ = I - Y * J.middle();
            
            // Create interval matrix I_minus_YJ_interval
            IntervalMatrix<N, N> I_minus_YJ_interval;
            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = 0; j < N; ++j)
                {
                    double val = I_minus_YJ(static_cast<int>(i), static_cast<int>(j));
                    I_minus_YJ_interval(i, j) = KaucherInterval(val, val);
                }
            }
            
            // Calculate box - c
            IntervalVector<N> box_minus_c;
            for (size_t i = 0; i < N; ++i)
            {
                double delta = (current[i].upper() - current[i].lower()) / 2.0;
                box_minus_c[i] = KaucherInterval(-delta, delta);
            }
            
            // Calculate term3: (I - YJ) * (box - c)
            IntervalVector<N> term3 = I_minus_YJ_interval * box_minus_c;
            
            // Calculate Krawczyk interval K = term1 + term3
            IntervalVector<N> K = term1 + term3;
            
            // Calculate intersection K ∩ current
            IntervalVector<N> next;
            for (size_t i = 0; i < N; ++i)
            {
                next[i] = current[i].meet(K[i]);
            }
            
            // Check if converged or diverged
            if (next[0].isEmpty() || next.maxWidth() > current.maxWidth() * 1000 || next.maxWidth() > 1e10)
            {
                return {false, current, iter, width, "Krawczyk iteration produced empty or diverging interval"};
            }
            
            current = next;
        }
        
        return {false, current, maxIterations_, current.maxWidth(), "Maximum iterations reached"};
    }
    
private:
    Function f_;
    JacobianFunction jacobian_;
    double tolerance_;
    int maxIterations_;
};

// Verification function: Check if f(X) is strictly contained within Y's dual
template<size_t N>
bool verify_inner_inclusion(const IntervalVector<N>& X, 
                           const IntervalVector<N>& Y, 
                           const typename FixedGeneralizedKrawczykSolver<N>::Function& f)
{
    // Calculate f(X)
    IntervalVector<N> fX = f(X);
    
    std::cout << "   f(X) = " << fX << std::endl;
    std::cout << "   Y    = " << Y << std::endl;
    
    // Check if fX is strictly contained within Y's dual interval
    // For improper Y intervals, the dual is proper [Y.upper(), Y.lower()]
    for (size_t i = 0; i < N; ++i)
    {
        if (!Y[i].isImproper())
        {
            std::cout << "   Error: Y[" << i << "] is not an improper interval" << std::endl;
            return false;
        }
        
        KaucherInterval Y_dual = Y[i].dual(); // Convert to proper interval [Y.upper(), Y.lower()]
        
        if (!fX[i].isProper())
        {
            std::cout << "   Error: f(X)[" << i << "] is not a proper interval" << std::endl;
            return false;
        }
        
        // Check if fX(i) is strictly contained within Y_dual
        bool contains = Y_dual.contains(fX[i]);
        std::cout << "   Is f(X)[" << i << "] contained in Y_dual[" << i << "]: " << (contains ? "Yes" : "No") << std::endl;
        
        if (!contains)
        {
            return false;
        }
    }
    
    return true;
}

int main()
{
    std::cout << "=== Tolerance Embedding Final Solution ===\n\n";
    
    // 1. Use Newton method to find high-precision center point x*
    std::cout << "1. Newton method pre-iteration to find center point...\n";
    double x1_init = 1.5;
    double x2_init = 1.5;
    auto [x1_star, x2_star] = newton_method_2d(x1_init, x2_init);
    
    std::cout << "   Initial guess: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   Newton method result: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // 2. Construct small initial interval
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.01; // Use smaller initial interval radius
    IntervalVector<2> X0;
    X0[0] = KaucherInterval(x1_star - delta, x1_star + delta);
    X0[1] = KaucherInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   Initial interval X0: " << X0 << std::endl;
    
    // 3. Define function f and Jacobian matrix
    std::cout << "\n3. Defining function f and Jacobian matrix...\n";
    
    // Define original function f
    auto f = [](const IntervalVector<2>& x) -> IntervalVector<2> {
        IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
    };
    
    // Define Jacobian matrix
    auto jacobian = [](const IntervalVector<2>& x) -> IntervalMatrix<2, 2> {
        IntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * KaucherInterval(2.0);
        J(0, 1) = x[1] * KaucherInterval(2.0);
        J(1, 0) = KaucherInterval(1.0);
        J(1, 1) = KaucherInterval(-1.0);
        return J;
    };
    
    // 4. Run generalized Krawczyk iteration
    std::cout << "\n4. Running generalized Krawczyk iteration...\n";
    FixedGeneralizedKrawczykSolver<2> solver(f, jacobian, 1e-8, 50);
    
    auto result = solver.solve(X0);
    
    if (result.success)
    {
        std::cout << "   ✓ Krawczyk iteration succeeded!\n";
        std::cout << "   Iterations: " << result.iterations << std::endl;
        std::cout << "   Solution interval X_final: " << result.solution << std::endl;
        std::cout << "   Final width: " << result.finalWidth << std::endl;
        
        // 5. Define target Y (improper interval)
        std::cout << "\n5. Defining target Y (improper interval)...\n";
        IntervalVector<2> Y;
        Y[0] = KaucherInterval(3.9, 4.1).dual(); // Improper interval [4.1, 3.9]
        Y[1] = KaucherInterval(-0.1, 0.1).dual(); // Improper interval [0.1, -0.1]
        
        std::cout << "   Y = " << Y << std::endl;
        
        // 6. Inner enclosure verification
        std::cout << "\n6. Inner enclosure verification...\n";
        bool inner_inclusion = verify_inner_inclusion(result.solution, Y, f);
        
        std::cout << "\n7. Final result: " << (inner_inclusion ? "✓ Inner enclosure verified" : "✗ Inner enclosure failed") << std::endl;
        
        // 7. Additional verification: Check four corner points
        std::cout << "\n8. Corner point verification...\n";
        
        // Define point function f
        auto point_f = [](double x1, double x2) -> std::array<double, 2> {
            return {
                x1 * x1 + x2 * x2,
                x1 - x2
            };
        };
        
        std::array<std::pair<double, double>, 4> corners;
        corners[0] = std::make_pair(result.solution[0].lower(), result.solution[1].lower());
        corners[1] = std::make_pair(result.solution[0].lower(), result.solution[1].upper());
        corners[2] = std::make_pair(result.solution[0].upper(), result.solution[1].lower());
        corners[3] = std::make_pair(result.solution[0].upper(), result.solution[1].upper());
        
        for (size_t i = 0; i < corners.size(); ++i)
        {
            auto [x1, x2] = corners[i];
            auto fx = point_f(x1, x2);
            
            // Check if within target Y (considering tolerance)
            bool in_y1 = (fx[0] >= 3.9 && fx[0] <= 4.1);
            bool in_y2 = (fx[1] >= -0.1 && fx[1] <= 0.1);
            
            std::cout << "   Corner " << i+1 << " (" << x1 << ", " << x2 << "): "
                      << "f = (" << fx[0] << ", " << fx[1] << ") - "
                      << (in_y1 && in_y2 ? "✓ Valid" : "✗ Invalid") << std::endl;
        }
    }
    else
    {
        std::cout << "   ✗ Krawczyk iteration failed!\n";
        std::cout << "   Message: " << result.message << std::endl;
        std::cout << "   Final interval: " << result.solution << std::endl;
        std::cout << "   Final width: " << result.finalWidth << std::endl;
    }
    
    std::cout << "\n=== Test Completed ===\n";
    
    return 0;
}