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
double newton_method_1d(double x_init, const std::function<double(double)>& f, const std::function<double(double)>& df, int max_iter = 10, double tol = 1e-12)
{
    double x = x_init;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        double fx = f(x);
        if (std::abs(fx) < tol)
        {
            break;
        }
        double dfx = df(x);
        if (std::abs(dfx) < 1e-12)
        {
            break;
        }
        double delta = -fx / dfx;
        x += delta;
        if (std::abs(delta) < tol)
        {
            break;
        }
    }
    return x;
}

// 2D Newton method pre-iteration function
auto newton_method_2d(double x1_init, double x2_init, 
                      const std::function<std::array<double, 2>(double, double)>& f, 
                      const std::function<std::array<std::array<double, 2>, 2>(double, double)>& df, 
                      int max_iter = 10, double tol = 1e-12)
{
    double x1 = x1_init;
    double x2 = x2_init;
    
    for (int iter = 0; iter < max_iter; ++iter)
    {
        // Calculate function values
        auto fx = f(x1, x2);
        double fx1 = fx[0];
        double fx2 = fx[1];
        
        // Check convergence
        if (std::abs(fx1) < tol && std::abs(fx2) < tol)
        {
            break;
        }
        
        // Calculate Jacobian matrix
        auto J = df(x1, x2);
        double J11 = J[0][0];
        double J12 = J[0][1];
        double J21 = J[1][0];
        double J22 = J[1][1];
        
        // Calculate determinant
        double det = J11 * J22 - J12 * J21;
        if (std::abs(det) < 1e-12)
        {
            break;
        }
        
        // Calculate inverse Jacobian matrix
        double invJ11 = J22 / det;
        double invJ12 = -J12 / det;
        double invJ21 = -J21 / det;
        double invJ22 = J11 / det;
        
        // Update x
        double delta1 = -(invJ11 * fx1 + invJ12 * fx2);
        double delta2 = -(invJ21 * fx1 + invJ22 * fx2);
        
        x1 += delta1;
        x2 += delta2;
        
        // Check step convergence
        if (std::abs(delta1) < tol && std::abs(delta2) < tol)
        {
            break;
        }
    }
    
    return std::make_pair(x1, x2);
}

// Fixed-dimensional generalized Krawczyk solver
template<size_t N>
class FixedGeneralizedKrawczykSolver
{
public:
    using Function = std::function<FixedIntervalVector<N>(const FixedIntervalVector<N>&)>;
    using JacobianFunction = std::function<FixedIntervalMatrix<N, N>(const FixedIntervalVector<N>&)>;
    
    FixedGeneralizedKrawczykSolver(Function f, JacobianFunction jacobian, 
                                  double tolerance = 1e-8, int maxIterations = 50)
        : f_(f), jacobian_(jacobian), tolerance_(tolerance), maxIterations_(maxIterations) {}
    
    struct Result
    {
        bool success;
        FixedIntervalVector<N> solution;
        int iterations;
        double finalWidth;
        std::string message;
    };
    
    Result solve(const FixedIntervalVector<N>& initialBox) const
    {
        FixedIntervalVector<N> current = initialBox;
        
        for (int iter = 0; iter < maxIterations_; ++iter)
        {
            double width = current.width();
            if (width < tolerance_)
            {
                return {true, current, iter, width, "Converged successfully"};
            }
            
            // Calculate midpoint c
            Eigen::VectorXd c = current.midpoint();
            
            // Create interval vector for midpoint
            FixedIntervalVector<N> c_vec;
            for (size_t i = 0; i < N; ++i)
            {
                c_vec[i] = ik::KaucherInterval(c(i), c(i));
            }
            
            // Calculate f(c)
            FixedIntervalVector<N> f_c = f_(c_vec);
            
            // Calculate Jacobian matrix J(x)
            FixedIntervalMatrix<N, N> J = jacobian_(current);
            
            // Calculate midpoint of Jacobian matrix
            Eigen::MatrixXd J_mid = J.midpoint();
            
            // Check if Jacobian matrix is singular
            double det = J_mid.determinant();
            if (std::abs(det) < 1e-12)
            {
                return {false, current, iter, width, "Singular Jacobian matrix"};
            }
            
            // Calculate inverse of Jacobian matrix
            Eigen::MatrixXd Y = J_mid.inverse();
            
            // Calculate term1: c - Y * f(c)
            FixedIntervalVector<N> term1;
            for (size_t i = 0; i < N; ++i)
            {
                double val = c(i);
                for (size_t j = 0; j < N; ++j)
                {
                    val -= Y(static_cast<int>(i), static_cast<int>(j)) * f_c[j].middle();
                }
                term1[i] = ik::KaucherInterval(val, val);
            }
            
            // Calculate I - Y * J(x)
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(static_cast<int>(N), static_cast<int>(N));
            Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
            
            // Create interval matrix I_minus_YJ_interval
            FixedIntervalMatrix<N, N> I_minus_YJ_interval;
            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = 0; j < N; ++j)
                {
                    double val = I_minus_YJ(static_cast<int>(i), static_cast<int>(j));
                    I_minus_YJ_interval(i, j) = ik::KaucherInterval(val, val);
                }
            }
            
            // Calculate box - c
            FixedIntervalVector<N> box_minus_c;
            for (size_t i = 0; i < N; ++i)
            {
                double delta = (current[i].upper() - current[i].lower()) / 2.0;
                box_minus_c[i] = ik::KaucherInterval(-delta, delta);
            }
            
            // Calculate term3: (I - YJ) * (box - c)
            FixedIntervalVector<N> term3 = I_minus_YJ_interval * box_minus_c;
            
            // Calculate Krawczyk interval K = term1 + term3
            FixedIntervalVector<N> K = term1 + term3;
            
            // Calculate intersection K ∩ current
            FixedIntervalVector<N> next;
            for (size_t i = 0; i < N; ++i)
            {
                next[i] = current[i].meet(K[i]);
            }
            
            // Check for convergence or divergence
            if (next[0].isEmpty() || next.width() > current.width() * 1000 || next.width() > 1e10)
            {
                return {false, current, iter, width, "Krawczyk iteration produced empty or diverging interval"};
            }
            
            current = next;
        }
        
        return {false, current, maxIterations_, current.width(), "Maximum iterations reached"};
    }
    
private:
    Function f_;
    JacobianFunction jacobian_;
    double tolerance_;
    int maxIterations_;
};

// Verification Function: Check if f(X) is Strictly Inside Y
template<size_t N>
bool verify_inner_inclusion(const FixedIntervalVector<N>& X, 
                           const FixedIntervalVector<N>& Y)
{
    // Calculate f(X)
    auto f = [](const FixedIntervalVector<N>& x) -> FixedIntervalVector<N> {
        FixedIntervalVector<N> result;
        if constexpr (N == 2)
        {
            result[0] = x[0] * x[0] + x[1] * x[1];
            result[1] = x[0] - x[1];
        }
        return result;
    };
    
    FixedIntervalVector<N> fX = f(X);
    
    std::cout << "   f(X) = " << fX << std::endl;
    std::cout << "   Y    = " << Y << std::endl;
    
    // Check if fX is strictly contained within the dual of Y
    // For improper Y intervals, the dual is a proper interval [Y.upper(), Y.lower()]
    for (size_t i = 0; i < N; ++i)
    {
        if (!Y[i].isImproper())
        {
            std::cout << "   Error: Y[" << i << "] is not an improper interval" << std::endl;
            return false;
        }
        
        ik::KaucherInterval Y_dual = Y[i].dual(); // Convert to proper interval [Y.upper(), Y.lower()]
        
        if (!fX[i].isProper())
        {
            std::cout << "   Error: f(X)[" << i << "] is not a proper interval" << std::endl;
            return false;
        }
        
        // Check if fX(i) is strictly contained within Y_dual
        bool contains = Y_dual.contains(fX[i]);
        std::cout << "   f(X)[" << i << "] contained in Y_dual[" << i << "]: " << (contains ? "Yes" : "No") << std::endl;
        
        if (!contains)
        {
            return false;
        }
    }
    
    return true;
}

// Homotopy Continuation Method for Tolerance Embedding Problem
void testToleranceEmbeddingHomotopy()
{
    std::cout << "=== Test: Tolerance Embedding Problem (2D) - Homotopy Continuation Method ===\n";
    
    // Define Problem:
    // f1(x1, x2) = x1^2 + x2^2 = 4 + e1
    // f2(x1, x2) = x1 - x2 = 0 + e2
    // where e1 ∈ [-0.1, 0.1], e2 ∈ [-0.1, 0.1]
    
    // Target Y: Using Improper Intervals for Dual Representation
    FixedIntervalVector<2> Y;
    Y[0] = ik::KaucherInterval(3.9, 4.1).dual(); // Improper interval [4.1, 3.9]
    Y[1] = ik::KaucherInterval(-0.1, 0.1).dual(); // Improper interval [0.1, -0.1]
    
    std::cout << "Target Y (improper): " << Y << std::endl;
    
    // 定义原函数f
    auto f = [](const FixedIntervalVector<2>& x) -> FixedIntervalVector<2> {
        FixedIntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
    };
    
    // 定义原函数f的点版本
    auto point_f = [](double x1, double x2) -> std::array<double, 2> {
        return {
            x1 * x1 + x2 * x2,
            x1 - x2
        };
    };
    
    // 定义雅可比矩阵
    auto jacobian = [](const FixedIntervalVector<2>& x) -> FixedIntervalMatrix<2, 2> {
        FixedIntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * ik::KaucherInterval(2.0);
        J(0, 1) = x[1] * ik::KaucherInterval(2.0);
        J(1, 0) = ik::KaucherInterval(1.0);
        J(1, 1) = ik::KaucherInterval(-1.0);
        return J;
    };
    
    // 1. Construct auxiliary system g(x) = x - x_start = 0
    std::cout << "\n1. Constructing auxiliary system and homotopy function...\n";
    std::pair<double, double> x_start = {2.5, 2.5}; // Selected center point
    std::cout << "   Auxiliary system center point: (" << x_start.first << ", " << x_start.second << ")" << std::endl;
    
    auto g = [x_start](const FixedIntervalVector<2>& x) -> FixedIntervalVector<2> {
        FixedIntervalVector<2> result;
        result[0] = x[0] - ik::KaucherInterval(x_start.first);
        result[1] = x[1] - ik::KaucherInterval(x_start.second);
        return result;
    };
    
    // 2. 定义同伦函数H(x, t) = t·f(x) + (1-t)·g(x)
    auto createHomotopyFunction = [f, g](double t) -> FixedGeneralizedKrawczykSolver<2>::Function {
        return [t, f, g](const FixedIntervalVector<2>& x) -> FixedIntervalVector<2> {
            return t * f(x) + (1.0 - t) * g(x);
        };
    };
    
    auto createHomotopyJacobian = [jacobian, x_start](double t) -> FixedGeneralizedKrawczykSolver<2>::JacobianFunction {
        return [t, jacobian](const FixedIntervalVector<2>& x) -> FixedIntervalMatrix<2, 2> {
            FixedIntervalMatrix<2, 2> J = jacobian(x);
            // 雅可比矩阵为 t*J_f + (1-t)*I
            FixedIntervalMatrix<2, 2> I = FixedIntervalMatrix<2, 2>::identity();
            return t * J + (1.0 - t) * I;
        };
    };
    
    // 3. Gradually increase t from 0 to 1
    std::cout << "\n2. Homotopy continuation process (t from 0 to 1 in 10 steps)...\n";
    int steps = 10;
    FixedIntervalVector<2> currentBox;
    bool converged = true;
    
    for (int step = 0; step <= steps; ++step)
    {
        double t = static_cast<double>(step) / steps;
        std::cout << "\n   Step " << step << "/" << steps << " (t = " << t << "):\n";
        
        // 构造当前t值下的同伦函数和雅可比
        auto H = createHomotopyFunction(t);
        auto H_jacobian = createHomotopyJacobian(t);
        
        // Construct initial interval
        FixedIntervalVector<2> initialBox;
        if (step == 0)
        {
            // First step t=0, using auxiliary system, interval centered at x_start
            double delta = 0.1; // Initial interval radius
            initialBox[0] = ik::KaucherInterval(x_start.first - delta, x_start.first + delta);
            initialBox[1] = ik::KaucherInterval(x_start.second - delta, x_start.second + delta);
        }
        else
        {
            // Subsequent steps, using previous solution as center
            double delta = 0.1; // Interval radius
            auto c = currentBox.midpoint();
            initialBox[0] = ik::KaucherInterval(c(0) - delta, c(0) + delta);
            initialBox[1] = ik::KaucherInterval(c(1) - delta, c(1) + delta);
        }
        
        std::cout << "   Initial interval: " << initialBox << std::endl;
        
        // 运行Krawczyk迭代
        FixedGeneralizedKrawczykSolver<2> solver(H, H_jacobian, 1e-8, 50);
        auto result = solver.solve(initialBox);
        
        if (result.success)
        {
            std::cout << "   ✓ Krawczyk iteration successful!\n";
            std::cout << "   Iterations: " << result.iterations << std::endl;
            std::cout << "   Solution interval: " << result.solution << std::endl;
            std::cout << "   Interval width: " << result.finalWidth << std::endl;
            currentBox = result.solution;
        }
        else
        {
            std::cout << "   ✗ Krawczyk iteration failed!\n";
            std::cout << "   Message: " << result.message << std::endl;
            std::cout << "   Final interval: " << result.solution << std::endl;
            converged = false;
            break;
        }
    }
    
    if (converged)
    {
        std::cout << "\n3. Homotopy continuation process completed successfully!\n";
        std::cout << "   Final solution interval: " << currentBox << std::endl;
        
        // 4. Verify inner inclusion
        std::cout << "\n4. Verifying inner inclusion...\n";
        bool inner_inclusion = verify_inner_inclusion(currentBox, Y);
        
        std::cout << "\n5. Final result: " << (inner_inclusion ? "✓ Inner inclusion verified" : "✗ Inner inclusion failed") << std::endl;
        
        // 5. Additional verification: Check four corners
        std::cout << "\n6. Corner verification...\n";
        std::array<std::pair<double, double>, 4> corners;
        corners[0] = std::make_pair(currentBox[0].lower(), currentBox[1].lower());
        corners[1] = std::make_pair(currentBox[0].lower(), currentBox[1].upper());
        corners[2] = std::make_pair(currentBox[0].upper(), currentBox[1].lower());
        corners[3] = std::make_pair(currentBox[0].upper(), currentBox[1].upper());
        
        for (size_t i = 0; i < corners.size(); ++i)
        {
            auto [x1, x2] = corners[i];
            auto fx = point_f(x1, x2);
            
            // Check if inside target Y (considering tolerance)
            bool in_y1 = (fx[0] >= 3.9 && fx[0] <= 4.1);
            bool in_y2 = (fx[1] >= -0.1 && fx[1] <= 0.1);
            
            std::cout << "   Corner " << i+1 << " (" << x1 << ", " << x2 << "): "
                      << "f = (" << fx[0] << ", " << fx[1] << ") - "
                      << (in_y1 && in_y2 ? "✓ Valid" : "✗ Invalid") << std::endl;
        }
    }
    else
    {
        std::cout << "\n3. Homotopy continuation process failed!\n";
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Tolerance Embedding Homotopy Continuation Implementation ===\n\n";
    
    // Demonstrate Basic Properties of Kaucher Interval Arithmetic
    std::cout << "=== Kaucher Interval Arithmetic Properties ===\n";
    ik::KaucherInterval proper(1.0, 2.0);
    ik::KaucherInterval improper(2.0, 1.0);
    
    std::cout << "Proper interval a = " << proper << std::endl;
    std::cout << "Improper interval b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "Inverse of a = " << (-proper) << std::endl;
    std::cout << "Dual of b = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // Run homotopy continuation test
    testToleranceEmbeddingHomotopy();
    
    std::cout << "=== Testing Complete ===\n";
    
    return 0;
}