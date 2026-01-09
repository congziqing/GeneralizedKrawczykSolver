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

// Generalized Krawczyk Solver (Fixed Dimension Version)
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
            
            // 计算中点c
            Eigen::VectorXd c = current.midpoint();
            
            // 创建中点的区间向量
            FixedIntervalVector<N> c_vec;
            for (size_t i = 0; i < N; ++i)
            {
                c_vec[i] = ik::KaucherInterval(c(i), c(i));
            }
            
            // 计算f(c)
            FixedIntervalVector<N> f_c = f_(c_vec);
            
            // 计算雅可比矩阵J(x)
            FixedIntervalMatrix<N, N> J = jacobian_(current);
            
            // 计算雅可比矩阵的中点
            Eigen::MatrixXd J_mid = J.midpoint();
            
            // 检查雅可比矩阵是否奇异
            double det = J_mid.determinant();
            if (std::abs(det) < 1e-12)
            {
                return {false, current, iter, width, "Singular Jacobian matrix"};
            }
            
            // 计算雅可比矩阵的逆
            Eigen::MatrixXd Y = J_mid.inverse();
            
            // 计算term1: c - Y * f(c)
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
            
            // 计算I - Y * J(x)
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(static_cast<int>(N), static_cast<int>(N));
            Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
            
            // 创建区间矩阵I_minus_YJ_interval
            FixedIntervalMatrix<N, N> I_minus_YJ_interval;
            for (size_t i = 0; i < N; ++i)
            {
                for (size_t j = 0; j < N; ++j)
                {
                    double val = I_minus_YJ(static_cast<int>(i), static_cast<int>(j));
                    I_minus_YJ_interval(i, j) = ik::KaucherInterval(val, val);
                }
            }
            
            // 计算box - c
            FixedIntervalVector<N> box_minus_c;
            for (size_t i = 0; i < N; ++i)
            {
                double delta = (current[i].upper() - current[i].lower()) / 2.0;
                box_minus_c[i] = ik::KaucherInterval(-delta, delta);
            }
            
            // 计算term3: (I - YJ) * (box - c)
            FixedIntervalVector<N> term3 = I_minus_YJ_interval * box_minus_c;
            
            // 计算Krawczyk区间K = term1 + term3
            FixedIntervalVector<N> K = term1 + term3;
            
            // 计算交集K ∩ current
            FixedIntervalVector<N> next;
            for (size_t i = 0; i < N; ++i)
            {
                next[i] = current[i].meet(K[i]);
            }
            
            // 检查是否收敛或发散
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
                           const FixedIntervalVector<N>& Y, 
                           const typename FixedGeneralizedKrawczykSolver<N>::Function& f)
{
    // 计算f(X)
    FixedIntervalVector<N> fX = f(X);
    
    // 检查fX是否严格包含在Y的内部
    for (size_t i = 0; i < N; ++i)
    {
        if (!fX[i].isProper() || !Y[i].isProper())
        {
            // 如果是improper区间，需要特殊处理
            // 对于内包围问题，Y应该是improper区间，fX应该是proper区间
            if (fX[i].isProper() && !Y[i].isProper())
            {
                // 检查fX是否包含在Y的对偶区间内
                ik::KaucherInterval Y_dual = Y[i].dual();
                if (!(fX[i].lower() > Y_dual.lower() && fX[i].upper() < Y_dual.upper()))
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        else
        {
            // proper区间的情况：fX应该严格包含在Y内部
            if (!(fX[i].lower() > Y[i].lower() && fX[i].upper() < Y[i].upper()))
            {
                return false;
            }
        }
    }
    
    return true;
}

// 测试：Tolerance Embedding Problem (2D)
void testToleranceEmbedding2D()
{
    std::cout << "=== Test: Tolerance Embedding Problem (2D) ===\n";
    
    // Define Problem:
    // f1(x1, x2) = x1^2 + x2^2 = 4 + e1
    // f2(x1, x2) = x1 - x2 = 0 + e2
    // where e1 ∈ [-0.1, 0.1], e2 ∈ [-0.1, 0.1]
    
    // Target Y: Using Improper Intervals for Dual Representation
    FixedIntervalVector<2> Y;
    Y[0] = ik::KaucherInterval(3.9, 4.1).dual(); // Improper interval [4.1, 3.9]
    Y[1] = ik::KaucherInterval(-0.1, 0.1).dual(); // Improper interval [0.1, -0.1]
    
    std::cout << "Target Y (improper): " << Y << std::endl;
    
    // 定义函数f：f(X) = Y
    auto f = [](const FixedIntervalVector<2>& x) -> FixedIntervalVector<2> {
        FixedIntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
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
    
    // 1. Use Newton Method Pre-iteration to Find High-Precision Center Point x*
    std::cout << "\n1. Newton method pre-iteration to find center point...\n";
    
    // 定义点函数f（用于牛顿法）
    auto point_f = [](double x1, double x2) -> std::array<double, 2> {
        return {
            x1 * x1 + x2 * x2 - 4.0,
            x1 - x2
        };
    };
    
    // 定义雅可比矩阵的点函数
    auto point_df = [](double x1, double x2) -> std::array<std::array<double, 2>, 2> {
        std::array<std::array<double, 2>, 2> result;
        result[0][0] = 2.0 * x1;
        result[0][1] = 2.0 * x2;
        result[1][0] = 1.0;
        result[1][1] = -1.0;
        return result;
    };
    
    // 初始猜测点
    double x1_init = 1.5;
    double x2_init = 1.5;
    
    // Run Newton Method
    auto [x1_star, x2_star] = newton_method_2d(x1_init, x2_init, point_f, point_df);
    
    std::cout << "   Initial guess: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   Newton result: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // Verify Center Point
    auto f_center = point_f(x1_star, x2_star);
    std::cout << "   f(x*) = (" << f_center[0] << ", " << f_center[1] << ")" << std::endl;
    
    // 2. Construct Small Initial Interval
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.1; // Small interval radius
    FixedIntervalVector<2> X0;
    X0[0] = ik::KaucherInterval(x1_star - delta, x1_star + delta);
    X0[1] = ik::KaucherInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   Initial interval X0: " << X0 << std::endl;
    std::cout << "   Interval width: " << X0[0].magnitude() << " × " << X0[1].magnitude() << std::endl;
    
    // 3. Run Generalized Krawczyk Iteration
    std::cout << "\n3. Running generalized Krawczyk iteration...\n";
    FixedGeneralizedKrawczykSolver<2> solver(f, jacobian, 1e-8, 50);
    
    auto result = solver.solve(X0);
    
    if (result.success)
    {
        std::cout << "   ✓ Krawczyk iteration successful!\n";
        std::cout << "   Iterations: " << result.iterations << std::endl;
        std::cout << "   Final interval X_final: " << result.solution << std::endl;
        std::cout << "   Final width: " << result.finalWidth << std::endl;
        
        // 4. Inner Inclusion Verification
        std::cout << "\n4. Verifying inner inclusion...\n";
        
        // Calculate f(X_final)
        FixedIntervalVector<2> fX = f(result.solution);
        std::cout << "   f(X_final) = " << fX << std::endl;
        std::cout << "   Y = " << Y << std::endl;
        
        // Verify Inclusion: Check if f(X_final) is Strictly Inside Y
        bool inner_inclusion = verify_inner_inclusion(result.solution, Y, f);
        std::cout << "   Inner inclusion result: " << (inner_inclusion ? "✓ Valid" : "✗ Invalid") << std::endl;
        
        // Additional Verification: Check Four Corners
        std::cout << "\n5. Corner verification...\n";
        std::array<std::pair<double, double>, 4> corners;
        corners[0] = std::make_pair(result.solution[0].lower(), result.solution[1].lower());
        corners[1] = std::make_pair(result.solution[0].lower(), result.solution[1].upper());
        corners[2] = std::make_pair(result.solution[0].upper(), result.solution[1].lower());
        corners[3] = std::make_pair(result.solution[0].upper(), result.solution[1].upper());
        
        for (size_t i = 0; i < corners.size(); ++i)
        {
            auto [x1, x2] = corners[i];
            auto fx = point_f(x1, x2);
            
            // Check if Inside Target Y (Considering Tolerance)
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
    
    std::cout << std::endl;
}

// Test: 1D Tolerance Embedding Problem (Simple Example)
void testToleranceEmbedding1D()
{
    std::cout << "=== Test: 1D Tolerance Embedding Problem ===\n";
    
    // Define Problem: f(x) = x^2 = 2 + e, where e ∈ [-0.1, 0.1]
    
    // Target Y: Using Improper Intervals for Dual Representation
    FixedIntervalVector<1> Y;
    Y[0] = ik::KaucherInterval(1.9, 2.1).dual(); // Improper interval [2.1, 1.9]
    
    std::cout << "Target Y (improper): " << Y << std::endl;
    
    // 定义函数f
    auto f = [](const FixedIntervalVector<1>& x) -> FixedIntervalVector<1> {
        FixedIntervalVector<1> result;
        result[0] = x[0] * x[0];
        return result;
    };
    
    // 定义雅可比矩阵
    auto jacobian = [](const FixedIntervalVector<1>& x) -> FixedIntervalMatrix<1, 1> {
        FixedIntervalMatrix<1, 1> J;
        J(0, 0) = x[0] * ik::KaucherInterval(2.0);
        return J;
    };
    
    // 1. Use Newton Method Pre-iteration
    std::cout << "\n1. Newton method pre-iteration to find center point...\n";
    
    auto f_1d = [](double x) -> double {
        return x * x - 2.0;
    };
    
    auto df_1d = [](double x) -> double {
        return 2.0 * x;
    };
    
    double x0 = 1.4;
    double x_star = newton_method_1d(x0, f_1d, df_1d);
    
    std::cout << "   Initial guess: " << x0 << std::endl;
    std::cout << "   Newton result: x* = " << x_star << std::endl;
    std::cout << "   f(x*) = " << f_1d(x_star) << std::endl;
    
    // 2. Construct Small Initial Interval
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.1;
    FixedIntervalVector<1> X0;
    X0[0] = ik::KaucherInterval(x_star - delta, x_star + delta);
    
    std::cout << "   Initial interval X0: " << X0 << std::endl;
    
    // 3. Run Generalized Krawczyk Iteration
    std::cout << "\n3. Running generalized Krawczyk iteration...\n";
    FixedGeneralizedKrawczykSolver<1> solver(f, jacobian, 1e-8, 50);
    
    auto result = solver.solve(X0);
    
    if (result.success)
    {
        std::cout << "   ✓ Krawczyk iteration successful!\n";
        std::cout << "   Iterations: " << result.iterations << std::endl;
        std::cout << "   Final interval X_final: " << result.solution << std::endl;
        
        // 4. Inner Inclusion Verification
        std::cout << "\n4. Verifying inner inclusion...\n";
        FixedIntervalVector<1> fX = f(result.solution);
        std::cout << "   f(X_final) = " << fX << std::endl;
        std::cout << "   Y = " << Y << std::endl;
        
        bool inner_inclusion = verify_inner_inclusion(result.solution, Y, f);
        std::cout << "   Inner inclusion result: " << (inner_inclusion ? "✓ Valid" : "✗ Invalid") << std::endl;
    }
    else
    {
        std::cout << "   ✗ Krawczyk iteration failed!\n";
        std::cout << "   Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Tolerance Embedding Fixed Dimension Implementation ===\n\n";
    
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
    
    // Run Tests
    testToleranceEmbedding1D(); // 1D test
    testToleranceEmbedding2D();   // 2D test
    
    std::cout << "=== Testing Complete ===\n";
    
    return 0;
}