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

// Newton method pre-iteration function for finding high-precision center point
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
        // 计算函数值
        auto fx = f(x1, x2);
        double fx1 = fx[0];
        double fx2 = fx[1];
        
        // 检查收敛
        if (std::abs(fx1) < tol && std::abs(fx2) < tol)
        {
            break;
        }
        
        // 计算雅可比矩阵
        auto J = df(x1, x2);
        double J11 = J[0][0];
        double J12 = J[0][1];
        double J21 = J[1][0];
        double J22 = J[1][1];
        
        // 计算行列式
        double det = J11 * J22 - J12 * J21;
        if (std::abs(det) < 1e-12)
        {
            break;
        }
        
        // 计算逆雅可比矩阵
        double invJ11 = J22 / det;
        double invJ12 = -J12 / det;
        double invJ21 = -J21 / det;
        double invJ22 = J11 / det;
        
        // 更新x
        double delta1 = -(invJ11 * fx1 + invJ12 * fx2);
        double delta2 = -(invJ21 * fx1 + invJ22 * fx2);
        
        x1 += delta1;
        x2 += delta2;
        
        // 检查步长收敛
        if (std::abs(delta1) < tol && std::abs(delta2) < tol)
        {
            break;
        }
    }
    
    return std::make_pair(x1, x2);
}

// Verification function: Check if f(X) is strictly inside Y
template<size_t N>
bool verify_inner_inclusion(const IntervalVector<N>& X, 
                           const IntervalVector<N>& Y)
{
    // 计算f(X)
    auto f = [](const IntervalVector<N>& x) -> IntervalVector<N> {
        IntervalVector<N> result;
        if constexpr (N == 1)
        {
            result[0] = x[0] * x[0];
        }
        else if constexpr (N == 2)
        {
            result[0] = x[0] * x[0] + x[1] * x[1];
            result[1] = x[0] - x[1];
        }
        return result;
    };
    
    IntervalVector<N> fX = f(X);
    
    // 检查fX是否严格包含在Y的对偶区间内
    // 对于Y是improper区间的情况，其对偶区间是proper区间[Y.upper(), Y.lower()]
    for (size_t i = 0; i < N; ++i)
    {
        if (!Y[i].isImproper())
        {
            std::cout << "   Error: Y[" << i << "] is not an improper interval" << std::endl;
            return false;
        }
        
        KaucherInterval Y_dual = Y[i].dual(); // 转换为proper区间 [Y.upper(), Y.lower()]
        
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

// Test: Tolerance Embedding Problem (2D) - Direct Verification
void testToleranceEmbeddingDirect()
{
    std::cout << "=== Test: Tolerance Embedding Problem (2D) - Direct Verification ===\n";
    
    // Define problem:
    // f1(x1, x2) = x1^2 + x2^2 = 4 + e1
    // f2(x1, x2) = x1 - x2 = 0 + e2
    // where e1 ∈ [-0.1, 0.1], e2 ∈ [-0.1, 0.1]
    
    // Target Y: Using improper intervals for dual representation
    IntervalVector<2> Y;
    Y[0] = KaucherInterval(3.9, 4.1).dual(); // 非正常区间 [4.1, 3.9]
    Y[1] = KaucherInterval(-0.1, 0.1).dual(); // 非正常区间 [0.1, -0.1]
    
    std::cout << "Target Y (improper): " << "[" << Y[0] << ", " << Y[1] << "]" << std::endl;
    
    // 1. Use Newton method pre-iteration to find high-precision center point x*
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
    
    // 运行牛顿法
    auto [x1_star, x2_star] = newton_method_2d(x1_init, x2_init, point_f, point_df);
    
    std::cout << "   Initial guess: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   Newton result: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // Verify center point
    auto f_center = point_f(x1_star, x2_star);
    std::cout << "   f(x*) = (" << f_center[0] << ", " << f_center[1] << ")" << std::endl;
    
    // 2. Construct small initial interval
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.05; // 小区间半径
    IntervalVector<2> X;
    X[0] = KaucherInterval(x1_star - delta, x1_star + delta);
    X[1] = KaucherInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   Initial interval X: " << "[" << X[0] << ", " << X[1] << "]" << std::endl;
    std::cout << "   Interval width: " << X[0].width() << " × " << X[1].width() << std::endl;
    
    // 3. Define function f
    auto f = [](const IntervalVector<2>& x) -> IntervalVector<2> {
        IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
    };
    
    // 4. Directly compute f(X) and verify inner inclusion
    std::cout << "\n3. Directly computing f(X) and verifying inner inclusion...\n";
    IntervalVector<2> fX = f(X);
    
    std::cout << "   f(X) = " << "[" << fX[0] << ", " << fX[1] << "]" << std::endl;
    
    // 5. Verify inner inclusion
    std::cout << "\n4. Verifying inner inclusion...\n";
    bool inner_inclusion = verify_inner_inclusion(X, Y);
    
    std::cout << "\n5. Final result: " << (inner_inclusion ? "✓ Inner inclusion verified" : "✗ Inner inclusion failed") << std::endl;
    
    // 6. Additional verification: Check four corners
    std::cout << "\n6. Corner verification...\n";
    std::array<std::pair<double, double>, 4> corners;
    corners[0] = std::make_pair(X[0].lower(), X[1].lower());
    corners[1] = std::make_pair(X[0].lower(), X[1].upper());
    corners[2] = std::make_pair(X[0].upper(), X[1].lower());
    corners[3] = std::make_pair(X[0].upper(), X[1].upper());
    
    for (size_t i = 0; i < corners.size(); ++i)
    {
        auto [x1, x2] = corners[i];
        auto fx = point_f(x1, x2);
        
        // 检查是否在目标Y内（考虑容差）
        bool in_y1 = (fx[0] >= 3.9 && fx[0] <= 4.1);
        bool in_y2 = (fx[1] >= -0.1 && fx[1] <= 0.1);
        
        std::cout << "   Corner " << i+1 << " (" << x1 << ", " << x2 << "): "
                  << "f = (" << fx[0] << ", " << fx[1] << ") - "
                  << (in_y1 && in_y2 ? "✓ Valid" : "✗ Invalid") << std::endl;
    }
    
    std::cout << std::endl;
}

// Test: 1D Tolerance Embedding Problem - Direct Verification
void testToleranceEmbedding1DDirect()
{
    std::cout << "=== Test: 1D Tolerance Embedding Problem - Direct Verification ===\n";
    
    // Define problem: f(x) = x^2 = 2 + e, where e ∈ [-0.1, 0.1]
    
    // Target Y: Using improper intervals for dual representation
    IntervalVector<1> Y;
    Y[0] = KaucherInterval(1.9, 2.1).dual(); // 非正常区间 [2.1, 1.9]
    
    std::cout << "Target Y (improper): " << "[" << Y[0] << "]" << std::endl;
    
    // 1. Use Newton method pre-iteration
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
    
    // 2. Construct small initial interval
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.05;
    IntervalVector<1> X;
    X[0] = KaucherInterval(x_star - delta, x_star + delta);
    
    std::cout << "   Initial interval X: " << "[" << X[0] << "]" << std::endl;
    
    // 3. Directly compute f(X) and verify inner inclusion
    std::cout << "\n3. Directly computing f(X) and verifying inner inclusion...\n";
    
    // 定义函数f
    auto f = [](const IntervalVector<1>& x) -> IntervalVector<1> {
        IntervalVector<1> result;
        result[0] = x[0] * x[0];
        return result;
    };
    
    IntervalVector<1> fX = f(X);
    std::cout << "   f(X) = " << "[" << fX[0] << "]" << std::endl;
    
    // 4. Verify inner inclusion
    std::cout << "\n4. Verifying inner inclusion...\n";
    bool inner_inclusion = verify_inner_inclusion(X, Y);
    
    std::cout << "\n5. Final result: " << (inner_inclusion ? "✓ Inner inclusion verified" : "✗ Inner inclusion failed") << std::endl;
    
    // 6. Additional verification: Check boundary points
    std::cout << "\n6. Boundary point verification...\n";
    std::array<double, 2> bounds = {X[0].lower(), X[0].upper()};
    
    for (size_t i = 0; i < bounds.size(); ++i)
    {
        double x = bounds[i];
        double fx = f_1d(x) + 2.0; // 因为f_1d(x) = x^2 - 2.0，所以x^2 = f_1d(x) + 2.0
        
        // 检查是否在目标Y内（考虑容差）
        bool in_y = (fx >= 1.9 && fx <= 2.1);
        
        std::cout << "   Boundary point " << i+1 << " (" << x << "): "
                  << "f = " << fx << " - "
                  << (in_y ? "✓ Valid" : "✗ Invalid") << std::endl;
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Tolerance Embedding Final Implementation ===\n\n";
    
    // Demonstrate basic properties of Kaucher interval arithmetic
    std::cout << "=== Kaucher Interval Arithmetic Properties ===\n";
    KaucherInterval proper(1.0, 2.0);
    KaucherInterval improper(2.0, 1.0);
    
    std::cout << "Proper interval a = " << proper << std::endl;
    std::cout << "Improper interval b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "Inverse of a = " << (-proper) << std::endl;
    std::cout << "Dual of b = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // Run tests
    testToleranceEmbedding1DDirect(); // 1D test
    testToleranceEmbeddingDirect();   // 2D test
    
    std::cout << "=== Testing Complete ===\n";
    
    return 0;
}