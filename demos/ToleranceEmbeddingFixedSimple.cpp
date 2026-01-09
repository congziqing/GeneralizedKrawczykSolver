#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <Eigen/Dense>
#include "PSGMDirectedInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"

// 牛顿法预迭代函数，用于找到高精度的中心点
auto newton_method_2d(double x1_init, double x2_init)
{
    // 定义函数f(x1, x2) = (x1^2 + x2^2 - 4, x1 - x2)
    auto f = [](double x1, double x2) -> std::array<double, 2> {
        return {
            x1 * x1 + x2 * x2 - 4.0,
            x1 - x2
        };
    };
    
    // 定义雅可比矩阵的点函数
    auto df = [](double x1, double x2) -> std::array<std::array<double, 2>, 2> {
        std::array<std::array<double, 2>, 2> result;
        result[0][0] = 2.0 * x1;
        result[0][1] = 2.0 * x2;
        result[1][0] = 1.0;
        result[1][1] = -1.0;
        return result;
    };
    
    // 牛顿迭代
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
    std::cout << "=== Tolerance Embedding 简单实现 ===\n\n";
    
    // 1. 使用牛顿法预迭代找到高精度中心点x*
    std::cout << "1. 牛顿法预迭代寻找中心点...\n";
    double x1_init = 1.5;
    double x2_init = 1.5;
    auto [x1_star, x2_star] = newton_method_2d(x1_init, x2_init);
    
    std::cout << "   初始猜测: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   牛顿法结果: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // 2. 构造小初始区间
    std::cout << "\n2. 构造小初始区间...\n";
    double delta = 0.01; // 使用更小的初始区间半径
    FixedIntervalVector<2> X;
    X[0] = PSGMDirectedInterval(x1_star - delta, x1_star + delta);
    X[1] = PSGMDirectedInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   初始区间X: " << X << std::endl;
    
    // 3. 定义函数f
    auto f = [](const FixedIntervalVector<2>& x) -> FixedIntervalVector<2> {
        FixedIntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
    };
    
    // 4. 计算f(X)
    std::cout << "\n3. 计算f(X)...\n";
    FixedIntervalVector<2> fX = f(X);
    
    std::cout << "   f(X) = " << fX << std::endl;
    
    // 5. 定义目标Y（improper区间）
    std::cout << "\n4. 定义目标Y（improper区间）...\n";
    FixedIntervalVector<2> Y;
    Y[0] = PSGMDirectedInterval(3.9, 4.1).dual(); // [4.1, 3.9]
    Y[1] = PSGMDirectedInterval(-0.1, 0.1).dual(); // [0.1, -0.1]
    
    std::cout << "   Y = " << Y << std::endl;
    
    // 6. 验证内包围
    std::cout << "\n5. 验证内包围...\n";
    
    // 检查fX是否严格包含在Y的对偶区间内
    bool inner_inclusion = true;
    for (size_t i = 0; i < 2; ++i)
    {
        PSGMDirectedInterval Y_dual = Y[i].dual(); // 转换为proper区间
        std::cout << "   Y_dual[" << i << "] = " << Y_dual << std::endl;
        std::cout << "   fX[" << i << "] = " << fX[i] << std::endl;
        
        bool contains = Y_dual.contains(fX[i]);
        std::cout << "   fX[" << i << "] 是否包含在 Y_dual[" << i << "] 中: " << (contains ? "是" : "否") << std::endl;
        
        if (!contains)
        {
            inner_inclusion = false;
        }
    }
    
    // 7. 最终结果
    std::cout << "\n6. 最终结果: " << (inner_inclusion ? "✓ 内包围验证通过" : "✗ 内包围验证失败") << std::endl;
    
    return 0;
}