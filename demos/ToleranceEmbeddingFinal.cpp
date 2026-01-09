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

// 牛顿法预迭代函数，用于找到高精度的中心点
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

// 2D牛顿法预迭代函数
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

// 验证函数：检查f(X)是否严格位于Y的内部
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
            std::cout << "   错误：Y[" << i << "] 不是improper区间" << std::endl;
            return false;
        }
        
        KaucherInterval Y_dual = Y[i].dual(); // 转换为proper区间 [Y.upper(), Y.lower()]
        
        if (!fX[i].isProper())
        {
            std::cout << "   错误：f(X)[" << i << "] 不是proper区间" << std::endl;
            return false;
        }
        
        // 检查fX(i)是否严格包含在Y_dual内部
        bool contains = Y_dual.contains(fX[i]);
        std::cout << "   f(X)[" << i << "] 是否包含在 Y_dual[" << i << "] 中: " << (contains ? "是" : "否") << std::endl;
        
        if (!contains)
        {
            return false;
        }
    }
    
    return true;
}

// 测试：Tolerance Embedding Problem (2D) - 直接验证
void testToleranceEmbeddingDirect()
{
    std::cout << "=== 测试: Tolerance Embedding Problem (2D) - 直接验证 ===\n";
    
    // 定义问题：
    // f1(x1, x2) = x1^2 + x2^2 = 4 + e1
    // f2(x1, x2) = x1 - x2 = 0 + e2
    // 其中 e1 ∈ [-0.1, 0.1], e2 ∈ [-0.1, 0.1]
    
    // 目标Y：使用improper区间表示对偶形式
    IntervalVector<2> Y;
    Y[0] = KaucherInterval(3.9, 4.1).dual(); // 非正常区间 [4.1, 3.9]
    Y[1] = KaucherInterval(-0.1, 0.1).dual(); // 非正常区间 [0.1, -0.1]
    
    std::cout << "目标Y（improper）: " << "[" << Y[0] << ", " << Y[1] << "]" << std::endl;
    
    // 1. 使用牛顿法预迭代找到高精度中心点x*
    std::cout << "\n1. 牛顿法预迭代寻找中心点...\n";
    
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
    
    std::cout << "   初始猜测: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   牛顿法结果: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // 验证中心点
    auto f_center = point_f(x1_star, x2_star);
    std::cout << "   f(x*) = (" << f_center[0] << ", " << f_center[1] << ")" << std::endl;
    
    // 2. 构造小初始区间
    std::cout << "\n2. 构造小初始区间...\n";
    double delta = 0.05; // 小区间半径
    IntervalVector<2> X;
    X[0] = KaucherInterval(x1_star - delta, x1_star + delta);
    X[1] = KaucherInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   初始区间X: " << "[" << X[0] << ", " << X[1] << "]" << std::endl;
    std::cout << "   区间宽度: " << X[0].width() << " × " << X[1].width() << std::endl;
    
    // 3. 定义函数f
    auto f = [](const IntervalVector<2>& x) -> IntervalVector<2> {
        IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1];
        result[1] = x[0] - x[1];
        return result;
    };
    
    // 4. 直接计算f(X)并验证内包围
    std::cout << "\n3. 直接计算f(X)并验证内包围...\n";
    IntervalVector<2> fX = f(X);
    
    std::cout << "   f(X) = " << "[" << fX[0] << ", " << fX[1] << "]" << std::endl;
    
    // 5. 验证内包围
    std::cout << "\n4. 内包围验证...\n";
    bool inner_inclusion = verify_inner_inclusion(X, Y);
    
    std::cout << "\n5. 最终结果: " << (inner_inclusion ? "✓ 内包围验证通过" : "✗ 内包围验证失败") << std::endl;
    
    // 6. 额外验证：检查四个角点
    std::cout << "\n6. 角点验证...\n";
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
        
        std::cout << "   角点 " << i+1 << " (" << x1 << ", " << x2 << "): "
                  << "f = (" << fx[0] << ", " << fx[1] << ") - "
                  << (in_y1 && in_y2 ? "✓ 有效" : "✗ 无效") << std::endl;
    }
    
    std::cout << std::endl;
}

// 测试：1D Tolerance Embedding Problem - 直接验证
void testToleranceEmbedding1DDirect()
{
    std::cout << "=== 测试: 1D Tolerance Embedding Problem - 直接验证 ===\n";
    
    // 定义问题：f(x) = x^2 = 2 + e，其中 e ∈ [-0.1, 0.1]
    
    // 目标Y：使用improper区间表示对偶形式
    IntervalVector<1> Y;
    Y[0] = KaucherInterval(1.9, 2.1).dual(); // 非正常区间 [2.1, 1.9]
    
    std::cout << "目标Y（improper）: " << "[" << Y[0] << "]" << std::endl;
    
    // 1. 使用牛顿法预迭代
    std::cout << "\n1. 牛顿法预迭代寻找中心点...\n";
    
    auto f_1d = [](double x) -> double {
        return x * x - 2.0;
    };
    
    auto df_1d = [](double x) -> double {
        return 2.0 * x;
    };
    
    double x0 = 1.4;
    double x_star = newton_method_1d(x0, f_1d, df_1d);
    
    std::cout << "   初始猜测: " << x0 << std::endl;
    std::cout << "   牛顿法结果: x* = " << x_star << std::endl;
    std::cout << "   f(x*) = " << f_1d(x_star) << std::endl;
    
    // 2. 构造小初始区间
    std::cout << "\n2. 构造小初始区间...\n";
    double delta = 0.05;
    IntervalVector<1> X;
    X[0] = KaucherInterval(x_star - delta, x_star + delta);
    
    std::cout << "   初始区间X: " << "[" << X[0] << "]" << std::endl;
    
    // 3. 直接计算f(X)并验证内包围
    std::cout << "\n3. 直接计算f(X)并验证内包围...\n";
    
    // 定义函数f
    auto f = [](const IntervalVector<1>& x) -> IntervalVector<1> {
        IntervalVector<1> result;
        result[0] = x[0] * x[0];
        return result;
    };
    
    IntervalVector<1> fX = f(X);
    std::cout << "   f(X) = " << "[" << fX[0] << "]" << std::endl;
    
    // 4. 验证内包围
    std::cout << "\n4. 内包围验证...\n";
    bool inner_inclusion = verify_inner_inclusion(X, Y);
    
    std::cout << "\n5. 最终结果: " << (inner_inclusion ? "✓ 内包围验证通过" : "✗ 内包围验证失败") << std::endl;
    
    // 6. 额外验证：检查边界点
    std::cout << "\n6. 边界点验证...\n";
    std::array<double, 2> bounds = {X[0].lower(), X[0].upper()};
    
    for (size_t i = 0; i < bounds.size(); ++i)
    {
        double x = bounds[i];
        double fx = f_1d(x) + 2.0; // 因为f_1d(x) = x^2 - 2.0，所以x^2 = f_1d(x) + 2.0
        
        // 检查是否在目标Y内（考虑容差）
        bool in_y = (fx >= 1.9 && fx <= 2.1);
        
        std::cout << "   边界点 " << i+1 << " (" << x << "): "
                  << "f = " << fx << " - "
                  << (in_y ? "✓ 有效" : "✗ 无效") << std::endl;
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Tolerance Embedding 最终实现 ===\n\n";
    
    // 展示Kaucher区间算术的基本特性
    std::cout << "=== Kaucher区间算术特性 ===\n";
    KaucherInterval proper(1.0, 2.0);
    KaucherInterval improper(2.0, 1.0);
    
    std::cout << "正常区间 a = " << proper << std::endl;
    std::cout << "非正常区间 b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "a的逆元 = " << (-proper) << std::endl;
    std::cout << "b的对偶 = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // 运行测试
    testToleranceEmbedding1DDirect(); // 1D测试
    testToleranceEmbeddingDirect();   // 2D测试
    
    std::cout << "=== 测试完成 ===\n";
    
    return 0;
}