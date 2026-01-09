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

// 广义Krawczyk求解器（固定维数版本）
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
                c_vec[i] = PSGMDirectedInterval(c(i), c(i));
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
                term1[i] = PSGMDirectedInterval(val, val);
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
                    I_minus_YJ_interval(i, j) = PSGMDirectedInterval(val, val);
                }
            }
            
            // 计算box - c
            FixedIntervalVector<N> box_minus_c;
            for (size_t i = 0; i < N; ++i)
            {
                double delta = (current[i].upper() - current[i].lower()) / 2.0;
                box_minus_c[i] = PSGMDirectedInterval(-delta, delta);
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

// 验证函数：检查f(X)是否严格位于Y的内部
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
                PSGMDirectedInterval Y_dual = Y[i].dual();
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
    std::cout << "=== 测试: Tolerance Embedding Problem (2D) ===\n";
    
    // 定义问题：
    // f1(x1, x2) = x1^2 + x2^2 = 4 + e1
    // f2(x1, x2) = x1 - x2 = 0 + e2
    // 其中 e1 ∈ [-0.1, 0.1], e2 ∈ [-0.1, 0.1]
    
    // 目标Y：使用improper区间表示对偶形式
    FixedIntervalVector<2> Y;
    Y[0] = PSGMDirectedInterval(3.9, 4.1).dual(); // 非正常区间 [4.1, 3.9]
    Y[1] = PSGMDirectedInterval(-0.1, 0.1).dual(); // 非正常区间 [0.1, -0.1]
    
    std::cout << "目标Y（improper）: " << Y << std::endl;
    
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
        J(0, 0) = x[0] * PSGMDirectedInterval(2.0);
        J(0, 1) = x[1] * PSGMDirectedInterval(2.0);
        J(1, 0) = PSGMDirectedInterval(1.0);
        J(1, 1) = PSGMDirectedInterval(-1.0);
        return J;
    };
    
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
    double delta = 0.1; // 小区间半径
    FixedIntervalVector<2> X0;
    X0[0] = PSGMDirectedInterval(x1_star - delta, x1_star + delta);
    X0[1] = PSGMDirectedInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   初始区间X0: " << X0 << std::endl;
    std::cout << "   区间宽度: " << X0[0].magnitude() << " × " << X0[1].magnitude() << std::endl;
    
    // 3. 运行广义Krawczyk迭代
    std::cout << "\n3. 运行广义Krawczyk迭代...\n";
    FixedGeneralizedKrawczykSolver<2> solver(f, jacobian, 1e-8, 50);
    
    auto result = solver.solve(X0);
    
    if (result.success)
    {
        std::cout << "   ✓ Krawczyk迭代成功!\n";
        std::cout << "   迭代次数: " << result.iterations << std::endl;
        std::cout << "   最终区间X_final: " << result.solution << std::endl;
        std::cout << "   最终宽度: " << result.finalWidth << std::endl;
        
        // 4. 内包围验证
        std::cout << "\n4. 内包围验证...\n";
        
        // 计算f(X_final)
        FixedIntervalVector<2> fX = f(result.solution);
        std::cout << "   f(X_final) = " << fX << std::endl;
        std::cout << "   Y = " << Y << std::endl;
        
        // 验证包含关系：检查f(X_final)是否严格位于Y的内部
        bool inner_inclusion = verify_inner_inclusion(result.solution, Y, f);
        std::cout << "   内包围验证结果: " << (inner_inclusion ? "✓ 有效" : "✗ 无效") << std::endl;
        
        // 额外验证：检查四个角点
        std::cout << "\n5. 角点验证...\n";
        std::array<std::pair<double, double>, 4> corners;
        corners[0] = std::make_pair(result.solution[0].lower(), result.solution[1].lower());
        corners[1] = std::make_pair(result.solution[0].lower(), result.solution[1].upper());
        corners[2] = std::make_pair(result.solution[0].upper(), result.solution[1].lower());
        corners[3] = std::make_pair(result.solution[0].upper(), result.solution[1].upper());
        
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
        
    }
    else
    {
        std::cout << "   ✗ Krawczyk迭代失败!\n";
        std::cout << "   消息: " << result.message << std::endl;
        std::cout << "   最终区间: " << result.solution << std::endl;
        std::cout << "   最终宽度: " << result.finalWidth << std::endl;
    }
    
    std::cout << std::endl;
}

// 测试：1D Tolerance Embedding Problem（简单示例）
void testToleranceEmbedding1D()
{
    std::cout << "=== 测试: 1D Tolerance Embedding Problem ===\n";
    
    // 定义问题：f(x) = x^2 = 2 + e，其中 e ∈ [-0.1, 0.1]
    
    // 目标Y：使用improper区间表示对偶形式
    FixedIntervalVector<1> Y;
    Y[0] = PSGMDirectedInterval(1.9, 2.1).dual(); // 非正常区间 [2.1, 1.9]
    
    std::cout << "目标Y（improper）: " << Y << std::endl;
    
    // 定义函数f
    auto f = [](const FixedIntervalVector<1>& x) -> FixedIntervalVector<1> {
        FixedIntervalVector<1> result;
        result[0] = x[0] * x[0];
        return result;
    };
    
    // 定义雅可比矩阵
    auto jacobian = [](const FixedIntervalVector<1>& x) -> FixedIntervalMatrix<1, 1> {
        FixedIntervalMatrix<1, 1> J;
        J(0, 0) = x[0] * PSGMDirectedInterval(2.0);
        return J;
    };
    
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
    double delta = 0.1;
    FixedIntervalVector<1> X0;
    X0[0] = PSGMDirectedInterval(x_star - delta, x_star + delta);
    
    std::cout << "   初始区间X0: " << X0 << std::endl;
    
    // 3. 运行广义Krawczyk迭代
    std::cout << "\n3. 运行广义Krawczyk迭代...\n";
    FixedGeneralizedKrawczykSolver<1> solver(f, jacobian, 1e-8, 50);
    
    auto result = solver.solve(X0);
    
    if (result.success)
    {
        std::cout << "   ✓ Krawczyk迭代成功!\n";
        std::cout << "   迭代次数: " << result.iterations << std::endl;
        std::cout << "   最终区间X_final: " << result.solution << std::endl;
        
        // 4. 内包围验证
        std::cout << "\n4. 内包围验证...\n";
        FixedIntervalVector<1> fX = f(result.solution);
        std::cout << "   f(X_final) = " << fX << std::endl;
        std::cout << "   Y = " << Y << std::endl;
        
        bool inner_inclusion = verify_inner_inclusion(result.solution, Y, f);
        std::cout << "   内包围验证结果: " << (inner_inclusion ? "✓ 有效" : "✗ 无效") << std::endl;
    }
    else
    {
        std::cout << "   ✗ Krawczyk迭代失败!\n";
        std::cout << "   消息: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Tolerance Embedding 固定维数实现 ===\n\n";
    
    // 展示Kaucher区间算术的基本特性
    std::cout << "=== Kaucher区间算术特性 ===\n";
    PSGMDirectedInterval proper(1.0, 2.0);
    PSGMDirectedInterval improper(2.0, 1.0);
    
    std::cout << "正常区间 a = " << proper << std::endl;
    std::cout << "非正常区间 b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "a的逆元 = " << (-proper) << std::endl;
    std::cout << "b的对偶 = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // 运行测试
    testToleranceEmbedding1D(); // 1D测试
    testToleranceEmbedding2D();   // 2D测试
    
    std::cout << "=== 测试完成 ===\n";
    
    return 0;
}