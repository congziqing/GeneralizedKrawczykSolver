#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"
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

// 固定维数广义Krawczyk求解器
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

// 验证函数：检查f(X)是否严格位于Y的内部
template<size_t N>
bool verify_inner_inclusion(const FixedIntervalVector<N>& X, 
                           const FixedIntervalVector<N>& Y)
{
    // 计算f(X)
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
    
    // 检查fX是否严格包含在Y的对偶区间内
    // 对于Y是improper区间的情况，其对偶区间是proper区间[Y.upper(), Y.lower()]
    for (size_t i = 0; i < N; ++i)
    {
        if (!Y[i].isImproper())
        {
            std::cout << "   错误：Y[" << i << "] 不是improper区间" << std::endl;
            return false;
        }
        
        ik::KaucherInterval Y_dual = Y[i].dual(); // Convert to proper interval [Y.upper(), Y.lower()]
        
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

// 同伦连续法求解Tolerance Embedding Problem
void testToleranceEmbeddingHomotopy()
{
    std::cout << "=== 测试: Tolerance Embedding Problem (2D) - 同伦连续法 ===\n";
    
    // 定义问题：
    // f1(x1, x2) = x1^2 + x2^2 = 4 + e1
    // f2(x1, x2) = x1 - x2 = 0 + e2
    // 其中 e1 ∈ [-0.1, 0.1], e2 ∈ [-0.1, 0.1]
    
    // 目标Y：使用improper区间表示对偶形式
    FixedIntervalVector<2> Y;
    Y[0] = ik::KaucherInterval(3.9, 4.1).dual(); // Improper interval [4.1, 3.9]
    Y[1] = ik::KaucherInterval(-0.1, 0.1).dual(); // Improper interval [0.1, -0.1]
    
    std::cout << "目标Y（improper）: " << Y << std::endl;
    
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
    
    // 1. 构造辅助系统g(x) = x - x_start = 0
    std::cout << "\n1. 构造辅助系统和同伦函数...\n";
    std::pair<double, double> x_start = {2.5, 2.5}; // 选择的中心点
    std::cout << "   辅助系统中心点: (" << x_start.first << ", " << x_start.second << ")" << std::endl;
    
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
    
    // 3. 将t从0逐步增加到1
    std::cout << "\n2. 同伦连续过程 (t从0到1，分10步)...\n";
    int steps = 10;
    FixedIntervalVector<2> currentBox;
    bool converged = true;
    
    for (int step = 0; step <= steps; ++step)
    {
        double t = static_cast<double>(step) / steps;
        std::cout << "\n   步骤 " << step << "/" << steps << " (t = " << t << "):\n";
        
        // 构造当前t值下的同伦函数和雅可比
        auto H = createHomotopyFunction(t);
        auto H_jacobian = createHomotopyJacobian(t);
        
        // 构造初始区间
        FixedIntervalVector<2> initialBox;
        if (step == 0)
        {
            // 第一步t=0，使用辅助系统，区间以x_start为中心
            double delta = 0.1; // 初始区间半径
            initialBox[0] = ik::KaucherInterval(x_start.first - delta, x_start.first + delta);
            initialBox[1] = ik::KaucherInterval(x_start.second - delta, x_start.second + delta);
        }
        else
        {
            // 后续步骤，使用上一步的解作为中心
            double delta = 0.1; // 区间半径
            auto c = currentBox.midpoint();
            initialBox[0] = ik::KaucherInterval(c(0) - delta, c(0) + delta);
            initialBox[1] = ik::KaucherInterval(c(1) - delta, c(1) + delta);
        }
        
        std::cout << "   初始区间: " << initialBox << std::endl;
        
        // 运行Krawczyk迭代
        FixedGeneralizedKrawczykSolver<2> solver(H, H_jacobian, 1e-8, 50);
        auto result = solver.solve(initialBox);
        
        if (result.success)
        {
            std::cout << "   ✓ Krawczyk迭代成功!\n";
            std::cout << "   迭代次数: " << result.iterations << std::endl;
            std::cout << "   解区间: " << result.solution << std::endl;
            std::cout << "   区间宽度: " << result.finalWidth << std::endl;
            currentBox = result.solution;
        }
        else
        {
            std::cout << "   ✗ Krawczyk迭代失败!\n";
            std::cout << "   消息: " << result.message << std::endl;
            std::cout << "   最终区间: " << result.solution << std::endl;
            converged = false;
            break;
        }
    }
    
    if (converged)
    {
        std::cout << "\n3. 同伦连续过程成功完成!\n";
        std::cout << "   最终解区间: " << currentBox << std::endl;
        
        // 4. 验证内包围
        std::cout << "\n4. 内包围验证...\n";
        bool inner_inclusion = verify_inner_inclusion(currentBox, Y);
        
        std::cout << "\n5. 最终结果: " << (inner_inclusion ? "✓ 内包围验证通过" : "✗ 内包围验证失败") << std::endl;
        
        // 5. 额外验证：检查四个角点
        std::cout << "\n6. 角点验证...\n";
        std::array<std::pair<double, double>, 4> corners;
        corners[0] = std::make_pair(currentBox[0].lower(), currentBox[1].lower());
        corners[1] = std::make_pair(currentBox[0].lower(), currentBox[1].upper());
        corners[2] = std::make_pair(currentBox[0].upper(), currentBox[1].lower());
        corners[3] = std::make_pair(currentBox[0].upper(), currentBox[1].upper());
        
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
        std::cout << "\n3. 同伦连续过程失败!\n";
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Tolerance Embedding 同伦连续实现 ===\n\n";
    
    // 展示Kaucher区间算术的基本特性
    std::cout << "=== Kaucher区间算术特性 ===\n";
    ik::KaucherInterval proper(1.0, 2.0);
    ik::KaucherInterval improper(2.0, 1.0);
    
    std::cout << "正常区间 a = " << proper << std::endl;
    std::cout << "非正常区间 b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "a的逆元 = " << (-proper) << std::endl;
    std::cout << "b的对偶 = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // 运行同伦连续测试
    testToleranceEmbeddingHomotopy();
    
    std::cout << "=== 测试完成 ===\n";
    
    return 0;
}