#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <Eigen/Dense>
#include "PSGMDirectedInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"

// Newton method pre-iteration function to find high-precision center point
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
    
    std::cout << "   f(X) = " << fX << std::endl;
    std::cout << "   Y    = " << Y << std::endl;
    
    // 检查fX是否严格包含在Y的对偶区间内
    // 对于Y是improper区间的情况，其对偶区间是proper区间[Y.upper(), Y.lower()]
    for (size_t i = 0; i < N; ++i)
    {
        if (!Y[i].isImproper())
        {
            std::cout << "   Error: Y[" << i << "] is not an improper interval" << std::endl;
            return false;
        }
        
        PSGMDirectedInterval Y_dual = Y[i].dual(); // 转换为proper区间 [Y.upper(), Y.lower()]
        
        if (!fX[i].isProper())
        {
            std::cout << "   Error: f(X)[" << i << "] is not a proper interval" << std::endl;
            return false;
        }
        
        // 检查fX(i)是否严格包含在Y_dual内部
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
    
    // 1. 使用牛顿法预迭代找到高精度中心点x*
    std::cout << "1. Newton method pre-iteration to find center point...\n";
    double x1_init = 1.5;
    double x2_init = 1.5;
    auto [x1_star, x2_star] = newton_method_2d(x1_init, x2_init);
    
    std::cout << "   初始猜测: (" << x1_init << ", " << x2_init << ")" << std::endl;
    std::cout << "   Newton method result: x* = (" << x1_star << ", " << x2_star << ")" << std::endl;
    
    // 2. 构造小初始区间
    std::cout << "\n2. Constructing small initial interval...\n";
    double delta = 0.01; // 使用更小的初始区间半径
    FixedIntervalVector<2> X0;
    X0[0] = PSGMDirectedInterval(x1_star - delta, x1_star + delta);
    X0[1] = PSGMDirectedInterval(x2_star - delta, x2_star + delta);
    
    std::cout << "   初始区间X0: " << X0 << std::endl;
    
    // 3. 定义函数f和雅可比矩阵
    std::cout << "\n3. Defining function f and Jacobian matrix...\n";
    
    // 定义原函数f
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
    
    // 4. 运行广义Krawczyk迭代
    std::cout << "\n4. Running generalized Krawczyk iteration...\n";
    FixedGeneralizedKrawczykSolver<2> solver(f, jacobian, 1e-8, 50);
    
    auto result = solver.solve(X0);
    
    if (result.success)
    {
        std::cout << "   ✓ Krawczyk iteration succeeded!\n";
        std::cout << "   Iterations: " << result.iterations << std::endl;
        std::cout << "   Solution interval X_final: " << result.solution << std::endl;
        std::cout << "   最终宽度: " << result.finalWidth << std::endl;
        
        // 5. 定义目标Y（improper区间）
        std::cout << "\n5. Defining target Y (improper interval)...\n";
        FixedIntervalVector<2> Y;
        Y[0] = PSGMDirectedInterval(3.9, 4.1).dual(); // 非正常区间 [4.1, 3.9]
        Y[1] = PSGMDirectedInterval(-0.1, 0.1).dual(); // 非正常区间 [0.1, -0.1]
        
        std::cout << "   Y = " << Y << std::endl;
        
        // 6. 内包围验证
        std::cout << "\n6. Inner enclosure verification...\n";
        bool inner_inclusion = verify_inner_inclusion(result.solution, Y, f);
        
        std::cout << "\n7. Final result: " << (inner_inclusion ? "✓ Inner enclosure verified" : "✗ Inner enclosure failed") << std::endl;
        
        // 7. 额外验证：检查四个角点
        std::cout << "\n8. Corner point verification...\n";
        
        // 定义点函数f
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
        std::cout << "   ✗ Krawczyk iteration failed!\n";
        std::cout << "   消息: " << result.message << std::endl;
        std::cout << "   最终区间: " << result.solution << std::endl;
        std::cout << "   最终宽度: " << result.finalWidth << std::endl;
    }
    
    std::cout << "\n=== Test Completed ===\n";
    
    return 0;
}