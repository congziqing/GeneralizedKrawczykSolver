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

// Test 1: Algebraic difference test X + [2, 3] = [5, 6]
void testAlgebraicDifference()
{
    std::cout << "=== Test 1: Algebraic Difference Test (X + [2, 3] = [5, 6]) ===\n\n";
    
    // 定义方程 f(X) = X + [2, 3] - [5, 6] = 0
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        // X + [2, 3] - [5, 6]
        result[0] = x[0] + KaucherInterval(2.0, 3.0) - KaucherInterval(5.0, 6.0);
        return result;
    };
    
    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        J(0, 0) = midPointInterval(1.0);
        return J;
    };
    
    // 使用 ik::IntervalVector 实现简单的广义 Krawczyk 迭代
    std::cout << "初始区间: [0, 10]\n";
    
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(0.0, 10.0);
    
    int maxIterations = 10;
    double tolerance = 1e-8;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // 1. 计算中点
        double c = xn[0].middle();
        
        // 2. 在中点处计算 f
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec);
        
        // 3. 计算 Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian(xn);
        Eigen::MatrixXd J_mid = J.midpoint();
        Eigen::MatrixXd Y = J_mid.inverse();
        
        // 4. 计算 term1: c - Y * f(c)
        double term1_val = c - Y(0, 0) * f_c[0].middle();
        ik::IntervalVector<1> term1;
        term1[0] = midPointInterval(term1_val);
        
        // 5. 计算 I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(1, 1);
        Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
        ik::IntervalMatrix<1, 1> I_minus_YJ_interval;
        I_minus_YJ_interval(0, 0) = midPointInterval(I_minus_YJ(0, 0));
        
        // 6. 计算 box - c
        double delta = (xn[0].upper() - xn[0].lower()) / 2.0;
        ik::IntervalVector<1> box_minus_c;
        box_minus_c[0] = KaucherInterval(-delta, delta);
        
        // 7. 计算 term3
        ik::IntervalVector<1> term3 = I_minus_YJ_interval * box_minus_c;
        
        // 8. 计算 K = term1 + term3
        ik::IntervalVector<1> K = term1 + term3;
        
        // 9. 计算 X_new = X ∩ K
        ik::IntervalVector<1> x_new = xn.intersection(K);
        
        std::cout << "Iteration " << iter + 1 << ":\n";
        std::cout << "  xn = " << xn << std::endl;
        std::cout << "  K = " << K << std::endl;
        std::cout << "  x_new = " << x_new << std::endl;
        std::cout << "  Width = " << x_new[0].width() << std::endl;
        
        // 检查收敛
        if (x_new[0].width() < tolerance)
        {
            std::cout << "\n✓ Converged to " << x_new << ", Iterations: " << iter + 1 << std::endl;
            break;
        }
        
        xn = x_new;
        
        // 如果区间为空，尝试使用非正常区间
        if (x_new[0].isEmpty())
        {
            std::cout << "  Proper interval intersection is empty, trying improper interval...\n";
            // 使用 Krawczyk 区间作为新的区间，允许非正常区间
            xn = K;
        }
        std::cout << std::endl;
    }
}

// Test 2: Reverse width iteration f(x) = x - tanh(2x) = 0
void testReverseWidthIteration()
{
    std::cout << "\n=== Test 2: Reverse Width Iteration f(x) = x - tanh(2x) = 0 ===\n\n";
    
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        // f(x) = x - tanh(2x)
        result[0] = x[0] - ::interval::tanh(x[0] * 2.0);
        return result;
    };
    
    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = 1 - 2 * (1 - tanh^2(2x))
        // 注意：tanh(2x) 的导数是 2(1 - tanh^2(2x))
        KaucherInterval tanh_2x = ::interval::tanh(x[0] * 2.0);
        J(0, 0) = midPointInterval(1.0) - midPointInterval(2.0) * (midPointInterval(1.0) - tanh_2x * tanh_2x);
        return J;
    };
    
    // 从非常大的初始区间开始
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(-10.0, 10.0);
    
    int maxIterations = 20;
    double tolerance = 1e-8;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // 1. 计算中点
        double c = xn[0].middle();
        
        // 2. 在中点处计算 f
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec);
        
        // 3. 计算 Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian(xn);
        Eigen::MatrixXd J_mid = J.midpoint();
        
        // 检查 Jacobian 是否可逆
        double det = J_mid(0, 0);
        if (std::abs(det) < 1e-12)
        {
            std::cout << "  Jacobian is singular, using approximate inverse...\n";
            det = (det > 0) ? 1e-12 : -1e-12;
        }
        Eigen::MatrixXd Y(1, 1);
        Y(0, 0) = 1.0 / det;
        
        // 4. 计算 term1: c - Y * f(c)
        double term1_val = c - Y(0, 0) * f_c[0].middle();
        ik::IntervalVector<1> term1;
        term1[0] = midPointInterval(term1_val);
        
        // 5. 计算 I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(1, 1);
        Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
        ik::IntervalMatrix<1, 1> I_minus_YJ_interval;
        I_minus_YJ_interval(0, 0) = midPointInterval(I_minus_YJ(0, 0));
        
        // 6. 计算 box - c
        double delta = (xn[0].upper() - xn[0].lower()) / 2.0;
        ik::IntervalVector<1> box_minus_c;
        box_minus_c[0] = KaucherInterval(-delta, delta);
        
        // 7. 计算 term3
        ik::IntervalVector<1> term3 = I_minus_YJ_interval * box_minus_c;
        
        // 8. 计算 K = term1 + term3
        ik::IntervalVector<1> K = term1 + term3;
        
        // 9. 计算 X_new = X ∩ K
        ik::IntervalVector<1> x_new = xn.intersection(K);
        
        std::cout << "迭代 " << iter + 1 << ":\n";
        std::cout << "  xn = " << xn << std::endl;
        std::cout << "  J = " << J(0, 0) << std::endl;
        std::cout << "  K = " << K << std::endl;
        std::cout << "  x_new = " << x_new << std::endl;
        std::cout << "  宽度 = " << x_new[0].width() << std::endl;
        
        // 检查收敛
        if (x_new[0].width() < tolerance)
        {
            std::cout << "\n✓ 收敛到 " << x_new << "，迭代次数: " << iter + 1 << std::endl;
            break;
        }
        
        xn = x_new;
        
        // 如果区间为空或宽度增加，尝试使用非正常区间
        if (x_new[0].isEmpty() || x_new[0].width() > xn[0].width())
        {
            std::cout << "  Proper interval method encountered issues, using Kaucher arithmetic...\n";
            // 使用 Krawczyk 区间作为新的区间，允许非正常区间
            xn = K;
        }
        std::cout << std::endl;
    }
}

int main()
{
    std::cout << "=== Kaucher Interval Arithmetic Properties Test ===\n\n";
    
    // Show basic properties of Kaucher interval arithmetic
    std::cout << "Kaucher Interval Arithmetic Basic Properties：\n";
    KaucherInterval proper(1.0, 2.0);
    KaucherInterval improper(2.0, 1.0);
    KaucherInterval a(2.0, 3.0);
    KaucherInterval b(5.0, 6.0);
    
    std::cout << "Proper interval: " << proper << std::endl;
    std::cout << "Improper interval: " << improper << std::endl;
    std::cout << "a = " << a << ", b = " << b << std::endl;
    std::cout << "a - b = " << (a - b) << " (经典差)\n";
    std::cout << "a.algebraicDiff(b) = " << a.algebraicDiff(b) << " (代数差)\n";
    std::cout << "a + (-b) = " << (a + (-b)) << " (使用逆元)\n\n";
    
    // 运行测试
    testAlgebraicDifference();
    testReverseWidthIteration();
    
    std::cout << "\n=== Tests Completed ===\n";
    return 0;
}