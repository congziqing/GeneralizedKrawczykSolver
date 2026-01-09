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

// Test 1: Simple linear expanding map x = 2x - [2, 4]
void testLinearExpandingMap()
{
    std::cout << "=== Test 1: Simple Linear Expanding Map (x = 2x - [2, 4]) ===\n\n";
    
    // 定义方程 f(x) = x - (2x - [2, 4]) = -x + [2, 4] = 0
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        // x - (2x - [2, 4]) = -x + [2, 4]
        result[0] = x[0] - (midPointInterval(2.0) * x[0] - KaucherInterval(2.0, 4.0));
        return result;
    };
    
    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = -1
        J(0, 0) = midPointInterval(-1.0);
        return J;
    };
    
    // 从初始区间 [2, 4] 开始
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(2.0, 4.0);
    
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
        
        std::cout << "Iteration " << iter + 1 << ":\n";
        std::cout << "  xn = " << xn << std::endl;
        std::cout << "  J = " << J(0, 0) << std::endl;
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

// Test 2: Quadratic expanding map f(x) = x² - x - [0, 2] = 0, finding solution near x=2
void testQuadraticExpandingMap()
{
    std::cout << "\n=== Test 2: Quadratic Expanding Map (f(x) = x² - x - [0, 2] = 0) ===\n\n";
    
    // 定义方程 f(x) = x² - x - [0, 2]
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        result[0] = x[0] * x[0] - x[0] - KaucherInterval(0.0, 2.0);
        return result;
    };
    
    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = 2x - 1
        J(0, 0) = midPointInterval(2.0) * x[0] - midPointInterval(1.0);
        return J;
    };
    
    // 从包含x=2的初始区间 [1, 3] 开始
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(1.0, 3.0);
    
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
            std::cout << "  Jacobian 奇异，使用近似逆...\n";
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
        std::cout << "  c = " << c << std::endl;
        std::cout << "  J = " << J(0, 0) << std::endl;
        std::cout << "  f(c) = " << f_c[0] << std::endl;
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
        
        // 如果区间为空或宽度增加，尝试使用 Kaucher 算术
        if (x_new[0].isEmpty() || x_new[0].width() > xn[0].width())
        {
            std::cout << "  正常区间方法遇到问题，使用 Kaucher 算术...\n";
            // 使用 Krawczyk 区间作为新的区间，允许非正常区间
            xn = K;
        }
        std::cout << std::endl;
    }
}

int main()
{
    std::cout << "=== Kaucher Expanding Map Fixed Point Test ===\n\n";
    
    // 运行测试
    testLinearExpandingMap();
    testQuadraticExpandingMap();
    
    std::cout << "\n=== Tests Completed ===\n";
    return 0;
}
