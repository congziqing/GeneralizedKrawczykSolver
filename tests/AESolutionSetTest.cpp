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

// Helper function: create dual interval (Kaucher dual)
KaucherInterval dualInterval(const KaucherInterval& interval)
{
    return KaucherInterval(interval.upper(), interval.lower());
}

// Test 1: AE-Solution Set Test - x² - p = 0, p ∈ [4, 9]
void testAESolutionSet()
{
    std::cout << "=== Test 1: AE-Solution Set Test (x² - p = 0, p ∈ [4, 9]) ===\n\n";
    
    // 定义方程 f(x, p) = x² - p = 0
    // 这里p是参数，使用对偶形式表示全称量词
    auto f = [](const ik::IntervalVector<1>& x, const KaucherInterval& p) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        result[0] = x[0] * x[0] - p;
        return result;
    };
    
    // 计算关于x的Jacobian
    auto jacobian_x = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = 2x
        J(0, 0) = midPointInterval(2.0) * x[0];
        return J;
    };
    
    // 参数p的正常区间：[4, 9]
    KaucherInterval p_normal(4.0, 9.0);
    // 参数p的对偶区间（表示全称量词）：[9, 4]
    KaucherInterval p_dual = dualInterval(p_normal);
    
    std::cout << "Parameter p proper interval: " << p_normal << std::endl;
    std::cout << "Parameter p dual interval (universal quantifier): " << p_dual << std::endl;
    
    // 初始区间X：[0, 10]
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(0.0, 10.0);
    
    int maxIterations = 20;
    double tolerance = 1e-8;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // 1. 计算中点
        double c = xn[0].middle();
        
        // 2. 在中点处计算 f
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec, p_dual);
        
        // 3. 计算 Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian_x(xn);
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
        std::cout << "  c = " << c << std::endl;
        std::cout << "  J = " << J(0, 0) << std::endl;
        std::cout << "  p_dual = " << p_dual << std::endl;
        std::cout << "  f(c) = " << f_c[0] << std::endl;
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
        
        // 如果区间为空或宽度增加，尝试使用 Kaucher 算术
        if (x_new[0].isEmpty() || x_new[0].width() > xn[0].width())
        {
            std::cout << "  Proper interval method encountered issues, using Kaucher arithmetic...\n";
            // 使用 Krawczyk 区间作为新的区间，允许非正常区间
            xn = K;
        }
        std::cout << std::endl;
    }
}

// Test 2: Complex AE-Solution Set Test - x - p + x³ = 0, p ∈ [-0.1, 0.1]
void testComplexAESolutionSet()
{
    std::cout << "\n=== Test 2: Complex AE-Solution Set Test (x - p + x³ = 0, p ∈ [-0.1, 0.1]) ===\n\n";
    
    // 定义方程 f(x, p) = x - p + x³ = 0
    // 这里p是参数，使用对偶形式表示全称量词
    auto f = [](const ik::IntervalVector<1>& x, const KaucherInterval& p) -> ik::IntervalVector<1>
    {
        ik::IntervalVector<1> result;
        result[0] = x[0] - p + x[0] * x[0] * x[0];
        return result;
    };
    
    // 计算关于x的Jacobian
    auto jacobian_x = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1>
    {
        ik::IntervalMatrix<1, 1> J;
        // J(x) = 1 + 3x²
        J(0, 0) = midPointInterval(1.0) + midPointInterval(3.0) * x[0] * x[0];
        return J;
    };
    
    // 参数p的正常区间：[-0.1, 0.1]
    KaucherInterval p_normal(-0.1, 0.1);
    // 参数p的对偶区间（表示全称量词）：[0.1, -0.1]
    KaucherInterval p_dual = dualInterval(p_normal);
    
    std::cout << "参数p的正常区间：" << p_normal << std::endl;
    std::cout << "参数p的对偶区间（全称量词）：" << p_dual << std::endl;
    
    // 初始区间X：[-1, 1]
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(-1.0, 1.0);
    
    int maxIterations = 20;
    double tolerance = 1e-8;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // 1. 计算中点
        double c = xn[0].middle();
        
        // 2. 在中点处计算 f
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec, p_dual);
        
        // 3. 计算 Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian_x(xn);
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
        std::cout << "  p_dual = " << p_dual << std::endl;
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
    std::cout << "=== Kaucher AE-Solution Set Test ===\n\n";
    
    // 运行测试
    testAESolutionSet();
    testComplexAESolutionSet();
    
    std::cout << "\n=== Tests Completed ===\n";
    return 0;
}
