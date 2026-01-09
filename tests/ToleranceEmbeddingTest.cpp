#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"

// Helper function: create midpoint interval
KaucherInterval midPointInterval(double x)
{
    return PSGMDirectedInterval(x, x);
}

// Helper function: create dual interval (Kaucher dual)
PSGMDirectedInterval dualInterval(const PSGMDirectedInterval& interval)
{
    return PSGMDirectedInterval(interval.upper(), interval.lower());
}

// 定义系统方程 f(x) = Y
// f(x) = [x1 + exp(x1 - x2), x1² + x2]^T
void toleranceEmbeddingTest()
{
    std::cout << "=== Tolerance Embedding Problem Solving (Generalized Krawczyk Iteration) ===\n\n";
    
    // 已知目标输出范围 Y = ([2.9, 3.1], [4.9, 5.1])^T
    // 为了求解内包围，需要将Y转换为对偶形式（Improper Interval）
    PSGMDirectedInterval y1_normal(2.9, 3.1);
    PSGMDirectedInterval y2_normal(4.9, 5.1);
    
    // 转换为对偶形式，触发内包围求解模式
    PSGMDirectedInterval y1_dual = dualInterval(y1_normal);
    PSGMDirectedInterval y2_dual = dualInterval(y2_normal);
    
    std::cout << "Target output range Y (proper interval):\n";
    std::cout << "  y1 = " << y1_normal << std::endl;
    std::cout << "  y2 = " << y2_normal << std::endl;
    
    std::cout << "\nTarget output range Y (dual form, Improper Interval):\n";
    std::cout << "  y1_dual = " << y1_dual << std::endl;
    std::cout << "  y2_dual = " << y2_dual << std::endl;
    
    // 初始区间 X0 = ([0, 5], [0, 5])^T
    ik::IntervalVector<2> xn;
    xn[0] = PSGMDirectedInterval(0.0, 5.0);
    xn[1] = PSGMDirectedInterval(0.0, 5.0);
    
    std::cout << "\nInitial interval X0:\n";
    std::cout << "  x1 = " << xn[0] << std::endl;
    std::cout << "  x2 = " << xn[1] << std::endl;
    
    int maxIterations = 20;
    double tolerance = 1e-4;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        std::cout << "\nIteration " << iter + 1 << ":\n";
        std::cout << "  Current interval xn = " << xn << std::endl;
        
        // 1. 计算中点
        double c1 = xn[0].middle();
        double c2 = xn[1].middle();
        std::cout << "  Midpoint c = (" << c1 << ", " << c2 << ")" << std::endl;
        
        // 2. 在中点处计算函数值和雅可比矩阵
        double dx1 = 0.01;
        double dx2 = 0.01;
        
        // 计算函数值 f(c)
        double f1_c = c1 + std::exp(c1 - c2);
        double f2_c = c1 * c1 + c2;
        
        // 计算雅可比矩阵的近似值（使用有限差分）
        double f1_c_dx1 = (c1 + dx1) + std::exp((c1 + dx1) - c2);
        double f1_c_dx2 = c1 + std::exp(c1 - (c2 + dx2));
        double f2_c_dx1 = (c1 + dx1) * (c1 + dx1) + c2;
        double f2_c_dx2 = c1 * c1 + (c2 + dx2);
        
        double J11 = (f1_c_dx1 - f1_c) / dx1;
        double J12 = (f1_c_dx2 - f1_c) / dx2;
        double J21 = (f2_c_dx1 - f2_c) / dx1;
        double J22 = (f2_c_dx2 - f2_c) / dx2;
        
        std::cout << "  Function value at midpoint: f(c) = (" << f1_c << ", " << f2_c << ")" << std::endl;
        std::cout << "  Jacobian matrix J = [[" << J11 << ", " << J12 << "], [" << J21 << ", " << J22 << "]]" << std::endl;
        
        // 3. 构建雅可比矩阵的中点矩阵
        Eigen::MatrixXd J_mid(2, 2);
        J_mid << J11, J12, J21, J22;
        
        // 4. 计算逆矩阵 Y
        Eigen::MatrixXd Y = J_mid.inverse();
        std::cout << "  Inverse matrix Y = \n" << Y << std::endl;
        
        // 5. 计算 term1: c - Y * (f(c) - Y_dual)
        // 注意：这里使用对偶区间 Y_dual
        Eigen::VectorXd c_eigen(2);
        c_eigen << c1, c2;
        
        Eigen::VectorXd f_c_eigen(2);
        f_c_eigen << f1_c, f2_c;
        
        Eigen::VectorXd Y_dual_eigen(2);
        Y_dual_eigen << y1_dual.middle(), y2_dual.middle();
        
        Eigen::VectorXd term1_eigen = c_eigen - Y * (f_c_eigen - Y_dual_eigen);
        std::cout << "  term1 = (" << term1_eigen(0) << ", " << term1_eigen(1) << ")" << std::endl;
        
        // 6. 计算 I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
        Eigen::MatrixXd I_minus_YJ = I - Y * J_mid;
        std::cout << "  I - YJ = \n" << I_minus_YJ << std::endl;
        
        // 7. 计算 box - c
        double delta1 = (xn[0].upper() - xn[0].lower()) / 2.0;
        double delta2 = (xn[1].upper() - xn[1].lower()) / 2.0;
        std::cout << "  box - c = ([" << -delta1 << ", " << delta1 << "], [" << -delta2 << ", " << delta2 << "])" << std::endl;
        
        // 8. 计算 term3 = (I - YJ) * (box - c)
        // 这里我们直接计算区间形式
        double term3_1_lower = I_minus_YJ(0, 0) * (-delta1) + I_minus_YJ(0, 1) * (-delta2);
        double term3_1_upper = I_minus_YJ(0, 0) * delta1 + I_minus_YJ(0, 1) * delta2;
        double term3_2_lower = I_minus_YJ(1, 0) * (-delta1) + I_minus_YJ(1, 1) * (-delta2);
        double term3_2_upper = I_minus_YJ(1, 0) * delta1 + I_minus_YJ(1, 1) * delta2;
        
        // 确保区间是正确的（lower <= upper）
        if (term3_1_lower > term3_1_upper) std::swap(term3_1_lower, term3_1_upper);
        if (term3_2_lower > term3_2_upper) std::swap(term3_2_lower, term3_2_upper);
        
        std::cout << "  term3 = ([" << term3_1_lower << ", " << term3_1_upper << "], [" << term3_2_lower << ", " << term3_2_upper << "])" << std::endl;
        
        // 9. 计算 Krawczyk 区间 K = term1 + term3
        PSGMDirectedInterval K1(term1_eigen(0) + term3_1_lower, term1_eigen(0) + term3_1_upper);
        PSGMDirectedInterval K2(term1_eigen(1) + term3_2_lower, term1_eigen(1) + term3_2_upper);
        
        ik::IntervalVector<2> K;
        K[0] = K1;
        K[1] = K2;
        
        std::cout << "  Krawczyk interval K = " << K << std::endl;
        
        // 10. 计算 X_new = X ∩ K
        ik::IntervalVector<2> x_new;
        x_new[0] = xn[0].meet(K[0]);
        x_new[1] = xn[1].meet(K[1]);
        
        std::cout << "  New interval x_new = " << x_new << std::endl;
        
        // 11. 确保新区间不为空
        if (x_new[0].isEmpty() || x_new[1].isEmpty())
        {
            std::cout << "  New interval is empty, using Krawczyk interval as new interval\n";
            x_new = K;
        }
        
        // 12. 计算区间宽度
        double width1 = x_new[0].width();
        double width2 = x_new[1].width();
        double max_width = std::max(width1, width2);
        std::cout << "  Interval widths: x1 = " << width1 << ", x2 = " << width2 << ", max = " << max_width << std::endl;
        
        // 13. 更新区间
        xn = x_new;
        
        // 14. 检查收敛
        if (max_width < tolerance)
        {
            std::cout << "\n✓ Converged to " << xn << ", Iterations: " << iter + 1 << std::endl;
            
            // 验证结果
            std::cout << "\n=== Result Verification ===\n";
            // 在收敛区间内取几个点进行验证
            double x1_test1 = xn[0].lower();
            double x2_test1 = xn[1].lower();
            double x1_test2 = xn[0].upper();
            double x2_test2 = xn[1].upper();
            double x1_test3 = xn[0].middle();
            double x2_test3 = xn[1].middle();
            
            // 计算函数值
            auto compute_f = [](double x1, double x2) -> std::pair<double, double>
            {
                double f1 = x1 + std::exp(x1 - x2);
                double f2 = x1 * x1 + x2;
                return {f1, f2};
            };
            
            auto [f1_val1, f2_val1] = compute_f(x1_test1, x2_test1);
            auto [f1_val2, f2_val2] = compute_f(x1_test2, x2_test2);
            auto [f1_val3, f2_val3] = compute_f(x1_test3, x2_test3);
            
            std::cout << "  At point x = (" << x1_test1 << ", " << x2_test1 << "):\n";
            std::cout << "    f1(x) = " << f1_val1 << ", should be in " << y1_normal << " range" << std::endl;
            std::cout << "    f2(x) = " << f2_val1 << ", should be in " << y2_normal << " range" << std::endl;
            
            std::cout << "  At point x = (" << x1_test2 << ", " << x2_test2 << "):\n";
            std::cout << "    f1(x) = " << f1_val2 << ", should be in " << y1_normal << " range" << std::endl;
            std::cout << "    f2(x) = " << f2_val2 << ", should be in " << y2_normal << " range" << std::endl;
            
            std::cout << "  At midpoint x = (" << x1_test3 << ", " << x2_test3 << "):\n";
            std::cout << "    f1(x) = " << f1_val3 << ", should be in " << y1_normal << " range" << std::endl;
            std::cout << "    f2(x) = " << f2_val3 << ", should be in " << y2_normal << " range" << std::endl;
            
            // 检查是否在目标范围内
            bool f1_in_range1 = (f1_val1 >= y1_normal.lower()) && (f1_val1 <= y1_normal.upper());
            bool f2_in_range1 = (f2_val1 >= y2_normal.lower()) && (f2_val1 <= y2_normal.upper());
            bool f1_in_range2 = (f1_val2 >= y1_normal.lower()) && (f1_val2 <= y1_normal.upper());
            bool f2_in_range2 = (f2_val2 >= y2_normal.lower()) && (f2_val2 <= y2_normal.upper());
            bool f1_in_range3 = (f1_val3 >= y1_normal.lower()) && (f1_val3 <= y1_normal.upper());
            bool f2_in_range3 = (f2_val3 >= y2_normal.lower()) && (f2_val3 <= y2_normal.upper());
            
            bool all_valid = f1_in_range1 && f2_in_range1 && f1_in_range2 && f2_in_range2 && f1_in_range3 && f2_in_range3;
            
            if (all_valid)
            {
                std::cout << "  ✓ All test points' function values are within target range\n";
            }
            else
            {
                std::cout << "  ✗ Some test points' function values are not within target range\n";
            }
            
            // 输出最终内包围解
            std::cout << "\n=== Final Inner Enclosure Solution ===\n";
            std::cout << "X_in = " << xn << std::endl;
            break;
        }
    }
}

int main()
{
    toleranceEmbeddingTest();
    return 0;
}
