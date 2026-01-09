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
    return KaucherInterval(x, x);
}

// Helper function: create dual interval (Kaucher dual)
KaucherInterval dualInterval(const KaucherInterval& interval)
{
    return KaucherInterval(interval.upper(), interval.lower());
}

// 定义系统方程 f(x) = Y
// f(x) = [x1 + exp(x1 - x2), x1² + x2]^T
void toleranceEmbeddingCorrect()
{
    std::cout << "=== Tolerance Embedding Problem Solving (Correct Generalized Krawczyk Iteration) ===\n\n";
    
    // 已知目标输出范围 Y = ([2.9, 3.1], [4.9, 5.1])^T
    // 为了求解内包围，需要将Y转换为对偶形式（Improper Interval）
    KaucherInterval y1_normal(2.9, 3.1);
    KaucherInterval y2_normal(4.9, 5.1);
    
    // Convert to dual form, triggering inner enclosure solving mode
    KaucherInterval y1_dual = dualInterval(y1_normal);
    KaucherInterval y2_dual = dualInterval(y2_normal);
    
    std::cout << "Target output range Y (proper interval):\n";
    std::cout << "  y1 = " << y1_normal << std::endl;
    std::cout << "  y2 = " << y2_normal << std::endl;
    
    std::cout << "\nTarget output range Y (dual form, Improper Interval):\n";
    std::cout << "  y1_dual = " << y1_dual << std::endl;
    std::cout << "  y2_dual = " << y2_dual << std::endl;
    
    // 定义系统方程 f(x) = Y
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2>
    {
        ik::IntervalVector<2> result;
        result[0] = x[0] + ::interval::exp(x[0] - x[1]);
        result[1] = x[0] * x[0] + x[1];
        return result;
    };
    
    // 计算区间雅可比矩阵 J(X)
    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2>
    {
        ik::IntervalMatrix<2, 2> J;
        
        // 计算偏导数
        // J[0][0] = 1 + exp(x1 - x2)
        J(0, 0) = midPointInterval(1.0) + ::interval::exp(x[0] - x[1]);
        
        // J[0][1] = -exp(x1 - x2)
        J(0, 1) = -::interval::exp(x[0] - x[1]);
        
        // J[1][0] = 2x1
        J(1, 0) = midPointInterval(2.0) * x[0];
        
        // J[1][1] = 1
        J(1, 1) = midPointInterval(1.0);
        
        return J;
    };
    
    // 初始区间 X0 = ([0, 5], [0, 5])^T
    ik::IntervalVector<2> xn;
    xn[0] = KaucherInterval(0.0, 5.0);
    xn[1] = KaucherInterval(0.0, 5.0);
    
    std::cout << "\nInitial interval X0:\n";
    std::cout << "  x1 = " << xn[0] << std::endl;
    std::cout << "  x2 = " << xn[1] << std::endl;
    
    int maxIterations = 50;
    double tolerance = 1e-2;
    
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        std::cout << "\nIteration " << iter + 1 << ":\n";
        std::cout << "  Current interval xn = " << xn << std::endl;
        
        // 1. 计算中点
        ik::IntervalVector<2> c_vec;
        c_vec[0] = midPointInterval(xn[0].middle());
        c_vec[1] = midPointInterval(xn[1].middle());
        std::cout << "  Midpoint c = " << c_vec << std::endl;
        
        // 2. 在中点处计算 f(c)
        ik::IntervalVector<2> f_c = f(c_vec);
        std::cout << "  f(c) = " << f_c << std::endl;
        
        // 3. 计算残差 r = f(c) - Y_dual（使用对偶区间Y_dual）
        ik::IntervalVector<2> r;
        r[0] = f_c[0] - y1_dual;
        r[1] = f_c[1] - y2_dual;
        std::cout << "  Residual r = " << r << std::endl;
        
        // 4. 计算区间雅可比矩阵 J(X)
        ik::IntervalMatrix<2, 2> J = jacobian(xn);
        std::cout << "  Jacobian matrix J:\n";
        std::cout << "    J(0,0) = " << J(0,0) << "\n    J(0,1) = " << J(0,1) << "\n    J(1,0) = " << J(1,0) << "\n    J(1,1) = " << J(1,1) << std::endl;
        
        // 5. 计算 J 的中点矩阵
        Eigen::MatrixXd J_mid = J.midpoint();
        std::cout << "  Jacobian midpoint matrix J_mid = \n" << J_mid << std::endl;
        
        // 6. 计算 J_mid 的逆矩阵 Y
        Eigen::MatrixXd Y = J_mid.inverse();
        std::cout << "  Inverse matrix Y = \n" << Y << std::endl;
        
        // 7. 计算 term1: c - Y * r (使用中点残差)
        Eigen::VectorXd c_eigen(2);
        c_eigen << c_vec[0].middle(), c_vec[1].middle();
        
        Eigen::VectorXd r_eigen(2);
        r_eigen << r[0].middle(), r[1].middle();
        
        Eigen::VectorXd term1_eigen = c_eigen - Y * r_eigen;
        ik::IntervalVector<2> term1;
        term1[0] = midPointInterval(term1_eigen(0));
        term1[1] = midPointInterval(term1_eigen(1));
        std::cout << "  term1 = " << term1 << std::endl;
        
        // 8. 计算 I - YJ
        // 首先创建单位区间矩阵
        ik::IntervalMatrix<2, 2> I;
        I(0, 0) = midPointInterval(1.0);
        I(0, 1) = midPointInterval(0.0);
        I(1, 0) = midPointInterval(0.0);
        I(1, 1) = midPointInterval(1.0);
        
        // 将 Y 转换为区间矩阵（中点矩阵）
        ik::IntervalMatrix<2, 2> Y_interval;
        for (size_t i = 0; i < 2; ++i)
        {
            for (size_t j = 0; j < 2; ++j)
            {
                Y_interval(i, j) = midPointInterval(Y(i, j));
            }
        }
        
        // 计算 Y * J（区间矩阵乘法）
        ik::IntervalMatrix<2, 2> YJ = Y_interval * J;
        std::cout << "  YJ = " << std::endl;
        std::cout << "    (0,0) = " << YJ(0,0) << "\n    (0,1) = " << YJ(0,1) << "\n    (1,0) = " << YJ(1,0) << "\n    (1,1) = " << YJ(1,1) << std::endl;
        
        // 计算 I - YJ
        ik::IntervalMatrix<2, 2> I_minus_YJ = I - YJ;
        std::cout << "  I - YJ = " << std::endl;
        std::cout << "    (0,0) = " << I_minus_YJ(0,0) << "\n    (0,1) = " << I_minus_YJ(0,1) << "\n    (1,0) = " << I_minus_YJ(1,0) << "\n    (1,1) = " << I_minus_YJ(1,1) << std::endl;
        
        // 9. 计算 box - c
        ik::IntervalVector<2> box_minus_c;
        double delta1 = (xn[0].upper() - xn[0].lower()) / 2.0;
        double delta2 = (xn[1].upper() - xn[1].lower()) / 2.0;
        box_minus_c[0] = KaucherInterval(-delta1, delta1);
        box_minus_c[1] = KaucherInterval(-delta2, delta2);
        std::cout << "  box - c = " << box_minus_c << std::endl;
        
        // 10. 计算 term3 = (I - YJ) * (box - c)（区间矩阵-向量乘法）
        ik::IntervalVector<2> term3 = I_minus_YJ * box_minus_c;
        std::cout << "  term3 = " << term3 << std::endl;
        
        // 11. 计算 K = term1 + term3
        ik::IntervalVector<2> K = term1 + term3;
        std::cout << "  Krawczyk interval K = " << K << std::endl;
        
        // 12. 计算 X_new = X ∩ K
        ik::IntervalVector<2> x_new = xn.intersection(K);
        std::cout << "  New interval x_new = " << x_new << std::endl;
        
        // 13. 确保新区间不为空
        if (x_new[0].isEmpty() || x_new[1].isEmpty())
        {
            std::cout << "  New interval is empty, using Krawczyk interval as new interval\n";
            x_new = K;
        }
        
        // 14. 计算区间宽度
        double width1 = x_new[0].width();
        double width2 = x_new[1].width();
        double max_width = std::max(width1, width2);
        std::cout << "  Interval widths: x1 = " << width1 << ", x2 = " << width2 << ", max = " << max_width << std::endl;
        
        // 15. 更新区间
        xn = x_new;
        
        // 16. 检查收敛
        if (max_width < tolerance)
        {
            std::cout << "\n✓ Converged to " << xn << ", Iterations: " << iter + 1 << std::endl;
            
            // 验证结果
            std::cout << "\n=== Result Verification ===\n";
            // 在收敛区间内取中点进行验证
            double x1_test = xn[0].middle();
            double x2_test = xn[1].middle();
            
            // 计算函数值
            double f1_val = x1_test + std::exp(x1_test - x2_test);
            double f2_val = x1_test * x1_test + x2_test;
            
            std::cout << "  At midpoint x = (" << x1_test << ", " << x2_test << "):\n";
            std::cout << "    f1(x) = " << f1_val << ", should be in " << y1_normal << " range" << std::endl;
            std::cout << "    f2(x) = " << f2_val << ", should be in " << y2_normal << " range" << std::endl;
            
            // 检查是否在目标范围内
            bool f1_in_range = (f1_val >= y1_normal.lower()) && (f1_val <= y1_normal.upper());
            bool f2_in_range = (f2_val >= y2_normal.lower()) && (f2_val <= y2_normal.upper());
            
            if (f1_in_range && f2_in_range)
            {
                std::cout << "  ✓ Midpoint function values are within target range\n";
            }
            else
            {
                std::cout << "  ✗ Midpoint function values are not within target range\n";
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
    toleranceEmbeddingCorrect();
    return 0;
}
