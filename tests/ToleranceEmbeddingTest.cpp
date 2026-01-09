#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <gtest/gtest.h>
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

// Define system equation f(x) = Y
// f(x) = [x1 + exp(x1 - x2), x1² + x2]^T
TEST(ToleranceEmbeddingTest, GeneralizedKrawczykIteration) {
    // Known target output range Y = ([2.9, 3.1], [4.9, 5.1])^T
    // To solve for inner enclosure, need to convert Y to dual form (Improper Interval)
    KaucherInterval y1_normal(2.9, 3.1);
    KaucherInterval y2_normal(4.9, 5.1);
    
    // Convert to dual form, triggering inner enclosure solving mode
    KaucherInterval y1_dual = dualInterval(y1_normal);
    KaucherInterval y2_dual = dualInterval(y2_normal);
    
    // Initial interval X0 = ([0, 5], [0, 5])^T
    ik::IntervalVector<2> xn;
    xn[0] = KaucherInterval(0.0, 5.0);
    xn[1] = KaucherInterval(0.0, 5.0);
    
    int maxIterations = 20;
    double tolerance = 1e-4;
    bool converged = false;
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        // 1. Calculate midpoint
        double c1 = xn[0].middle();
        double c2 = xn[1].middle();
        
        // 2. Calculate function value and Jacobian at midpoint
        double dx1 = 0.01;
        double dx2 = 0.01;
        
        // Calculate function value f(c)
        double f1_c = c1 + std::exp(c1 - c2);
        double f2_c = c1 * c1 + c2;
        
        // Calculate approximate Jacobian using finite differences
        double f1_c_dx1 = (c1 + dx1) + std::exp((c1 + dx1) - c2);
        double f1_c_dx2 = c1 + std::exp(c1 - (c2 + dx2));
        double f2_c_dx1 = (c1 + dx1) * (c1 + dx1) + c2;
        double f2_c_dx2 = c1 * c1 + (c2 + dx2);
        
        double J11 = (f1_c_dx1 - f1_c) / dx1;
        double J12 = (f1_c_dx2 - f1_c) / dx2;
        double J21 = (f2_c_dx1 - f2_c) / dx1;
        double J22 = (f2_c_dx2 - f2_c) / dx2;
        
        // 3. Construct midpoint matrix of Jacobian
        Eigen::MatrixXd J_mid(2, 2);
        J_mid << J11, J12, J21, J22;
        
        // 4. Calculate inverse matrix Y
        Eigen::MatrixXd Y = J_mid.inverse();
        
        // 5. Calculate term1: c - Y * (f(c) - Y_dual)
        // Note: Using dual interval Y_dual here
        Eigen::VectorXd c_eigen(2);
        c_eigen << c1, c2;
        
        Eigen::VectorXd f_c_eigen(2);
        f_c_eigen << f1_c, f2_c;
        
        Eigen::VectorXd Y_dual_eigen(2);
        Y_dual_eigen << y1_dual.middle(), y2_dual.middle();
        
        Eigen::VectorXd term1_eigen = c_eigen - Y * (f_c_eigen - Y_dual_eigen);
        
        // 6. Calculate I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);
        Eigen::MatrixXd I_minus_YJ = I - Y * J_mid;
        
        // 7. Calculate box - c
        double delta1 = (xn[0].upper() - xn[0].lower()) / 2.0;
        double delta2 = (xn[1].upper() - xn[1].lower()) / 2.0;
        
        // 8. Calculate term3 = (I - YJ) * (box - c)
        // Directly calculate in interval form here
        double term3_1_lower = I_minus_YJ(0, 0) * (-delta1) + I_minus_YJ(0, 1) * (-delta2);
        double term3_1_upper = I_minus_YJ(0, 0) * delta1 + I_minus_YJ(0, 1) * delta2;
        double term3_2_lower = I_minus_YJ(1, 0) * (-delta1) + I_minus_YJ(1, 1) * (-delta2);
        double term3_2_upper = I_minus_YJ(1, 0) * delta1 + I_minus_YJ(1, 1) * delta2;
        
        // Ensure interval is valid (lower <= upper)
        if (term3_1_lower > term3_1_upper) std::swap(term3_1_lower, term3_1_upper);
        if (term3_2_lower > term3_2_upper) std::swap(term3_2_lower, term3_2_upper);
        
        // 9. Calculate Krawczyk interval K = term1 + term3
        KaucherInterval K1(term1_eigen(0) + term3_1_lower, term1_eigen(0) + term3_1_upper);
        KaucherInterval K2(term1_eigen(1) + term3_2_lower, term1_eigen(1) + term3_2_upper);
        
        ik::IntervalVector<2> K;
        K[0] = K1;
        K[1] = K2;
        
        // 10. Calculate X_new = X ∩ K
        ik::IntervalVector<2> x_new;
        x_new[0] = xn[0].meet(K[0]);
        x_new[1] = xn[1].meet(K[1]);
        
        // 11. Ensure new interval is not empty
        if (x_new[0].isEmpty() || x_new[1].isEmpty()) {
            x_new = K;
        }
        
        // 12. 计算区间宽度
        double width1 = x_new[0].width();
        double width2 = x_new[1].width();
        double max_width = std::max(width1, width2);
        
        // 13. 更新区间
        xn = x_new;
        
        // 14. 检查收敛
        if (max_width < tolerance) {
            converged = true;
            
            // 验证结果
            // 在收敛区间内取几个点进行验证
            double x1_test1 = xn[0].lower();
            double x2_test1 = xn[1].lower();
            double x1_test2 = xn[0].upper();
            double x2_test2 = xn[1].upper();
            double x1_test3 = xn[0].middle();
            double x2_test3 = xn[1].middle();
            
            // 计算函数值
            auto compute_f = [](double x1, double x2) -> std::pair<double, double> {
                double f1 = x1 + std::exp(x1 - x2);
                double f2 = x1 * x1 + x2;
                return {f1, f2};
            };
            
            auto [f1_val1, f2_val1] = compute_f(x1_test1, x2_test1);
            auto [f1_val2, f2_val2] = compute_f(x1_test2, x2_test2);
            auto [f1_val3, f2_val3] = compute_f(x1_test3, x2_test3);
            
            // 检查是否在目标范围内
            bool f1_in_range1 = (f1_val1 >= y1_normal.lower()) && (f1_val1 <= y1_normal.upper());
            bool f2_in_range1 = (f2_val1 >= y2_normal.lower()) && (f2_val1 <= y2_normal.upper());
            bool f1_in_range2 = (f1_val2 >= y1_normal.lower()) && (f1_val2 <= y1_normal.upper());
            bool f2_in_range2 = (f2_val2 >= y2_normal.lower()) && (f2_val2 <= y2_normal.upper());
            bool f1_in_range3 = (f1_val3 >= y1_normal.lower()) && (f1_val3 <= y1_normal.upper());
            bool f2_in_range3 = (f2_val3 >= y2_normal.lower()) && (f2_val3 <= y2_normal.upper());
            
            bool all_valid = f1_in_range1 && f2_in_range1 && f1_in_range2 && f2_in_range2 && f1_in_range3 && f2_in_range3;
            
            // All test points should be within target range
            EXPECT_TRUE(all_valid) << "Some test points' function values are not within target range";
            break;
        }
    }
    
    // Test should converge within max iterations
    EXPECT_TRUE(converged) << "Tolerance Embedding test did not converge within maximum iterations";
}
