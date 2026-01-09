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
TEST(ToleranceEmbeddingCorrectTest, GeneralizedKrawczykIteration) {
    // Known target output range Y = ([2.9, 3.1], [4.9, 5.1])^T
    // To solve for inner enclosure, need to convert Y to dual form (Improper Interval)
    KaucherInterval y1_normal(2.9, 3.1);
    KaucherInterval y2_normal(4.9, 5.1);
    
    // Convert to dual form, triggering inner enclosure solving mode
    KaucherInterval y1_dual = dualInterval(y1_normal);
    KaucherInterval y2_dual = dualInterval(y2_normal);
    
    // 定义系统方程 f(x) = Y
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2> {
        ik::IntervalVector<2> result;
        result[0] = x[0] + ::interval::exp(x[0] - x[1]);
        result[1] = x[0] * x[0] + x[1];
        return result;
    };
    
    // Calculate interval Jacobian matrix J(X)
    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2> {
        ik::IntervalMatrix<2, 2> J;
        
        // Calculate partial derivatives
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
    
    // Initial interval X0 = ([0, 5], [0, 5])^T
    ik::IntervalVector<2> xn;
    xn[0] = KaucherInterval(0.0, 5.0);
    xn[1] = KaucherInterval(0.0, 5.0);
    
    int maxIterations = 50;
    double tolerance = 1e-2;
    bool converged = false;
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        // 1. Calculate midpoint
        ik::IntervalVector<2> c_vec;
        c_vec[0] = midPointInterval(xn[0].middle());
        c_vec[1] = midPointInterval(xn[1].middle());
        
        // 2. Calculate f(c) at midpoint
        ik::IntervalVector<2> f_c = f(c_vec);
        
        // 3. Calculate residual r = f(c) - Y_dual (using dual interval Y_dual)
        ik::IntervalVector<2> r;
        r[0] = f_c[0] - y1_dual;
        r[1] = f_c[1] - y2_dual;
        
        // 4. Calculate interval Jacobian matrix J(X)
        ik::IntervalMatrix<2, 2> J = jacobian(xn);
        
        // 5. Calculate midpoint matrix of J
        Eigen::MatrixXd J_mid = J.midpoint();
        
        // 6. Calculate inverse matrix Y of J_mid
        Eigen::MatrixXd Y = J_mid.inverse();
        
        // 7. Calculate term1: c - Y * r (using midpoint residual)
        Eigen::VectorXd c_eigen(2);
        c_eigen << c_vec[0].middle(), c_vec[1].middle();
        
        Eigen::VectorXd r_eigen(2);
        r_eigen << r[0].middle(), r[1].middle();
        
        Eigen::VectorXd term1_eigen = c_eigen - Y * r_eigen;
        ik::IntervalVector<2> term1;
        term1[0] = midPointInterval(term1_eigen(0));
        term1[1] = midPointInterval(term1_eigen(1));
        
        // 8. Calculate I - YJ
        // First create identity interval matrix
        ik::IntervalMatrix<2, 2> I;
        I(0, 0) = midPointInterval(1.0);
        I(0, 1) = midPointInterval(0.0);
        I(1, 0) = midPointInterval(0.0);
        I(1, 1) = midPointInterval(1.0);
        
        // Convert Y to interval matrix (midpoint matrix)
        ik::IntervalMatrix<2, 2> Y_interval;
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                Y_interval(i, j) = midPointInterval(Y(i, j));
            }
        }
        
        // Calculate Y * J (interval matrix multiplication)
        ik::IntervalMatrix<2, 2> YJ = Y_interval * J;
        
        // Calculate I - YJ
        ik::IntervalMatrix<2, 2> I_minus_YJ = I - YJ;
        
        // 9. Calculate box - c
        ik::IntervalVector<2> box_minus_c;
        double delta1 = (xn[0].upper() - xn[0].lower()) / 2.0;
        double delta2 = (xn[1].upper() - xn[1].lower()) / 2.0;
        box_minus_c[0] = KaucherInterval(-delta1, delta1);
        box_minus_c[1] = KaucherInterval(-delta2, delta2);
        
        // 10. Calculate term3 = (I - YJ) * (box - c)（区间矩阵-向量乘法）
        ik::IntervalVector<2> term3 = I_minus_YJ * box_minus_c;
        
        // 11. Calculate K = term1 + term3
        ik::IntervalVector<2> K = term1 + term3;
        
        // 12. Calculate X_new = X ∩ K
        ik::IntervalVector<2> x_new = xn.intersection(K);
        
        // 13. Ensure new interval is not empty
        if (x_new[0].isEmpty() || x_new[1].isEmpty()) {
            x_new = K;
        }
        
        // 14. Calculate interval widths
        double width1 = x_new[0].width();
        double width2 = x_new[1].width();
        double max_width = std::max(width1, width2);
        
        // 15. Update interval
        xn = x_new;
        
        // 16. Check convergence
        if (max_width < tolerance) {
            converged = true;
            
            // Verify result
            // Test midpoint in converged interval
            double x1_test = xn[0].middle();
            double x2_test = xn[1].middle();
            
            // Calculate function values
            double f1_val = x1_test + std::exp(x1_test - x2_test);
            double f2_val = x1_test * x1_test + x2_test;
            
            // Check if within target range
            bool f1_in_range = (f1_val >= y1_normal.lower()) && (f1_val <= y1_normal.upper());
            bool f2_in_range = (f2_val >= y2_normal.lower()) && (f2_val <= y2_normal.upper());
            
            // Both values should be within target range
            EXPECT_TRUE(f1_in_range) << "f1(x) = " << f1_val << " is not in range " << y1_normal;
            EXPECT_TRUE(f2_in_range) << "f2(x) = " << f2_val << " is not in range " << y2_normal;
            
            break;
        }
    }
    
    // Test should converge within max iterations
    EXPECT_TRUE(converged) << "Tolerance Embedding Correct test did not converge within maximum iterations";
}
