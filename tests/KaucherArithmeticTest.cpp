#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
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

// Test 1: Algebraic difference test X + [2, 3] = [5, 6]
TEST(KaucherArithmeticTest, AlgebraicDifference) {
    // Define equation f(X) = X + [2, 3] - [5, 6] = 0
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1> {
        ik::IntervalVector<1> result;
        result[0] = x[0] + KaucherInterval(2.0, 3.0) - KaucherInterval(5.0, 6.0);
        return result;
    };
    
    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1> {
        ik::IntervalMatrix<1, 1> J;
        J(0, 0) = midPointInterval(1.0);
        return J;
    };
    
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(0.0, 10.0);
    
    int maxIterations = 10;
    double tolerance = 1e-8;
    bool converged = false;
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        // 1. Calculate midpoint
        double c = xn[0].middle();
        
        // 2. Calculate f at midpoint
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec);
        
        // 3. Calculate Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian(xn);
        Eigen::MatrixXd J_mid = J.midpoint();
        Eigen::MatrixXd Y = J_mid.inverse();
        
        // 4. Calculate term1: c - Y * f(c)
        double term1_val = c - Y(0, 0) * f_c[0].middle();
        ik::IntervalVector<1> term1;
        term1[0] = midPointInterval(term1_val);
        
        // 5. Calculate I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(1, 1);
        Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
        ik::IntervalMatrix<1, 1> I_minus_YJ_interval;
        I_minus_YJ_interval(0, 0) = midPointInterval(I_minus_YJ(0, 0));
        
        // 6. Calculate box - c
        double delta = (xn[0].upper() - xn[0].lower()) / 2.0;
        ik::IntervalVector<1> box_minus_c;
        box_minus_c[0] = KaucherInterval(-delta, delta);
        
        // 7. Calculate term3
        ik::IntervalVector<1> term3 = I_minus_YJ_interval * box_minus_c;
        
        // 8. Calculate K = term1 + term3
        ik::IntervalVector<1> K = term1 + term3;
        
        // 9. Calculate X_new = X ∩ K
        ik::IntervalVector<1> x_new = xn.intersection(K);
        
        // Check convergence
        if (x_new[0].width() < tolerance) {
            converged = true;
            break;
        }
        
        xn = x_new;
        
        // If interval is empty, try using improper interval
        if (x_new[0].isEmpty()) {
            xn = K;
        }
    }
    
    ASSERT_TRUE(converged) << "Algebraic difference test did not converge";
}

// Test 2: Reverse width iteration f(x) = x - tanh(2x) = 0
TEST(KaucherArithmeticTest, ReverseWidthIteration) {
    auto f = [](const ik::IntervalVector<1>& x) -> ik::IntervalVector<1> {
        ik::IntervalVector<1> result;
        result[0] = x[0] - ::interval::tanh(x[0] * 2.0);
        return result;
    };
    
    auto jacobian = [](const ik::IntervalVector<1>& x) -> ik::IntervalMatrix<1, 1> {
        ik::IntervalMatrix<1, 1> J;
        KaucherInterval tanh_2x = ::interval::tanh(x[0] * 2.0);
        J(0, 0) = midPointInterval(1.0) - midPointInterval(2.0) * (midPointInterval(1.0) - tanh_2x * tanh_2x);
        return J;
    };
    
    // Start with very large initial interval
    ik::IntervalVector<1> xn;
    xn[0] = KaucherInterval(-10.0, 10.0);
    
    int maxIterations = 20;
    double tolerance = 1e-8;
    bool converged = false;
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        // 1. Calculate midpoint
        double c = xn[0].middle();
        
        // 2. Calculate f at midpoint
        ik::IntervalVector<1> c_vec;
        c_vec[0] = midPointInterval(c);
        ik::IntervalVector<1> f_c = f(c_vec);
        
        // 3. Calculate Jacobian
        ik::IntervalMatrix<1, 1> J = jacobian(xn);
        Eigen::MatrixXd J_mid = J.midpoint();
        
        // Check if Jacobian is invertible
        double det = J_mid(0, 0);
        if (std::abs(det) < 1e-12) {
            det = (det > 0) ? 1e-12 : -1e-12;
        }
        Eigen::MatrixXd Y(1, 1);
        Y(0, 0) = 1.0 / det;
        
        // 4. Calculate term1: c - Y * f(c)
        double term1_val = c - Y(0, 0) * f_c[0].middle();
        ik::IntervalVector<1> term1;
        term1[0] = midPointInterval(term1_val);
        
        // 5. Calculate I - YJ
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(1, 1);
        Eigen::MatrixXd I_minus_YJ = I - Y * J.midpoint();
        ik::IntervalMatrix<1, 1> I_minus_YJ_interval;
        I_minus_YJ_interval(0, 0) = midPointInterval(I_minus_YJ(0, 0));
        
        // 6. Calculate box - c
        double delta = (xn[0].upper() - xn[0].lower()) / 2.0;
        ik::IntervalVector<1> box_minus_c;
        box_minus_c[0] = KaucherInterval(-delta, delta);
        
        // 7. Calculate term3
        ik::IntervalVector<1> term3 = I_minus_YJ_interval * box_minus_c;
        
        // 8. Calculate K = term1 + term3
        ik::IntervalVector<1> K = term1 + term3;
        
        // 9. Calculate X_new = X ∩ K
        ik::IntervalVector<1> x_new = xn.intersection(K);
        
        // Check convergence
        if (x_new[0].width() < tolerance) {
            converged = true;
            break;
        }
        
        xn = x_new;
        
        // If interval is empty or width increases, try using improper interval
        if (x_new[0].isEmpty() || x_new[0].width() > xn[0].width()) {
            xn = K;
        }
    }
    
    ASSERT_TRUE(converged) << "Reverse width iteration did not converge";
}

// Test 3: Basic Kaucher interval properties
TEST(KaucherArithmeticTest, BasicProperties) {
    KaucherInterval proper(1.0, 2.0);
    KaucherInterval improper(2.0, 1.0);
    KaucherInterval a(2.0, 3.0);
    KaucherInterval b(5.0, 6.0);
    
    // Test proper interval properties
    EXPECT_LE(proper.lower(), proper.upper());
    EXPECT_FALSE(proper.isImproper());
    
    // Test improper interval properties
    EXPECT_GT(improper.lower(), improper.upper());
    EXPECT_TRUE(improper.isImproper());
    
    // Test dual operation
    KaucherInterval dual = improper.dual();
    EXPECT_LE(dual.lower(), dual.upper());
    EXPECT_FALSE(dual.isImproper());
    
    // Test algebraic difference
    KaucherInterval algebraic_diff = a.algebraicDiff(b);
    KaucherInterval classic_diff = a - b;
    KaucherInterval inverse_sum = a + (-b);
    
    // These should all be different for Kaucher intervals
    EXPECT_NE(algebraic_diff.lower(), classic_diff.lower());
    EXPECT_NE(algebraic_diff.upper(), classic_diff.upper());
}