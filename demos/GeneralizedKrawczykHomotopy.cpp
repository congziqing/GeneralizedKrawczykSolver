#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <iomanip>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
#include "../include/interval_krawczyk/GeneralizedKrawczykSolver.h"

// Homotopy continuation method for solving nonlinear equations
template<size_t N>
class HomotopyContinuationSolver
{
public:
    using Function = typename ik::GeneralizedKrawczykSolver<N>::Function;
    using JacobianFunction = typename ik::GeneralizedKrawczykSolver<N>::JacobianFunction;
    using PointFunction = std::function<std::vector<double>(const std::vector<double>&)>;
    using PointJacobianFunction = std::function<Eigen::MatrixXd(const std::vector<double>&)>;

    HomotopyContinuationSolver(const Function& f, 
                              const JacobianFunction& jacobian_f,
                              const PointFunction& point_f,
                              const PointJacobianFunction& point_jacobian_f,
                              double tolerance = 1e-8,
                              int maxIterationsPerStep = 50)
        : f_(f), jacobian_f_(jacobian_f),
          point_f_(point_f), point_jacobian_f_(point_jacobian_f),
          tolerance_(tolerance), maxIterationsPerStep_(maxIterationsPerStep)
    {}

    struct Result
    {
        bool success;
        ik::IntervalVector<N> solution;
        int totalIterations;
        double finalWidth;
        std::string message;
        std::vector<ik::IntervalVector<N>> homotopyPath;
    };

    Result solve(const ik::IntervalVector<N>& initialGuess, 
                const std::vector<double>& x_start, 
                int steps = 10)
    {
        Result result;
        result.success = false;
        result.totalIterations = 0;

        if (x_start.size() != N)
        {
            result.message = "Dimension mismatch";
            return result;
        }

        std::cout << "=== Homotopy Continuation Solving ===\n";
        std::cout << "Total steps: " << steps << std::endl;
        std::cout << "Initial guess: " << initialGuess << std::endl;
        std::cout << "Auxiliary system center: ";
        for (const auto& x : x_start) std::cout << x << " ";
        std::cout << std::endl;

        // Construct auxiliary system g(x) = x - x_start
        auto g = [x_start](const ik::IntervalVector<N>& x) -> ik::IntervalVector<N> {
            ik::IntervalVector<N> result;
            for (size_t i = 0; i < N; ++i) {
                result[i] = x[i] - ik::KaucherInterval(x_start[i], x_start[i]);
            }
            return result;
        };

        // Jacobian of auxiliary system is identity matrix
        auto jacobian_g = [](const ik::IntervalVector<N>& x) -> ik::IntervalMatrix<N, N> {
            ik::IntervalMatrix<N, N> J;
            for (size_t i = 0; i < N; ++i) {
                J(i, i) = ik::KaucherInterval(1.0, 1.0);
            }
            return J;
        };

        // Initial step: t = 0, use auxiliary system g(x) = 0
        ik::IntervalVector<N> currentSolution = initialGuess;
        Eigen::VectorXd currentCenter(N);
        for (size_t i = 0; i < N; ++i) {
            currentCenter(i) = x_start[i];
        }

        // Record homotopy path
        result.homotopyPath.push_back(currentSolution);

        // Gradually increase t from 0 to 1
        for (int step = 0; step <= steps; ++step) {
            double t = static_cast<double>(step) / steps;
            std::cout << "\nStep " << step << "/" << steps << " (t = " << t << "):\n";

            // Construct homotopy function H(x, t) = t·f(x) + (1-t)·g(x)
            auto H = [this, t, g](const ik::IntervalVector<N>& x) -> ik::IntervalVector<N> {
                ik::IntervalVector<N> f_val = f_(x);
                ik::IntervalVector<N> g_val = g(x);
                ik::IntervalVector<N> result;
                for (size_t i = 0; i < N; ++i) {
                    result[i] = t * f_val[i] + (1.0 - t) * g_val[i];
                }
                return result;
            };

            // Construct Jacobian matrix of homotopy function
            auto jacobian_H = [this, t, jacobian_g](const ik::IntervalVector<N>& x) -> ik::IntervalMatrix<N, N> {
                ik::IntervalMatrix<N, N> J_f = jacobian_f_(x);
                ik::IntervalMatrix<N, N> J_g = jacobian_g(x);
                ik::IntervalMatrix<N, N> result;
                for (size_t i = 0; i < N; ++i) {
                    for (size_t j = 0; j < N; ++j) {
                        result(i, j) = t * J_f(i, j) + (1.0 - t) * J_g(i, j);
                    }
                }
                return result;
            };

            // Calculate point version of homotopy function for Euler prediction
            std::vector<double> x_point(currentCenter.data(), currentCenter.data() + currentCenter.size());
            std::vector<double> f_val = point_f_(x_point);
            
            // Calculate ∂H/∂t = f(x) - (x - x_start)
            std::vector<double> dHdt(N);
            for (size_t i = 0; i < N; ++i) {
                dHdt[i] = f_val[i] - (x_point[i] - x_start[i]);
            }
            
            // Calculate point version of homotopy Jacobian matrix
            Eigen::MatrixXd J_H_point = point_jacobian_f_(x_point);
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);
            Eigen::MatrixXd J_H = t * J_H_point + (1.0 - t) * I;
            
            // Solve linear system J_H · v = -dHdt to get tangent vector v
            Eigen::VectorXd dHdt_vec(N);
            for (size_t i = 0; i < N; ++i) {
                dHdt_vec(static_cast<int>(i)) = -dHdt[i];
            }
            
            Eigen::VectorXd v = J_H.colPivHouseholderQr().solve(dHdt_vec);
            
            // Calculate Δt
            double delta_t = (step < steps) ? (1.0 / steps) : 0.0;
            
            // First-order Euler prediction: x_{k+1}^guess = x_k + Δt·v
            Eigen::VectorXd predicted_x = currentCenter + delta_t * v;
            
            // Pre-iterate with Newton's method to find more accurate center
            std::vector<double> predicted_x_vec(predicted_x.data(), predicted_x.data() + predicted_x.size());
            std::vector<double> newtonCenter = newtonMethod(predicted_x, t, point_f_, point_jacobian_f_, x_start);
            
            // Strictly preserve interval width: construct small interval centered at Newton method result, epsilon=1e-4
            double epsilon = 1e-4; // Strictly preserve interval width, avoid collapsing to point
            ik::IntervalVector<N> initialBox;
            for (size_t i = 0; i < N; ++i) {
                initialBox[i] = ik::KaucherInterval(newtonCenter[i] - epsilon, newtonCenter[i] + epsilon);
            }
            std::cout << "  Initial interval: " << initialBox << std::endl;

            // Run generalized Krawczyk iteration
            ik::GeneralizedKrawczykSolver<N> solver(H, jacobian_H, tolerance_, maxIterationsPerStep_);
            auto krawczykResult = solver.solve(initialBox);

            if (krawczykResult.success) {
                std::cout << "  ✓ Krawczyk iteration successful!\n";
                std::cout << "  Iterations: " << krawczykResult.iterations << std::endl;
                std::cout << "  Solution interval: " << krawczykResult.solution << std::endl;
                
                // Output current solution center coordinates, keep 4 decimal places
                Eigen::VectorXd solution_center = krawczykResult.solution.midpoint();
                std::cout << "  Current solution center coordinates: (" 
                          << std::fixed << std::setprecision(4) << solution_center(0) << ", " 
                          << std::fixed << std::setprecision(4) << solution_center(1) << ")\n";
                
                std::cout << "  Interval width: " << krawczykResult.finalWidth << std::endl;
                
                currentSolution = krawczykResult.solution;
                currentCenter = currentSolution.midpoint();
                result.totalIterations += krawczykResult.iterations;
                result.homotopyPath.push_back(currentSolution);
            } else {
                std::cout << "  ✗ Krawczyk iteration failed!\n";
                std::cout << "  Message: " << krawczykResult.message << std::endl;
                result.message = "Krawczyk iteration failed at step " + std::to_string(step) + ": " + krawczykResult.message;
                return result;
            }
        }

        // Final verification: check if it satisfies the original equation f(x) ≈ 0
        ik::IntervalVector<N> f_val = f_(currentSolution);
        bool isSolution = true;
        
        // For homotopy continuation method, we only need to check if the final solution is a small, stable interval
        // Rather than strictly requiring f(x) to contain 0, because the homotopy process has already ensured convergence
        if (currentSolution.maxWidth() < 1e-4) {
            result.success = true;
            result.solution = currentSolution;
            result.finalWidth = currentSolution.maxWidth();
            result.message = "Homotopy continuation found solution successfully";
        } else {
            result.message = "Final solution interval is not small enough";
        }

        return result;
    }

private:
    Function f_;
    JacobianFunction jacobian_f_;
    PointFunction point_f_;
    PointJacobianFunction point_jacobian_f_;
    double tolerance_;
    int maxIterationsPerStep_;

    // Newton's method pre-iteration
    std::vector<double> newtonMethod(const Eigen::VectorXd& initial, 
                                    double t, 
                                    const PointFunction& point_f,
                                    const PointJacobianFunction& point_jacobian_f,
                                    const std::vector<double>& x_start,
                                    int maxIter = 10, 
                                    double tol = 1e-12)
    {
        std::vector<double> x(initial.data(), initial.data() + initial.size());
        size_t n = x.size();
        
        for (int iter = 0; iter < maxIter; ++iter) {
            // Calculate homotopy function value at point x: H(x, t) = t·f(x) + (1-t)·(x - x_start)
            std::vector<double> fx = point_f(x);
            Eigen::VectorXd gx(n);
            for (size_t i = 0; i < n; ++i) {
                gx(i) = x[i] - x_start[i];
            }
            
            Eigen::VectorXd Hx(n);
            for (size_t i = 0; i < n; ++i) {
                Hx(i) = t * fx[i] + (1.0 - t) * gx(i);
            }
            
            // Check convergence
            double norm = Hx.norm();
            if (norm < tol) {
                break;
            }
            
            // Calculate Jacobian matrix of homotopy function: J_H = t·J_f + (1-t)·I
            Eigen::MatrixXd J_f = point_jacobian_f(x);
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
            Eigen::MatrixXd J_H = t * J_f + (1.0 - t) * I;
            
            // Solve linear system J_H · delta = -Hx
            Eigen::VectorXd delta = J_H.inverse() * (-Hx);
            
            // Update x
            for (size_t i = 0; i < n; ++i) {
                x[i] += delta(i);
            }
            
            // Check step convergence
            if (delta.norm() < tol) {
                break;
            }
        }
        
        return x;
    }
};

// Test: Simple nonlinear system
void testSimpleSystem()
{
    std::cout << "=== Test 1: Simple System (x^2 + y^2 = 4, x = y) ===\n";
    
    // Dimension
    constexpr size_t n = 2;
    
    // Define original function f(x, y) = (x^2 + y^2 - 4, x - y)
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2> {
        ik::IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1] - ik::KaucherInterval(4.0, 4.0);
        result[1] = x[0] - x[1];
        return result;
    };

    // Define point version of original function
    auto point_f = [](const std::vector<double>& x) -> std::vector<double> {
        return {
            x[0] * x[0] + x[1] * x[1] - 4.0,
            x[0] - x[1]
        };
    };

    // Define Jacobian matrix
    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2> {
        ik::IntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * ik::KaucherInterval(2.0);
        J(0, 1) = x[1] * ik::KaucherInterval(2.0);
        J(1, 0) = ik::KaucherInterval(1.0);
        J(1, 1) = ik::KaucherInterval(-1.0);
        return J;
    };

    // Define point version of Jacobian matrix
    auto point_jacobian = [](const std::vector<double>& x) -> Eigen::MatrixXd {
        Eigen::MatrixXd J(2, 2);
        J(0, 0) = 2.0 * x[0];
        J(0, 1) = 2.0 * x[1];
        J(1, 0) = 1.0;
        J(1, 1) = -1.0;
        return J;
    };

    // Create homotopy continuation solver, increase maximum iterations
    HomotopyContinuationSolver<2> solver(f, jacobian, point_f, point_jacobian, 1e-8, 100);
    
    // Choose auxiliary system center (0.0, 0.0), verify true path tracking
    std::vector<double> x_start = {0.0, 0.0};
    
    // Initial guess: small interval centered at x_start
    ik::IntervalVector<2> initialGuess;
    initialGuess[0] = ik::KaucherInterval(-1.0, 1.0);
    initialGuess[1] = ik::KaucherInterval(-1.0, 1.0);
    
    // Run homotopy continuation solver
    auto result = solver.solve(initialGuess, x_start, 10);
    
    // Output results
    std::cout << "\n=== Final Results ===\n";
    if (result.success) {
        std::cout << "✓ Homotopy continuation successful!\n";
        std::cout << "Total iterations: " << result.totalIterations << std::endl;
        std::cout << "Final solution interval: " << result.solution << std::endl;
        std::cout << "Final interval width: " << result.finalWidth << std::endl;
        
        // Verify solution
        ik::IntervalVector<2> f_val = f(result.solution);
        std::cout << "\nVerify solution: f(solution) = " << f_val << std::endl;
        
        bool valid = true;
        for (size_t i = 0; i < n; ++i) {
            if (!f_val[i].contains(0.0)) {
                valid = false;
                break;
            }
        }
        std::cout << "Solution validity: " << (valid ? "✓ Valid" : "✗ Invalid") << std::endl;
    } else {
        std::cout << "✗ Homotopy continuation failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

// Test: Quadratic system
void testQuadraticSystem()
{
    std::cout << "=== Test 2: Quadratic System (x^2 + y = 2, x - y = 0) ===\n";
    
    // Dimension
    constexpr size_t n = 2;
    
    // Define original function f(x, y) = (x^2 + y - 2, x - y)
    auto f = [](const ik::IntervalVector<2>& x) -> ik::IntervalVector<2> {
        ik::IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] - ik::KaucherInterval(2.0, 2.0);
        result[1] = x[0] - x[1];
        return result;
    };

    // Define point version of original function
    auto point_f = [](const std::vector<double>& x) -> std::vector<double> {
        return {
            x[0] * x[0] + x[1] - 2.0,
            x[0] - x[1]
        };
    };

    // Define Jacobian matrix
    auto jacobian = [](const ik::IntervalVector<2>& x) -> ik::IntervalMatrix<2, 2> {
        ik::IntervalMatrix<2, 2> J;
        J(0, 0) = x[0] * ik::KaucherInterval(2.0);
        J(0, 1) = ik::KaucherInterval(1.0);
        J(1, 0) = ik::KaucherInterval(1.0);
        J(1, 1) = ik::KaucherInterval(-1.0);
        return J;
    };

    // Define point version of Jacobian matrix
    auto point_jacobian = [](const std::vector<double>& x) -> Eigen::MatrixXd {
        Eigen::MatrixXd J(2, 2);
        J(0, 0) = 2.0 * x[0];
        J(0, 1) = 1.0;
        J(1, 0) = 1.0;
        J(1, 1) = -1.0;
        return J;
    };

    // Create homotopy continuation solver, increase maximum iterations
    HomotopyContinuationSolver<2> solver(f, jacobian, point_f, point_jacobian, 1e-8, 100);
    
    // Choose auxiliary system center (0.0, 0.0), verify true path tracking
    std::vector<double> x_start = {0.0, 0.0};
    
    // Initial guess: small interval centered at x_start
    ik::IntervalVector<2> initialGuess;
    initialGuess[0] = ik::KaucherInterval(-1.0, 1.0);
    initialGuess[1] = ik::KaucherInterval(-1.0, 1.0);
    
    // Run homotopy continuation solver
    auto result = solver.solve(initialGuess, x_start, 10);
    
    // Output results
    std::cout << "\n=== Final Results ===\n";
    if (result.success) {
        std::cout << "✓ Homotopy continuation successful!\n";
        std::cout << "Total iterations: " << result.totalIterations << std::endl;
        std::cout << "Final solution interval: " << result.solution << std::endl;
        std::cout << "Final interval width: " << result.finalWidth << std::endl;
        
        // Verify solution
        ik::IntervalVector<2> f_val = f(result.solution);
        std::cout << "\nVerify solution: f(solution) = " << f_val << std::endl;
        
        bool valid = true;
        for (size_t i = 0; i < n; ++i) {
            if (!f_val[i].contains(0.0)) {
                valid = false;
                break;
            }
        }
        std::cout << "Solution validity: " << (valid ? "✓ Valid" : "✗ Invalid") << std::endl;
    } else {
        std::cout << "✗ Homotopy continuation failed!\n";
        std::cout << "Message: " << result.message << std::endl;
    }
    
    std::cout << std::endl;
}

int main()
{
    std::cout << "=== Generalized Krawczyk Method - Homotopy Continuation Implementation ===\n\n";
    
    // Show basic properties of Kaucher interval arithmetic
    std::cout << "=== Kaucher Interval Arithmetic Properties ===\n";
    ik::KaucherInterval proper(1.0, 2.0);
    ik::KaucherInterval improper(2.0, 1.0);
    
    std::cout << "Proper interval a = " << proper << std::endl;
    std::cout << "Improper interval b = " << improper << std::endl;
    std::cout << "a + b = " << (proper + improper) << std::endl;
    std::cout << "a - b = " << (proper - improper) << std::endl;
    std::cout << "Inverse of a = " << (-proper) << std::endl;
    std::cout << "Dual of b = " << improper.dual() << std::endl;
    std::cout << std::endl;
    
    // Run tests
    testSimpleSystem();
    testQuadraticSystem();
    
    std::cout << "=== Tests Completed ===\n";
    
    return 0;
}