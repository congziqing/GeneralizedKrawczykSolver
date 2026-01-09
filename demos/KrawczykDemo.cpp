#include <iostream>
#include <iomanip>
#include "PSGMDirectedInterval.h"
#include "IntervalVector.h"
#include "IntervalMatrix.h"
#include "KrawczykSolver.h"

void printResult(const KrawczykResult& result)
{
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "求解结果:" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    if (result.success)
    {
        std::cout << "✓ 成功找到解!" << std::endl;
        std::cout << "迭代次数: " << result.iterations << std::endl;
        std::cout << "最终区间宽度: " << std::scientific << std::setprecision(6)
                  << result.finalWidth << std::endl;
        std::cout << "\n解区间:" << std::endl;

        for (size_t i = 0; i < result.solution.size(); ++i)
        {
            std::cout << "  x[" << i << "] = [" << std::fixed << std::setprecision(10)
                      << result.solution[i].lower() << ", "
                      << result.solution[i].upper() << "]" << std::endl;
            std::cout << "      中点: " << result.solution[i].middle() << std::endl;
        }

        double mid = result.solution[0].middle();
        if (result.solution.size() > 1)
            mid = result.solution[1].middle();

        std::cout << "\n消息: " << result.message << std::endl;
    }
    else
    {
        std::cout << "✗ 求解失败" << std::endl;
        std::cout << "状态码: " << static_cast<int>(result.status) << std::endl;
        std::cout << "消息: " << result.message << std::endl;
    }
    std::cout << std::string(60, '=') << std::endl;
}

void example1_simpleSystem()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "示例1: 简单非线性系统" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "方程组: " << std::endl;
    std::cout << "  f1(x,y) = x^2 + y^2 - 4 = 0" << std::endl;
    std::cout << "  f2(x,y) = x - y = 0" << std::endl;
    std::cout << "预期解: x = ±√2, y = ±√2" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = x[0] * x[0] + x[1] * x[1] - PSGMDirectedInterval(4.0, 4.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = x[0] * PSGMDirectedInterval(2.0);
        J[0][1] = x[1] * PSGMDirectedInterval(2.0);
        J[1][0] = PSGMDirectedInterval(1.0);
        J[1][1] = PSGMDirectedInterval(-1.0);
        return J;
    };

    std::cout << "\n测试1.1: 初始区间 [-3, 3] × [-3, 3]" << std::endl;
    IntervalVector initial1(2);
    initial1[0] = PSGMDirectedInterval(-3.0, 3.0);
    initial1[1] = PSGMDirectedInterval(-3.0, 3.0);

    KrawczykSolver solver(2, f, jacobian, 1e-8, 50, 30);
    KrawczykResult result1 = solver.solve(initial1);
    printResult(result1);

    std::cout << "\n测试1.2: 初始区间 [0, 3] × [0, 3]" << std::endl;
    IntervalVector initial2(2);
    initial2[0] = PSGMDirectedInterval(0.0, 3.0);
    initial2[1] = PSGMDirectedInterval(0.0, 3.0);

    KrawczykResult result2 = solver.solve(initial2);
    printResult(result2);
}

void example2_polynomial()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "示例2: 多项式方程组" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "方程组: " << std::endl;
    std::cout << "  f1(x,y) = x^2 + y - 3 = 0" << std::endl;
    std::cout << "  f2(x,y) = x*y - 1 = 0" << std::endl;
    std::cout << "预期解: (1,2) 和 (-1,-2)" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = x[0] * x[0] + x[1] - PSGMDirectedInterval(3.0, 3.0);
        result[1] = x[0] * x[1] - PSGMDirectedInterval(1.0, 1.0);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = x[0] * PSGMDirectedInterval(2.0);
        J[0][1] = PSGMDirectedInterval(1.0);
        J[1][0] = x[1];
        J[1][1] = x[0];
        return J;
    };

    std::cout << "\n测试2.1: 初始区间 [-3, 3] × [-3, 3]" << std::endl;
    IntervalVector initial(2);
    initial[0] = PSGMDirectedInterval(-3.0, 3.0);
    initial[1] = PSGMDirectedInterval(-3.0, 3.0);

    KrawczykSolver solver(2, f, jacobian, 1e-8, 50, 30);
    KrawczykResult result = solver.solve(initial);
    printResult(result);

    std::cout << "\n测试2.2: 搜索所有解" << std::endl;
    auto allResults = solver.solveAll(initial, 5);
    std::cout << "找到 " << allResults.size() << " 个解" << std::endl;
    for (size_t i = 0; i < allResults.size(); ++i)
    {
        std::cout << "\n解 #" << (i + 1) << ":" << std::endl;
        for (size_t j = 0; j < allResults[i].solution.size(); ++j)
        {
            std::cout << "  x[" << j << "] = [" << std::fixed << std::setprecision(8)
                      << allResults[i].solution[j].lower() << ", "
                      << allResults[i].solution[j].upper() << "]" << std::endl;
        }
    }
}

void example3_cubic()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "示例3: 立方方程" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "方程: x^3 - x - 1 = 0" << std::endl;
    std::cout << "预期解: x ≈ 1.3247" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(1);
        result[0] = x[0] * x[0] * x[0] - x[0] - PSGMDirectedInterval(1.0, 1.0);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(1, 1);
        J[0][0] = x[0] * x[0] * PSGMDirectedInterval(3.0) - PSGMDirectedInterval(1.0);
        return J;
    };

    std::cout << "\n测试3.1: 初始区间 [0, 2]" << std::endl;
    IntervalVector initial(1);
    initial[0] = PSGMDirectedInterval(0.0, 2.0);

    KrawczykSolver solver(1, f, jacobian, 1e-10, 100, 20);
    KrawczykResult result = solver.solve(initial);
    printResult(result);
}

void example4_trigonometric()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "示例4: 三角函数方程" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "方程: sin(x) - x/2 = 0" << std::endl;
    std::cout << "预期解: x ≈ 0" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(1);
        result[0] = interval::sin(x[0]) - x[0] * PSGMDirectedInterval(0.5);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(1, 1);
        J[0][0] = interval::cos(x[0]) - PSGMDirectedInterval(0.5);
        return J;
    };

    std::cout << "\n测试4.1: 初始区间 [-π, π]" << std::endl;
    IntervalVector initial(1);
    initial[0] = PSGMDirectedInterval(-3.14159, 3.14159);

    KrawczykSolver solver(1, f, jacobian, 1e-10, 50, 20);
    KrawczykResult result = solver.solve(initial);
    printResult(result);
}

void example5_exponential()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "示例5: 指数方程组" << std::endl;
    std::cout << std::string(60, '*') << std::endl;
    std::cout << "方程组: " << std::endl;
    std::cout << "  f1(x,y) = e^x + y - 2 = 0" << std::endl;
    std::cout << "  f2(x,y) = x - e^y = 0" << std::endl;
    std::cout << "预期解: x ≈ 0.5, y ≈ 0.5" << std::endl;

    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = interval::exp(x[0]) + x[1] - PSGMDirectedInterval(2.0, 2.0);
        result[1] = x[0] - interval::exp(x[1]);
        return result;
    };

    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = interval::exp(x[0]);
        J[0][1] = PSGMDirectedInterval(1.0);
        J[1][0] = PSGMDirectedInterval(1.0);
        J[1][1] = -interval::exp(x[1]);
        return J;
    };

    std::cout << "\n测试5.1: 初始区间 [0, 1] × [0, 1]" << std::endl;
    IntervalVector initial(2);
    initial[0] = PSGMDirectedInterval(0.0, 1.0);
    initial[1] = PSGMDirectedInterval(0.0, 1.0);

    KrawczykSolver solver(2, f, jacobian, 1e-10, 50, 20);
    KrawczykResult result = solver.solve(initial);
    printResult(result);
}

void demonstrateIntervalArithmetic()
{
    std::cout << "\n" << std::string(60, '*') << std::endl;
    std::cout << "区间算术演示" << std::endl;
    std::cout << std::string(60, '*') << std::endl;

    PSGMDirectedInterval a(1.0, 3.0);
    PSGMDirectedInterval b(2.0, 4.0);

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;

    std::cout << "\n基本运算:" << std::endl;
    std::cout << "a + b = " << (a + b) << std::endl;
    std::cout << "a - b = " << (a - b) << std::endl;
    std::cout << "a * b = " << (a * b) << std::endl;
    std::cout << "a / b = " << (a / b) << std::endl;

    std::cout << "\n常用函数:" << std::endl;
    PSGMDirectedInterval x(0.5, 2.0);
    std::cout << "x = " << x << std::endl;
    std::cout << "x² = " << (x * x) << std::endl;
    std::cout << "√x = " << interval::sqrt(x) << std::endl;
    std::cout << "e^x = " << interval::exp(x) << std::endl;
    std::cout << "sin(x) = " << interval::sin(x) << std::endl;
    std::cout << "cos(x) = " << interval::cos(x) << std::endl;
}

int main()
{
    std::cout << std::string(70, '#') << std::endl;
    std::cout << "#" << std::setw(68) << std::left << " Interval-Krawczyk 方法求解非线性方程组" << "#" << std::endl;
    std::cout << std::string(70, '#') << std::endl;

    demonstrateIntervalArithmetic();

    try
    {
        example1_simpleSystem();
        example2_polynomial();
        example3_cubic();
        example4_trigonometric();
        example5_exponential();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n错误: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n" << std::string(70, '#') << std::endl;
    std::cout << "#" << std::setw(68) << std::left << " 所有示例完成" << "#" << std::endl;
    std::cout << std::string(70, '#') << std::endl;

    return 0;
}
