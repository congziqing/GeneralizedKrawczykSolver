# 区间-Krawczyk方法求解非线性方程组库

## 1. 库介绍

本库实现了基于区间算术和Krawczyk方法的非线性方程组求解器，结合同伦连续法解决复杂问题。支持Kaucher区间算术（广义区间算术），适用于固定维和动态维数的问题。

### 1.1 主要特性

- ✅ **区间算术**：支持Kaucher区间算术，处理proper和improper区间
- ✅ **动态维数支持**：支持动态大小的区间向量和矩阵
- ✅ **广义Krawczyk方法**：实现广义Krawczyk迭代，支持非线性方程组求解
- ✅ **同伦连续法**：结合同伦连续法，解决复杂方程组和避免迭代停滞
- ✅ **多起点求解**：支持从多个起点求解，找到所有根
- ✅ **精确结果**：提供区间形式的解，保证结果的可靠性
- ✅ **Eigen集成**：与Eigen库无缝集成，提供高效的线性代数运算

### 1.2 应用领域

- 非线性方程组求解
- 复数多项式求根
- 机器人运动学
- 控制系统设计
- 优化问题
- 科学计算

## 2. 安装与配置

### 2.1 依赖项

- C++17或更高版本
- [Eigen 3.4.0](http://eigen.tuxfamily.org/)：用于线性代数运算

### 2.2 安装

1. 克隆或下载本项目
2. 确保Eigen库位于项目根目录下的`eigen-3.4.0`文件夹
3. 编译时添加Eigen库路径：`-I./eigen-3.4.0`

### 2.3 目录结构

```
├── include/            # 核心头文件
│   └── interval_krawczyk/   # 区间-Krawczyk相关实现
│       ├── PSGMDirectedInterval.h       # Kaucher区间算术类
│       ├── IntervalVector.h             # 区间向量类
│       ├── IntervalMatrix.h             # 区间矩阵类
│       └── GeneralizedKrawczykSolver.h  # 广义Krawczyk求解器
├── src/                # 其他实现文件（已精简）
├── demos/              # 演示代码
│   ├── ComplexPolynomialDemo.cpp    # 复数多项式求解演示
│   ├── GeneralizedKrawczykDemo.cpp  # 广义Krawczyk方法演示
│   └── ...
├── tests/              # 测试代码
│   ├── KaucherArithmeticTest.cpp    # Kaucher算术测试
│   └── ...
└── eigen-3.4.0/        # Eigen库
```

## 3. 核心概念

### 3.1 区间算术

区间算术是一种数学工具，用区间表示不确定性。对于实数区间`[a, b]`，其中`a ≤ b`，支持基本算术运算：

- 加法：`[a, b] + [c, d] = [a+c, b+d]`
- 减法：`[a, b] - [c, d] = [a-d, b-c]`
- 乘法：`[a, b] × [c, d] = [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]`
- 除法：`[a, b] / [c, d] = [a, b] × [1/d, 1/c]`（当0不在`[c, d]`中）

### 3.2 Kaucher区间算术

Kaucher区间算术扩展了传统区间算术，允许`a > b`的improper区间。这对于内包围求解和验证非常有用。

### 3.3 Krawczyk方法

Krawczyk方法是一种用于求解非线性方程组的迭代方法，基于区间算术和牛顿法。对于方程组`f(x) = 0`，Krawczyk迭代公式为：

```
K(X) = c - Yf(c) + (I - YJ(X))(X - c)
```

其中：
- `c`是区间`X`的中心点
- `Y`是`J(c)`的近似逆矩阵
- `J(X)`是`f`在`X`上的区间雅可比矩阵
- `I`是单位矩阵

### 3.4 同伦连续法

同伦连续法通过构造连接简单问题和目标问题的同伦函数，逐步求解。对于目标函数`f(x)`和辅助函数`g(x)`，同伦函数定义为：

```
H(x, t) = t·f(x) + (1-t)·g(x)
```

将`t`从0逐步增加到1，跟踪解的路径，最终得到目标问题的解。

## 4. 快速入门

### 4.1 简单非线性方程组求解

```cpp
#include "../include/interval_krawczyk/PSGMDirectedInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
#include "../include/interval_krawczyk/GeneralizedKrawczykSolver.h"

using namespace std;

int main() {
    // 定义方程组：x² + y = 2, x - y = 0
    auto f = [](const IntervalVector& x) -> IntervalVector {
        IntervalVector result(2);
        result[0] = x[0] * x[0] + x[1] - PSGMDirectedInterval(2.0, 2.0);
        result[1] = x[0] - x[1];
        return result;
    };

    // 定义雅可比矩阵
    auto jacobian = [](const IntervalVector& x) -> IntervalMatrix {
        IntervalMatrix J(2, 2);
        J[0][0] = PSGMDirectedInterval(2.0) * x[0];
        J[0][1] = PSGMDirectedInterval(1.0, 1.0);
        J[1][0] = PSGMDirectedInterval(1.0, 1.0);
        J[1][1] = PSGMDirectedInterval(-1.0, -1.0);
        return J;
    };

    // 创建Krawczyk求解器
    GeneralizedKrawczykSolver solver(2, f, jacobian, 1e-8, 100);

    // 初始区间
    IntervalVector initialBox(2);
    initialBox[0] = PSGMDirectedInterval(0.0, 2.0);
    initialBox[1] = PSGMDirectedInterval(0.0, 2.0);

    // 求解
    auto result = solver.solve(initialBox);

    if (result.success) {
        cout << "求解成功！" << endl;
        cout << "解区间：" << result.solution << endl;
        cout << "区间宽度：" << result.finalWidth << endl;
    } else {
        cout << "求解失败：" << result.message << endl;
    }

    return 0;
}
```

### 4.2 复数多项式求解

```cpp
#include "../include/interval_krawczyk/PSGMDirectedInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
#include "../include/interval_krawczyk/GeneralizedKrawczykSolver.h"

// 将复数多项式z³ - 1 = 0转换为实系统
// 1. x³ - 3xy² - 1 = 0
// 2. 3x²y - y³ = 0

auto f = [](const IntervalVector& x) -> IntervalVector {
    IntervalVector result(2);
    PSGMDirectedInterval x_val = x[0];
    PSGMDirectedInterval y_val = x[1];
    
    result[0] = x_val * x_val * x_val - PSGMDirectedInterval(3.0) * x_val * y_val * y_val - PSGMDirectedInterval(1.0, 1.0);
    result[1] = PSGMDirectedInterval(3.0) * x_val * x_val * y_val - y_val * y_val * y_val;
    return result;
};

// 使用同伦连续法求解...
```

## 5. 进阶用法

### 5.1 同伦连续法

```cpp
#include "../include/interval_krawczyk/PSGMDirectedInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
#include "../include/interval_krawczyk/GeneralizedKrawczykSolver.h"

// 创建同伦连续求解器
SimpleHomotopySolver solver(n, f, jacobian, point_f, point_jacobian);

// 三个不同的起点
vector<vector<double>> start_points = {
    {1.0, 0.0},
    {-0.5, 0.8},
    {-0.5, -0.8}
};

// 求解
for (const auto& start : start_points) {
    auto result = solver.solve(start, 10);
    // 处理结果...
}
```



## 6. API参考

### 6.1 核心类

#### 6.1.1 PSGMDirectedInterval

区间算术的核心类，支持Kaucher区间算术。

```cpp
class PSGMDirectedInterval {
public:
    // 构造函数
    PSGMDirectedInterval(double lower, double upper);
    PSGMDirectedInterval(double value);
    
    // 算术运算
    PSGMDirectedInterval operator+(const PSGMDirectedInterval& other) const;
    PSGMDirectedInterval operator-(const PSGMDirectedInterval& other) const;
    PSGMDirectedInterval operator*(const PSGMDirectedInterval& other) const;
    PSGMDirectedInterval operator/(const PSGMDirectedInterval& other) const;
    
    // 辅助方法
    bool isEmpty() const;
    bool isProper() const;
    double mid() const;
    double width() const;
    double lower() const;
    double upper() const;
};
```

#### 6.1.2 IntervalVector

区间向量类，支持动态维数。

```cpp
class IntervalVector {
public:
    // 构造函数
    IntervalVector(size_t size);
    
    // 访问器
    PSGMDirectedInterval& operator[](size_t index);
    const PSGMDirectedInterval& operator[](size_t index) const;
    
    // 辅助方法
    size_t size() const;
    bool empty() const;
    double width() const;
    Eigen::VectorXd midpoint() const;
    
    // 算术运算
    IntervalVector operator+(const IntervalVector& other) const;
    IntervalVector operator-(const IntervalVector& other) const;
    IntervalVector operator*(const PSGMDirectedInterval& scalar) const;
};
```

#### 6.1.3 IntervalMatrix

区间矩阵类，支持动态维数。

```cpp
class IntervalMatrix {
public:
    // 构造函数
    IntervalMatrix(size_t rows, size_t cols);
    
    // 访问器
    PSGMDirectedInterval& operator()(size_t row, size_t col);
    const PSGMDirectedInterval& operator()(size_t row, size_t col) const;
    PSGMDirectedInterval* operator[](size_t row);
    const PSGMDirectedInterval* operator[](size_t row) const;
    
    // 辅助方法
    size_t rows() const;
    size_t cols() const;
    
    // 算术运算
    IntervalMatrix operator+(const IntervalMatrix& other) const;
    IntervalMatrix operator-(const IntervalMatrix& other) const;
    IntervalMatrix operator*(const IntervalMatrix& other) const;
    IntervalVector operator*(const IntervalVector& vec) const;
};
```

#### 6.1.4 GeneralizedKrawczykSolver

广义Krawczyk求解器。

```cpp
class GeneralizedKrawczykSolver {
public:
    struct Result {
        bool success;
        IntervalVector solution;
        int iterations;
        double finalWidth;
        std::string message;
    };
    
    // 构造函数
    GeneralizedKrawczykSolver(size_t dimension, 
                             const std::function<IntervalVector(const IntervalVector&)>& f, 
                             const std::function<IntervalMatrix(const IntervalVector&)>& jacobian, 
                             double tolerance = 1e-8, 
                             int maxIterations = 100);
    
    // 求解方法
    Result solve(const IntervalVector& initialBox);
};
```



## 7. 最佳实践

### 7.1 初始区间选择

- 对于简单问题，选择包含解的合理初始区间
- 对于复杂问题，使用同伦连续法，从简单问题逐步过渡到目标问题
- 对于多根问题，使用多个不同的初始起点

### 7.2 收敛条件设置

- 容差（tolerance）：根据问题精度要求设置，建议1e-8到1e-12
- 最大迭代次数：根据问题复杂度设置，建议100到1000
- 对于同伦连续法，步数建议10到50步

### 7.3 性能优化

- 预计算常数和矩阵，减少重复计算
- 优化雅可比矩阵计算，尽量使用区间算术的高效实现

### 7.4 结果验证

- 检查解区间是否包含实际解
- 对于关键应用，使用内包围验证
- 比较不同方法和参数的结果，提高可靠性

## 8. 演示程序

### 8.1 复数多项式求解

```bash
g++ -std=c++17 -I./eigen-3.4.0 ./demos/ComplexPolynomialDemo.cpp -o ComplexPolynomialDemo.exe
./ComplexPolynomialDemo.exe
```

### 8.2 同伦连续法演示

```bash
g++ -std=c++17 -I./eigen-3.4.0 ./demos/HomotopyDemoSimple.cpp -o HomotopyDemoSimple.exe
./HomotopyDemoSimple.exe
```

### 8.3 广义Krawczyk方法演示

```bash
g++ -std=c++17 -I./eigen-3.4.0 ./demos/GeneralizedKrawczykDemo.cpp -o GeneralizedKrawczykDemo.exe
./GeneralizedKrawczykDemo.exe
```

## 9. 测试程序

### 9.1 Kaucher算术测试

```bash
g++ -std=c++17 -I./eigen-3.4.0 ./tests/KaucherArithmeticTest.cpp -o KaucherArithmeticTest.exe
./KaucherArithmeticTest.exe
```

### 9.2 扩张映射测试

```bash
g++ -std=c++17 -I./eigen-3.4.0 ./tests/ExpandingMapTest.cpp -o ExpandingMapTest.exe
./ExpandingMapTest.exe
```

### 9.3 区间向量和矩阵测试

```bash
g++ -std=c++17 -I./eigen-3.4.0 ./tests/FixedIntervalTest.cpp -o FixedIntervalTest.exe
./FixedIntervalTest.exe
```

## 10. 常见问题解答

### 10.1 Krawczyk迭代不收敛怎么办？

- 检查初始区间是否包含解
- 尝试使用同伦连续法
- 调整容差和最大迭代次数
- 验证雅可比矩阵计算是否正确

### 10.2 如何处理多根问题？

- 使用多个不同的初始起点
- 结合同伦连续法，从不同起点追踪不同的解路径
- 对于多项式问题，可以使用复数平面的不同区域作为起点

### 10.3 如何提高求解速度？

- 使用固定维数类模板
- 优化雅可比矩阵计算
- 减少同伦连续法的步数
- 调整收敛条件，适当降低精度要求

### 10.4 如何验证解的正确性？

- 将解代入原方程组，检查是否满足f(x) ≈ 0
- 使用内包围验证，确保解严格位于期望区间内
- 比较不同方法和参数的结果

## 11. 扩展阅读

- [Interval Arithmetic: Introduction and Applications](https://doi.org/10.1007/978-3-642-16216-1)
- [Kaucher Interval Arithmetic](https://en.wikipedia.org/wiki/Kaucher_interval_arithmetic)
- [Krawczyk_method](https://en.wikipedia.org/wiki/Krawczyk_method)
- [Homotopy Continuation](https://en.wikipedia.org/wiki/Homotopy_continuation_method)
- [Eigen Library Documentation](https://eigen.tuxfamily.org/dox/)

## 12. 贡献与反馈

欢迎提交issue和pull request！如有问题或建议，请联系：

- 项目地址：[GitHub Repository]
- 邮箱：[your-email@example.com]

## 13. 许可证

本项目采用MIT许可证，详见LICENSE文件。

---

**从入门到精通，掌握区间-Krawczyk方法求解非线性方程组！** 🚀