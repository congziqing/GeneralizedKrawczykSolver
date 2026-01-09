
# Interval Krawczyk: 基于区间算术的广义求解器

![Language](https://img.shields.io/badge/language-C%2B%2B17-blue.svg)
![License](https://img.shields.io/badge/license-0BSD-green.svg)
![Dependencies](https://img.shields.io/badge/dependencies-Eigen3-orange.svg)
![Build](https://img.shields.io/badge/build-passing-brightgreen.svg)

**一个基于广义 Krawczyk 算子和 Kaucher 区间算术的高性能 C++ 头文件库，专为严谨数值分析设计。**

---

## 📖 为什么选择 Interval Krawczyk?

在标准的数值计算中，浮点误差和对初值敏感的特性往往会导致解不准确甚至无法收敛。**Interval Krawczyk** 通过处理**数集（区间）**而非单点数值，从根本上解决了这一问题。

本库专为需要**数学保证（Guarantees）**的工程师和研究人员设计。它能助您实现：
*   **严格证明**在给定范围内解的存在性和唯一性。
*   计算不确定性控制系统中的 **AE-解集**（包含全称量词 $\forall$ 和存在量词 $\exists$）。
*   得益于零动态内存分配，可直接在 **嵌入式系统** 或实时环境上运行。

### 核心概念

-   **Kaucher 区间算术 (Kaucher Interval Arithmetic):** 与经典区间算术不同，Kaucher 区间（或称定向区间）允许“非正常”区间的存在（即下界大于上界）。这构建了一个完整的代数群结构，赋予了更好的代数性质，并支持求解包含复杂量词的问题（例如 $\forall x \exists y: f(x,y)=0$）。
-   **广义 Krawczyk 算子 (Generalized Krawczyk Operator):** 这是 Newton-Kantorovich 方法的区间扩展。它将一个区间向量映射到另一个区间向量。如果映射后的图像被包含在原定义域内，即可在数学上严格证明该区域内存在解。

---

## ✨ 功能特性

*   🚀 **高性能:** 基于 **Eigen 3.4** 构建，尽可能利用 SIMD 指令集加速。
*   💾 **零动态内存分配:** 使用定长模板类 `IntervalVector` 和 `IntervalMatrix`。确定的内存使用模式使其非常适合实时系统和嵌入式应用。
*   🧮 **Kaucher 算术:** 完整支持正常和非正常区间运算，实现了复杂的代数补全。
*   🔍 **AE-解集计算:** 能够求解混合量词（$\forall$ 和 $\exists$）问题，这对鲁棒控制和公差分析至关重要。
*   🔄 **同伦延拓法 (Homotopy Continuation):** 内置方法用于寻找多重根或追踪随参数变化的解路径。
*   ✅ **验证求解:** 实现了广义 Krawczyk 方法，用于严谨非线性系统求解。

---

## 🛠️ 技术栈

*   **编程语言:** C++17 (现代语义与优化)
*   **数学后端:** Eigen 3.4.0 (线性代数支持)
*   **架构设计:** Header-only (纯头文件)，大量使用定长模板元编程。

---

## 📦 安装指南

由于这是一个依赖 Eigen 的纯头文件库，集成非常简单。

```bash
# 1. 克隆仓库
git clone https://github.com/congziqing/GeneralizedKrawczykSolver.git
cd GeneralizedKrawczykSolver

# 2. 编译并运行测试（确保环境配置正确）
cd tests
# 请确保 Eigen 在您的包含路径中
g++ -std=c++17 -I../include -I../eigen-3.4.0 -o FixedIntervalTest FixedIntervalTest.cpp
./FixedIntervalTest
```

---

## 🚀 快速开始

下面是一个求解非线性系统的完整示例。
我们将求解以下方程组的交点：
1. $x^2 + y^2 = 4$ (圆)
2. $x - y = 0$ (直线)

求解器将严格地找到交点区间。

```cpp
#include <iostream>
#include "include/interval_krawczyk/KaucherInterval.h"
#include "include/interval_krawczyk/IntervalVector.h"
#include "include/interval_krawczyk/GeneralizedKrawczykSolver.h"

using namespace ik;

int main() {
    // 1. 定义非线性系统 F(x)
    auto f = [](const IntervalVector<2>& x) -> IntervalVector<2> {
        IntervalVector<2> result;
        // 方程 1: x^2 + y^2 - 4 = 0
        result[0] = x[0] * x[0] + x[1] * x[1] - KaucherInterval(4.0, 4.0);
        // 方程 2: x - y = 0
        result[1] = x[0] - x[1];
        return result;
    };

    // 2. 定义雅可比矩阵 J(x)
    auto jacobian = [](const IntervalVector<2>& x) -> IntervalMatrix<2, 2> {
        IntervalMatrix<2, 2> J;
        // 偏导数计算
        J(0, 0) = KaucherInterval(2.0) * x[0]; // d(eq1)/dx
        J(0, 1) = KaucherInterval(2.0) * x[1]; // d(eq1)/dy
        J(1, 0) = KaucherInterval(1.0);        // d(eq2)/dx
        J(1, 1) = KaucherInterval(-1.0);       // d(eq2)/dy
        return J;
    };

    // 3. 初始化求解器和搜索框(Box)
    GeneralizedKrawczykSolver<2> solver(f, jacobian);
    
    IntervalVector<2> initialBox;
    // 我们猜测解在 (1.4, 1.4) 附近，因此在 [1, 2] 区间内搜索
    initialBox[0] = KaucherInterval(1.0, 2.0);
    initialBox[1] = KaucherInterval(1.0, 2.0);

    // 4. 执行求解
    auto result = solver.solve(initialBox);

    if (result.success) {
        std::cout << "✅ 已验证解存在于以下区间内:\n" << result.solution << std::endl;
    } else {
        std::cout << "在该搜索框内未验证到唯一解。" << std::endl;
    }

    return 0;
}
```

---

## 📂 项目结构

```text
Interval-Krawczyk/
├── include/
│   └── interval_krawczyk/
│       ├── KaucherInterval.h           # 定向区间(Directed Interval)的核心算术实现
│       ├── IntervalVector.h            # 定长向量封装类
│       ├── IntervalMatrix.h            # 定长矩阵封装类
│       └── GeneralizedKrawczykSolver.h # 主要求解器算法
├── tests/                              # 单元测试与验证用例
├── demos/                              # 高级演示程序
├── eigen-3.4.0/                        # Eigen 依赖库
└── README.md
```

---

## 🤝 贡献指南

我们非常欢迎来自社区的贡献！无论是添加新的测试用例、优化算术性能，还是改进文档。

1.  Fork 本项目
2.  创建您的特性分支 (`git checkout -b feature/AmazingFeature`)
3.  提交您的修改 (`git commit -m 'Add some AmazingFeature'`)
4.  推送到分支 (`git push origin feature/AmazingFeature`)
5.  提交 Pull Request

请确保所有新代码遵循现有的编码风格，并通过 `tests/` 目录下的所有测试。

---

## 📄 开源协议

本项目基于 **0BSD License** 分发。这意味着您可以将其用于任何目的（包括商业或个人用途），且没有任何限制。详情请参阅 `LICENSE` 文件。

---

## 📧 联系方式与致谢

**维护者:** [Cong Ziqing](mailto:congziqing@126.com)  
**项目链接:** [https://github.com/congziqing/GeneralizedKrawczykSolver](https://github.com/congziqing/GeneralizedKrawczykSolver)

**特别致谢:**
*   [Eigen Library](https://eigen.tuxfamily.org): 提供了健壮的线性代数后端支持。
