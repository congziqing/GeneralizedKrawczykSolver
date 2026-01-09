
# Interval Krawczyk: Generalized Solver with Interval Arithmetic

![Language](https://img.shields.io/badge/language-C%2B%2B17-blue.svg)
![License](https://img.shields.io/badge/license-0BSD-green.svg)
![Dependencies](https://img.shields.io/badge/dependencies-Eigen3-orange.svg)
![Build](https://img.shields.io/badge/build-passing-brightgreen.svg)

**A high-performance, header-only C++ library for rigorous numerical analysis using Generalized Krawczyk Operator and Kaucher Interval Arithmetic.**

---

## 📖 Why Interval Krawczyk?

In standard numerical computing, floating-point errors and sensitive initial guesses can lead to incorrect solutions or convergence failures. **Interval Krawczyk** solves this by working with *sets* of numbers rather than single points.

This library is designed for engineers and researchers who need **guarantees**. It allows you to:
*   **Prove the existence and uniqueness** of a solution within a given box.
*   Compute **AE-Solution Sets** (Universal and Existential quantifiers) for control systems under uncertainty.
*   Run on **embedded systems** thanks to zero dynamic memory allocation.

### Key Concepts

-   **Kaucher Interval Arithmetic:** Unlike classical interval arithmetic, Kaucher intervals (or directed intervals) allow "improper" intervals where the lower bound is greater than the upper bound. This creates a group structure, enabling better algebraic properties and the solution of problems involving quantifiers (e.g., $\forall x \exists y: f(x,y)=0$).
-   **Generalized Krawczyk Operator:** An extension of the Newton-Kantorovich method. It maps an interval vector to another. If the image is contained within the domain, it mathematically proves a solution exists inside.

---

## ✨ Features

*   🚀 **High Performance:** Built on **Eigen 3.4**, utilizing SIMD vectorization where possible.
*   💾 **Zero Dynamic Allocation:** Uses fixed-size `IntervalVector` and `IntervalMatrix` templates. deterministic memory usage makes it ideal for real-time and embedded applications.
*   🧮 **Kaucher Arithmetic:** Full support for both proper and improper intervals, enabling complex algebraic completions.
*   🔍 **AE-Solution Sets:** Capable of solving problems with mixed quantifiers ($\forall$ and $\exists$), critical for robust control and tolerance analysis.
*   🔄 **Homotopy Continuation:** Integrated methods to find multiple roots or track solutions as parameters change.
*   ✅ **Verified Solving:** Implementation of the Generalized Krawczyk method for rigorous nonlinear system solving.

---

## 🛠️ Technology Stack

*   **Language:** C++17 (Modern semantics and optimization)
*   **Math Backend:** Eigen 3.4.0 (Linear algebra)
*   **Architecture:** Header-only design (mostly) with fixed-size template metaprogramming.

---

## 📦 Installation

Since this is a header-only library dependent on Eigen, integration is straightforward.

```bash
# 1. Clone the repository
git clone https://github.com/congziqing/GeneralizedKrawczykSolver.git
cd GeneralizedKrawczykSolver

# 2. Run the tests to ensure your environment is ready
cd tests
# Make sure Eigen is in your include path
g++ -std=c++17 -I../include -I../eigen-3.4.0 -o FixedIntervalTest FixedIntervalTest.cpp
./FixedIntervalTest
```

---

## 🚀 Quick Start

Here is a complete example of solving a nonlinear system.
We are solving for:
1. $x^2 + y^2 = 4$ (Circle)
2. $x - y = 0$ (Line)

The solver will rigorously find the intersection points.

```cpp
#include <iostream>
#include "include/interval_krawczyk/KaucherInterval.h"
#include "include/interval_krawczyk/IntervalVector.h"
#include "include/interval_krawczyk/GeneralizedKrawczykSolver.h"

using namespace ik;

int main() {
    // 1. Define the Nonlinear System F(x)
    auto f = [](const IntervalVector<2>& x) -> IntervalVector<2> {
        IntervalVector<2> result;
        // Equation 1: x^2 + y^2 - 4 = 0
        result[0] = x[0] * x[0] + x[1] * x[1] - KaucherInterval(4.0, 4.0);
        // Equation 2: x - y = 0
        result[1] = x[0] - x[1];
        return result;
    };

    // 2. Define the Jacobian Matrix J(x)
    auto jacobian = [](const IntervalVector<2>& x) -> IntervalMatrix<2, 2> {
        IntervalMatrix<2, 2> J;
        // Partial derivatives
        J(0, 0) = KaucherInterval(2.0) * x[0]; // d(eq1)/dx
        J(0, 1) = KaucherInterval(2.0) * x[1]; // d(eq1)/dy
        J(1, 0) = KaucherInterval(1.0);        // d(eq2)/dx
        J(1, 1) = KaucherInterval(-1.0);       // d(eq2)/dy
        return J;
    };

    // 3. Initialize Solver and Search Box
    GeneralizedKrawczykSolver<2> solver(f, jacobian);
    
    IntervalVector<2> initialBox;
    // We suspect a solution is near (1.4, 1.4), so we search in [1, 2]
    initialBox[0] = KaucherInterval(1.0, 2.0);
    initialBox[1] = KaucherInterval(1.0, 2.0);

    // 4. Solve
    auto result = solver.solve(initialBox);

    if (result.success) {
        std::cout << "✅ Verified Solution found within:\n" << result.solution << std::endl;
    } else {
        std::cout << "No unique solution verified in this box." << std::endl;
    }

    return 0;
}
```

---

## 📂 Project Structure

```text
Interval-Krawczyk/
├── include/
│   └── interval_krawczyk/
│       ├── KaucherInterval.h           # Core arithmetic for Directed Intervals
│       ├── IntervalVector.h            # Static size vector wrapper
│       ├── IntervalMatrix.h            # Static size matrix wrapper
│       └── GeneralizedKrawczykSolver.h # The main solver algorithm
├── tests/                              # Unit tests and validation
├── demos/                              # Advanced usage examples
├── eigen-3.4.0/                        # Eigen dependency
└── README.md
```

---

## 🤝 Contribution

We welcome contributions from the community! Whether it's adding new test cases, optimizing the arithmetic, or improving documentation.

1.  Fork the Project
2.  Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3.  Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4.  Push to the Branch (`git push origin feature/AmazingFeature`)
5.  Open a Pull Request

Please ensure all new code follows the existing style and passes all tests in the `tests/` directory.

---

## 📄 License

Distributed under the **0BSD License**. This means you can use it for whatever you want, commercial or personal, without restriction. See `LICENSE` for more information.

---

## 📧 Contact & Acknowledgments

**Maintainer:** [Cong Ziqing](mailto:congziqing@126.com)  
**Project Link:** [https://github.com/congziqing/GeneralizedKrawczykSolver](https://github.com/congziqing/GeneralizedKrawczykSolver)

**Special Thanks:**
*   [Eigen Library](https://eigen.tuxfamily.org) for the robust linear algebra backend.
