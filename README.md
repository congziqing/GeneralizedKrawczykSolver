# Interval Krawczyk - Generalized Krawczyk Solver Based on Interval Arithmetic

An efficient C++ library implementing the generalized Krawczyk method with interval arithmetic and Kaucher intervals for solving nonlinear equations and AE solution sets.

## Features

- ✅ Complete arithmetic operations for both proper and improper intervals (Kaucher intervals)
- ✅ Generalized Krawczyk method for solving nonlinear equations
- ✅ Support for AE solution set computation with universal and existential quantifiers
- ✅ Fixed-size interval vector and matrix classes with no dynamic memory allocation
- ✅ Homotopy continuation method for finding multiple solutions
- ✅ Rich test cases and demonstration programs

## Technology Stack

- **Programming Language**: C++17
- **Core Library**: Eigen 3.4.0
- **Build System**: GCC/Clang/MSVC compilers
- **Testing Framework**: Custom test cases

## Installation and Usage

```bash
# Clone the repository
git clone https://github.com/yourusername/Interval-Krawczyk.git
cd Interval-Krawczyk

# Compile test program
cd tests
g++ -std=c++17 -I../include -I../eigen-3.4.0 -o FixedIntervalTest FixedIntervalTest.cpp
./FixedIntervalTest
```

## Quick Start Example

```cpp
#include <iostream>
#include "include/interval_krawczyk/KaucherInterval.h"
#include "include/interval_krawczyk/IntervalVector.h"
#include "include/interval_krawczyk/GeneralizedKrawczykSolver.h"

using namespace ik;

int main() {
    // Define nonlinear system and Jacobian
    auto f = [](const IntervalVector<2>& x) -> IntervalVector<2> {
        IntervalVector<2> result;
        result[0] = x[0] * x[0] + x[1] * x[1] - KaucherInterval(4.0, 4.0);
        result[1] = x[0] - x[1];
        return result;
    };

    auto jacobian = [](const IntervalVector<2>& x) -> IntervalMatrix<2, 2> {
        IntervalMatrix<2, 2> J;
        J(0, 0) = KaucherInterval(2.0) * x[0];
        J(0, 1) = KaucherInterval(2.0) * x[1];
        J(1, 0) = KaucherInterval(1.0);
        J(1, 1) = KaucherInterval(-1.0);
        return J;
    };

    // Create solver and solve
    GeneralizedKrawczykSolver<2> solver(f, jacobian);
    IntervalVector<2> initialBox;
    initialBox[0] = KaucherInterval(1.0, 2.0);
    initialBox[1] = KaucherInterval(1.0, 2.0);

    auto result = solver.solve(initialBox);
    if (result.success) {
        std::cout << "Solution found: " << result.solution << std::endl;
    }

    return 0;
}
```

## Project Structure

```
Interval-Krawczyk/
├── include/
│   └── interval_krawczyk/
│       ├── KaucherInterval.h        # Kaucher interval class
│       ├── IntervalVector.h         # Fixed-size interval vector class
│       ├── IntervalMatrix.h         # Fixed-size interval matrix class
│       └── GeneralizedKrawczykSolver.h  # Generalized Krawczyk solver
├── tests/                          # Test cases
├── demos/                          # Demonstration programs
├── eigen-3.4.0/                    # Eigen library dependency
├── LICENSE                         # 0BSD License
└── README.md                       # Project documentation
```

## Contribution Guidelines

Contributions are welcome! Please ensure:
1. Follow the project's code style
2. Add appropriate test cases
3. Update relevant documentation
4. Run all tests before submission

## License

This project is licensed under the 0BSD License - see the [LICENSE](LICENSE) file for details.

## Contact Information

- Project Link: https://github.com/yourusername/Interval-Krawczyk