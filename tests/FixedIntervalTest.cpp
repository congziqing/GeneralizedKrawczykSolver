#define _USE_MATH_DEFINES

#include <iostream>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"
#include "../include/interval_krawczyk/IntervalMatrix.h"
using namespace ik;

int main()
{
    std::cout << "=== Fixed-Dimension Interval Classes Test ===\n\n";

    std::cout << "1. IntervalVector Test:\n";
    IntervalVector<3> v1;
    v1[0] = KaucherInterval(1.0, 2.0);
    v1[1] = KaucherInterval(3.0, 4.0);
    v1[2] = KaucherInterval(5.0, 6.0);
    std::cout << "v1 = [[" << v1[0] << ", " << v1[1] << ", " << v1[2] << "]]" << std::endl;
    std::cout << "v1.midpoint() = " << v1.midpoint().transpose() << std::endl;
    std::cout << "v1.maxWidth() = " << v1.maxWidth() << std::endl;

    IntervalVector<3> v2 = IntervalVector<3>::constants(1.0);
    std::cout << "v2 (constants 1.0) = [[" << v2[0] << ", " << v2[1] << ", " << v2[2] << "]]" << std::endl;

    auto v3 = v1 + v2;
    std::cout << "v1 + v2 = [[" << v3[0] << ", " << v3[1] << ", " << v3[2] << "]]" << std::endl;

    auto v4 = v1 * 2.0;
    std::cout << "v1 * 2.0 = [[" << v4[0] << ", " << v4[1] << ", " << v4[2] << "]]" << std::endl;

    std::cout << "\n2. IntervalMatrix Test:\n";
    IntervalMatrix<2, 2> M1;
    M1(0, 0) = KaucherInterval(1.0, 2.0);
    M1(0, 1) = KaucherInterval(3.0, 4.0);
    M1(1, 0) = KaucherInterval(5.0, 6.0);
    M1(1, 1) = KaucherInterval(7.0, 8.0);
    std::cout << "M1 = " << std::endl;
    std::cout << "[[" << M1(0, 0) << ", " << M1(0, 1) << "]" << std::endl;
    std::cout << " [" << M1(1, 0) << ", " << M1(1, 1) << "]]" << std::endl;
    std::cout << "M1.midpoint() = " << std::endl;
    std::cout << M1.midpoint() << std::endl;

    // 手动创建单位矩阵
    IntervalMatrix<2, 2> M2;
    M2(0, 0) = KaucherInterval(1.0, 1.0);
    M2(0, 1) = KaucherInterval(0.0, 0.0);
    M2(1, 0) = KaucherInterval(0.0, 0.0);
    M2(1, 1) = KaucherInterval(1.0, 1.0);
    std::cout << "M2 (identity) = " << std::endl;
    std::cout << "[[" << M2(0, 0) << ", " << M2(0, 1) << "]" << std::endl;
    std::cout << " [" << M2(1, 0) << ", " << M2(1, 1) << "]]" << std::endl;

    auto M3 = M1 + M2;
    std::cout << "M1 + M2 = " << std::endl;
    std::cout << "[[" << M3(0, 0) << ", " << M3(0, 1) << "]" << std::endl;
    std::cout << " [" << M3(1, 0) << ", " << M3(1, 1) << "]]" << std::endl;

    std::cout << "\n3. Matrix-Vector Multiplication Test:\n";
    IntervalVector<2> vec;
    vec[0] = KaucherInterval(1.0, 1.0);
    vec[1] = KaucherInterval(2.0, 2.0);
    std::cout << "vec = [[" << vec[0] << ", " << vec[1] << "]]" << std::endl;

    auto result = M1 * vec;
    std::cout << "M1 * vec = [[" << result[0] << ", " << result[1] << "]]" << std::endl;

    std::cout << "\n4. Performance Comparison:\n";
    std::cout << "sizeof(KaucherInterval) = " << sizeof(KaucherInterval) << " bytes\n";
    std::cout << "sizeof(IntervalVector<4>) = " << sizeof(IntervalVector<4>) << " bytes\n";
    std::cout << "sizeof(IntervalMatrix<4,4>) = " << sizeof(IntervalMatrix<4,4>) << " bytes\n";
    std::cout << "No dynamic allocation - all on stack!\n";

    std::cout << "\n=== Test Complete ===\n";

    return 0;
}