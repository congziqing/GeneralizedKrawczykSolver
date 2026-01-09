#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"
#include "../include/interval_krawczyk/IntervalVector.h"

using namespace ik;

int main() {
    IntervalVector<3> v1;
    v1[0] = KaucherInterval(1.0, 2.0);
    v1[1] = KaucherInterval(3.0, 4.0);
    v1[2] = KaucherInterval(5.0, 6.0);
    
    std::cout << "v1 = " << v1 << std::endl;
    std::cout << "v1.midpoint() = " << v1.midpoint().transpose() << std::endl;
    std::cout << "v1.maxWidth() = " << v1.maxWidth() << std::endl;
    
    return 0;
}