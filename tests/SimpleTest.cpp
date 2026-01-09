#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "../include/interval_krawczyk/KaucherInterval.h"

using namespace ik;

int main() {
    KaucherInterval a(1.0, 2.0);
    std::cout << "a = " << a << std::endl;
    return 0;
}