#include "include/interval_krawczyk/KaucherInterval.h"
#include <iostream>

int main() {
    // 创建一个简单的KaucherInterval对象
    ik::KaucherInterval a(1.0, 2.0);
    ik::KaucherInterval b(3.0, 4.0);
    
    // 测试基本运算
    ik::KaucherInterval c = a + b;
    ik::KaucherInterval d = a * b;
    
    // 输出结果
    std::cout << "a = [" << a.lower() << ", " << a.upper() << "]" << std::endl;
    std::cout << "b = [" << b.lower() << ", " << b.upper() << "]" << std::endl;
    std::cout << "a + b = [" << c.lower() << ", " << c.upper() << "]" << std::endl;
    std::cout << "a * b = [" << d.lower() << ", " << d.upper() << "]" << std::endl;
    
    // 测试非正常区间
    ik::KaucherInterval improper(2.0, 1.0);
    std::cout << "Improper interval: [" << improper.lower() << ", " << improper.upper() << "]" << std::endl;
    std::cout << "Dual: [" << improper.dual().lower() << ", " << improper.dual().upper() << "]" << std::endl;
    
    return 0;
}