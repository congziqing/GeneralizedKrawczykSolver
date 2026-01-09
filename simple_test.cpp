#include "include/interval_krawczyk/KaucherInterval.h"
#include <iostream>

int main() {
    // Create a simple KaucherInterval object
    ik::KaucherInterval a(1.0, 2.0);
    ik::KaucherInterval b(3.0, 4.0);
    
    // Test basic operations
    ik::KaucherInterval c = a + b;
    ik::KaucherInterval d = a * b;
    
    // Output results
    std::cout << "a = [" << a.lower() << ", " << a.upper() << "]" << std::endl;
    std::cout << "b = [" << b.lower() << ", " << b.upper() << "]" << std::endl;
    std::cout << "a + b = [" << c.lower() << ", " << c.upper() << "]" << std::endl;
    std::cout << "a * b = [" << d.lower() << ", " << d.upper() << "]" << std::endl;
    
    // Test improper interval
    ik::KaucherInterval improper(2.0, 1.0);
    std::cout << "Improper interval: [" << improper.lower() << ", " << improper.upper() << "]" << std::endl;
    std::cout << "Dual: [" << improper.dual().lower() << ", " << improper.dual().upper() << "]" << std::endl;
    
    return 0;
}