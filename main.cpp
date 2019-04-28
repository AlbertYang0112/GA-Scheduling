#include <iostream>
#include <iomanip>
#include "NaiveGA.h"

void NaiveGADemo() {
    /* Population Size 100
     * Cross Rate: 10%
     * Mutation Rate: 10%
     */
    NaiveGA tooYoung(10000, 0.1, 0.1);
    uint32_t bestGene;
    double_t maxFit;

    tooYoung.evaluate(1000000, bestGene, maxFit);

    std::cout << "Best Gene "
              << std::hex << std::setfill('0') << std::setw(8)
              << bestGene << std::endl;
    std::cout << "Best Fitness "
              << std::dec
              << maxFit << std::endl;
}

int main() {
    NaiveGADemo();

    return 0;
}