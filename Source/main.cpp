#include <iostream>
#include <iomanip>
#include "Descriptor.h"
#include "NaiveGA.h"
#include "SchGA.h"

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

void SchGADemo() {
    // Task and flight descriptor
    TASK_PARAMETER taskParameter[2][3];
    TASK_QUEUE taskQueue[2];
    TASK task;

    taskParameter[0][0].x = 1;
    taskParameter[0][0].y = 1;
    taskParameter[0][0].deg = 0;
    taskParameter[0][1].x = 1;
    taskParameter[0][1].y = 1;
    taskParameter[0][1].deg = 0;
    taskParameter[0][2].x = 1;
    taskParameter[0][2].y = 1;
    taskParameter[0][2].deg = 0;

    taskParameter[1][0].x = 1000;
    taskParameter[1][0].y = 1000;
    taskParameter[1][0].deg = 0;
    taskParameter[1][1].x = 1000;
    taskParameter[1][1].y = 1000;
    taskParameter[1][1].deg = 0;
    taskParameter[1][2].x = 1000;
    taskParameter[1][2].y = 1000;
    taskParameter[1][2].deg = 0;

    taskQueue[0].num = 3;
    taskQueue[0].tasks = taskParameter[0];
    taskQueue[1].num = 3;
    taskQueue[1].tasks = taskParameter[1];

    task.queueNum = 2;
    task.totalNum = 6;
    task.taskQueue = taskQueue;
    
    SINGLE_FLIGHT_STATE singleFlightState[1];
    FLIGHT_STATE flightState;
    singleFlightState[0].x = 1010;
    singleFlightState[0].y= 1010;
    singleFlightState[0].deg= 0;
    flightState.num = 1;
    flightState.flightState = singleFlightState;

    // Task scheduler
    SchGA schGA(100, &flightState, &task, 0.01, 0.999, 1-0.0001);
    uint32_t aaa;
    double_t bestFitness;
    schGA.evaluate(100, aaa, bestFitness);

    const uint32_t* bestGene = schGA.getBestGene();
    uint32_t geneLength = schGA.getGeneLength();
    std::cout << "Best Fitness: " << bestFitness << std::endl;
    std::cout << "Best Gene: ";
    for(uint32_t i = 0; i < geneLength; i++) {
        std::cout << bestGene[i] << ' ';
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "Start" << std::endl;
    SchGADemo();
    std::cout << "Done" << std::endl;

    return 0;
}
