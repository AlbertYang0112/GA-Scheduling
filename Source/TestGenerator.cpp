#include <random>
#include "TestGenerator.h"
#include "Descriptor.h"

void randomFlightState(
    FLIGHT_STATE* flightState, 
    double_t xUpperBound, 
    double_t yUpperBound, 
    double_t degUpperBound,
    double_t xLowerBound, 
    double_t yLowerBound,
    double_t degLowerBound) {

    // Random Number Generator
    std::random_device rd;
    std::mt19937 rng;
    rng.seed(rd());
    std::uniform_real_distribution<double_t> xRng(xLowerBound, xUpperBound);
    std::uniform_real_distribution<double_t> yRng(yLowerBound, yUpperBound);
    std::uniform_real_distribution<double_t> degRng(degLowerBound, degUpperBound);

    for(uint32_t i = 0; i < flightState->num; i++) {
        flightState->flightState[i].x = xRng(rng);
        flightState->flightState[i].y = yRng(rng);
        flightState->flightState[i].deg = degRng(rng);
    }
}

void randomTask(
    TASK* task, 
    double_t xUpperBound, 
    double_t yUpperBound, 
    double_t degUpperBound,
    double_t xLowerBound, 
    double_t yLowerBound,
    double_t degLowerBound) {

    // Random Number Generator
    std::random_device rd;
    std::mt19937 rng;
    rng.seed(rd());
    std::uniform_real_distribution<double_t> xRng(xLowerBound, xUpperBound);
    std::uniform_real_distribution<double_t> yRng(yLowerBound, yUpperBound);
    std::uniform_real_distribution<double_t> degRng(degLowerBound, degUpperBound);
    for(uint32_t queueNo = 0; queueNo < task->queueNum; queueNo++) {
        TASK_QUEUE *pQueue = task->taskQueue + queueNo;
        for(uint32_t taskNo = 0; taskNo < pQueue->num; taskNo++) {
            pQueue->tasks[taskNo].x = xRng(rng);
            pQueue->tasks[taskNo].y = yRng(rng);
            pQueue->tasks[taskNo].deg = degRng(rng);
        }
    }
}