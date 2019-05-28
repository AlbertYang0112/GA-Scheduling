#include <iostream>
#include <iomanip>
#include <exception>
#include <cfloat>
#include "Descriptor.h"
#include "NaiveGA.h"
#include "SchGA.h"
#include "TestGenerator.h"
#include "debug.h"

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

void SchGADemoStatic() {
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
    SchGA schGA(6000, &flightState, &task, 0.01, 0.5, 0.3);
    uint32_t aaa;
    double_t bestFitness;
    schGA.evaluate(5000, aaa, bestFitness);

    const uint32_t* bestGene = schGA.getBestGene();
    uint32_t geneLength = schGA.getGeneLength();
    std::cout << "Best Fitness: " << bestFitness << std::endl;
    std::cout << "Best Gene: ";
    for(uint32_t i = 0; i < geneLength; i++) {
        std::cout << bestGene[i] << ' ';
    }
    std::cout << std::endl;
}

void SchGADemoRandom() {
    // Configurations
    // Task and flight
    const uint32_t MAX_FLIGHT_NUM = 3;
    const uint32_t MIN_FLIGHT_NUM = 3;
    const uint32_t MAX_TASK_QUEUE_NUM = 100;
    const uint32_t MIN_TASK_QUEUE_NUM = 100;
    const uint32_t MAX_TASK_QUEUE_LEN = 3;
    const uint32_t MIN_TASK_QUEUE_LEN = 3;
    const double_t RHO = 5;
    // GA
    const uint32_t POPULATION = 2000;
    const uint32_t ITERATIONS = 10000;
    const double_t MUTATION_RATE = 0.3;
    const double_t CROSS_RATE = 0.7;
    TASK task;
    FLIGHT_STATE flightState;

    // Random Number Generator
    std::random_device rd;
    std::mt19937 rng;
    rng.seed(rd());

    // Generate the random task and flight configuration
    uint32_t queueNum;
    uint32_t* queueLen;
    uint32_t taskSum = 0;
    uint32_t flightNum;
    char acceptRecv;

    while (true) {
        queueNum = rng() % (MAX_TASK_QUEUE_NUM - MIN_TASK_QUEUE_NUM + 1) + MIN_TASK_QUEUE_NUM;
        queueLen = new uint32_t[queueNum];
        for(uint32_t queueNo = 0; queueNo < queueNum; queueNo++) {
            queueLen[queueNo] = rng() % (MAX_TASK_QUEUE_LEN - MIN_TASK_QUEUE_LEN + 1) + MIN_TASK_QUEUE_LEN;
            taskSum += queueLen[queueNo];
        }
        flightNum = rng() % (MAX_FLIGHT_NUM - MIN_FLIGHT_NUM + 1) + MIN_FLIGHT_NUM;

        printf(
            "Test configuration:\n" \
            "\tFlight Num: %d\n" \
            "\tTask Queue Num: %d\tTotal Task Num: %d\n" \
            "\tTask Queue Len: ",
            flightNum, queueNum, taskSum
        );
        for(uint32_t queueNo = 0; queueNo < queueNum; queueNo++) {
            printf("%d ", queueLen[queueNo]);
        }
        break;
        printf("\n\tAccept?(y/n) ");
        scanf("%c", &acceptRecv);
        if(acceptRecv == 'Y' | acceptRecv == 'y') break;
        delete [] queueLen;
    }
    printf(
        "Task Scheduler Configuration:\n"\
        "\tGA Population: %d\t Iteration: %d\n" \
        "\tMutation Rate: %.2f\n" \
        "\tCross Rate: %.2f\n",
        POPULATION, ITERATIONS,
        static_cast<float>(MUTATION_RATE),
        static_cast<float>(CROSS_RATE)
    );

    createFlight(&flightState, flightNum);
    createTask(&task, queueNum, queueLen);

    randomFlightState(&flightState);
    randomTask(&task);

    /*
    printf("Flight Data:\n");
    for(uint32_t i = 0; i < flightNum; i++) {
        printf("%f %f %f\n", flightState.flightState[i].x, flightState.flightState[i].y, flightState.flightState[i].deg);
    }
    printf("Task Data:\n");
    for(uint32_t i = 0; i < queueNum; i++) {
        printf("%f %f %f\n", task.taskQueue[i].tasks[0].x, task.taskQueue[i].tasks[0].y, task.taskQueue[i].tasks[0].deg);
    }
    */

    delete [] queueLen;

    // The Task Scheduler
    SchGA schGA(POPULATION, &flightState, &task, RHO, CROSS_RATE, MUTATION_RATE);
    uint32_t fillBlank;
    double_t bestFitness;
    // Run
    // schGA.evaluate(ITERATIONS, fillBlank, bestFitness);
    try {
        schGA.evaluate(ITERATIONS, fillBlank, bestFitness);
    } catch(std::exception& e) {
        printf("Exception: %s\n", e.what());
        destructTask(&task);
        destructFlight(&flightState);
        return;
    }

    // Collect the Result
    const uint32_t* bestGene = schGA.getBestGene();
    double_t avg = schGA.fitnessAverage();
    if(avg != DBL_MAX) {
        uint32_t geneLength = schGA.getGeneLength();
        printf("Best Fitness: %f\nBest Task Sequence: ", static_cast<float>(bestFitness));
        for (uint32_t i = 0; i < geneLength; i++) {
            if (bestGene[i] < taskSum) {
                printf("%d ", bestGene[i]);
            } else {
                printf("S ");
            }
        }
        printf("S\n");
    } else {
        printf("No Feasible Solution\n");
    }
    destructTask(&task);
    destructFlight(&flightState);
}

int main() {
    std::cout << "Start" << std::endl;
    SchGADemoRandom();
    //SchGADemoStatic();
    std::cout << "Done" << std::endl;

    return 0;
}
