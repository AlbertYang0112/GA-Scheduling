#ifndef GAHW_SCHGA_H
#define GAHW_SCHGA_H

#include "GA.h"
#include "dubins.h"

class SchGA: public GA {
public:
    typedef struct tTASK_PARAMETER {
        double_t deg;
        double_t x;
        double_t y;
        uint32_t taskID;
    } TASK_PARAMETER;
    typedef struct tTASK_QUEUE {
        uint32_t num;
        TASK_PARAMETER *tasks;
    } TASK_QUEUE;
    typedef struct tTASK {
        uint32_t queueNum;
        uint32_t totalNum;
        TASK_QUEUE *taskQueue;
    } TASK;
    typedef struct tSINGLE_FLIGHT_STATE {
        double_t deg;
        double_t x;
        double_t y;
    } SINGLE_FLIGHT_STATE;
    typedef struct tFLIGHT_STATE {
        uint32_t num;
        SINGLE_FLIGHT_STATE *filghtState;
    } FLIGHT_STATE;

    void evaluate(uint32_t iterations,
                  uint32_t &bestGene, double &bestFitness) override;
    SchGA(uint32_t population, FLIGHT_STATE* flights, TASK* taskTable, bool copy=false);

private:
    void _generateInitGene() override;
    void _fitnessCal() override;
    void _cross(uint32_t parentA, uint32_t parentB, uint32_t &childA, uint32_t &childB) override;
    uint32_t _mutation(uint32_t child) override;

    TASK *_taskTable;
    FLIGHT_STATE *_initialFlightState;
    uint32_t _population;
    uint32_t * _gene;
    DubinsPath* _path;
    uint32_t _numPath;
};
#endif //GAHW_SCHGA_H
