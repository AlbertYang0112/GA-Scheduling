#ifndef GAHW_SCHGA_H
#define GAHW_SCHGA_H

#include "GA.h"
#include <random>
extern "C" {
#include "dubins.h"
};

class SchGA: public GA {
public:
    typedef struct tTASK_PARAMETER {
        double_t deg;
        double_t x;
        double_t y;
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
        SINGLE_FLIGHT_STATE *flightState;
    } FLIGHT_STATE;

    void evaluate(uint32_t iterations,
                  uint32_t &bestGene, double_t &bestFitness) override;
    SchGA(uint32_t population, FLIGHT_STATE* flights, TASK* taskTable, double_t rho, double_t crossRate, double_t mutationRate);
    ~SchGA();
    void mutationTest();
    const uint32_t *getBestGene();
    uint32_t getGeneLength();

private:
    TASK_PARAMETER* _visitTask(uint32_t task);

    /*
     * Generate the initial gene sequences.
     * Note: Make sure the _gene points to a valid memory
     *          with size (population * (taskNum + flightNum))
     */
    void _generateInitGene() override;
    void _fitnessCal() override;
    void _cross(uint32_t parentA, uint32_t parentB, uint32_t &childA, uint32_t &childB) override;
    uint32_t _mutation(uint32_t child) override;
    void _selectParents(uint32_t* parentsNo, uint32_t num);

    TASK *_taskTable;
    FLIGHT_STATE *_initialFlightState;
    uint32_t _population;
    uint32_t *_gene;
    uint32_t *_nextGene;
    uint32_t _geneLength;
    uint32_t* _bestGene;
    double_t _bestFitness;
    uint32_t _feasibleGeneCnt;
    uint32_t _crossRate;
    uint32_t _mutationRate;
    double_t *_fitness;
    double_t _rho;
    std::mt19937 _rng;

    DubinsPath* _path;
    uint32_t _numPath;
};

inline SchGA::TASK_PARAMETER* SchGA::_visitTask(uint32_t task) {
    if(task > _taskTable->totalNum) {
        return nullptr;
    }
    TASK_QUEUE *pQueue = _taskTable->taskQueue;
    while(task >= pQueue->num) {
        task -= pQueue->num;
        pQueue++;
    }
    return (pQueue->tasks + task);
}

inline const uint32_t* SchGA::getBestGene() {
    return const_cast<const uint32_t*>(_bestGene);
}

inline uint32_t SchGA::getGeneLength() {
    return _geneLength;
}

#endif //GAHW_SCHGA_H
