#ifndef GAHW_SCHGA_H
#define GAHW_SCHGA_H

#include "GA.h"
#include <random>
#include "Descriptor.h"
extern "C" {
#include "dubins.h"
};

class SchGA: public GA {
public:
    void evaluate(uint32_t iterations,
                  uint32_t &bestGene, double_t &bestFitness) override;
    SchGA(uint32_t population, FLIGHT_STATE* flights, TASK* taskTable, double_t rho, double_t crossRate, double_t mutationRate);
    ~SchGA();
    const uint32_t *getBestGene();
    uint32_t getGeneLength();
    double_t fitnessAverage();
    double_t fitnessVar();

private:
    TASK_PARAMETER* _visitTask(uint32_t task);

    /*
     * Generate the initial gene sequences.
     * Note: Make sure the _gene points to a valid memory
     *          with size (population * (taskNum + flightNum))
     */
    void _generateInitGene() override;
    void _fitnessCal() override;
    double_t _singleFitnessCal(uint32_t* pGene);
    void _cross(uint32_t parentA, uint32_t parentB, uint32_t &childA, uint32_t &childB) override;
    uint32_t _mutation(uint32_t child) override;
    void _selectParents(uint32_t* parentsNo, uint32_t num);
    double_t _timeCompute(double_t* timeStamp, uint32_t* filghtNo);

    TASK *_taskTable;
    TASK_PARAMETER *_taskParamTable;
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

inline TASK_PARAMETER* SchGA::_visitTask(uint32_t task) {
    if(task > _taskTable->totalNum) {
        return nullptr;
    }
    //TASK_QUEUE *pQueue = _taskTable->taskQueue;
    //while(task >= pQueue->num) {
    //    task -= pQueue->num;
    //    pQueue++;
    //}
    //return (pQueue->tasks + task);
    return _taskParamTable + task;
}

inline const uint32_t* SchGA::getBestGene() {
    return const_cast<const uint32_t*>(_bestGene);
}

inline uint32_t SchGA::getGeneLength() {
    return _geneLength;
}

#endif //GAHW_SCHGA_H
