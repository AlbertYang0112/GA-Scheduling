#ifndef GAHW_SCHGA_H
#define GAHW_SCHGA_H

#include "GA.h"
#include <random>
#include <vector>
#include <fstream>
#include "Descriptor.h"
#include "SearchEngine.h"
extern "C" {
#include "dubins.h"
};


class SchGA: public GA {
public:
    friend class SearchEngine;
    void evaluate(uint32_t iterations,
                  uint32_t &bestGene, double_t &bestFitness) override;
    SchGA(uint32_t population, FLIGHT_STATE* flights, TASK* taskTable, 
        double_t rho, double_t crossRate, double_t mutationRate, 
        const char* recorderName
        );
    ~SchGA();
    const uint32_t *getBestGene();
    uint32_t getGeneLength();
    double_t fitnessAverage();
    double_t fitnessVar();
    double_t distantToBestAverage();
    double_t distantToBestVar();

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

    /*
     * Resolve the sequence conflict
     * Param:
     *  timeStamp: the task execution time, len = task->totalNum；
     *  flightNo: the flight number which executed the according task;
     * Return: The task time, if the conflict is unresolvable, return DBL_MAX
     */
    double_t _timeCompute(double_t* timeStamp, uint32_t* flightNo, uint32_t max_step=100);
    uint32_t _distance(uint32_t* pGeneA, uint32_t* pGeneB);
    uint32_t _selectByDistant(uint32_t** geneList, uint32_t minDist, uint32_t num);

    // Gene Storage: Ping-pong buffer
    uint32_t *_gene;
    uint32_t *_nextGene;

    // Fitness Buffer
    double_t *_fitness;

    // Path Buffer
    DubinsPath* _path;
    uint32_t _numPath;

    // Parameters
    TASK *_taskTable;
    TASK_PARAMETER *_taskParamTable;
    FLIGHT_STATE *_initialFlightState;
    uint32_t _population;
    uint32_t _geneLength;
    uint32_t _crossRate;
    uint32_t _mutationRate;
    double_t _rho;
#ifdef ADVANCED_FEATURE
    static const uint32_t PRESERVED_SLOT = 6;
    static const uint32_t SEARCH_ENGINE_NUM = 5;
    static const uint32_t TABOO_LIST_LEN = 5;
#else
    static const uint32_t PRESERVED_SLOT = 6;
#endif

    // Signals
    bool _bestFitnessUpdated;

    // Statistic Data
    uint32_t* _bestGene;
    double_t _bestFitness;
    uint32_t _feasibleGeneCnt;
    uint32_t _saveCnt;

    // Submodules
    std::mt19937 _rng;
#ifdef ADVANCED_FEATURE
    std::vector<SearchEngine> _searchEngines;
#endif
    std::ofstream _recorder;

    // Preallocated Buffer for time computation
    uint32_t* __timeSeparator;
    bool* __timeCircleDetect;
    uint32_t* __timeStampSeq;
};


inline TASK_PARAMETER* SchGA::_visitTask(uint32_t task) {
    if(task > _taskTable->totalNum) {
        return nullptr;
    }
    return _taskParamTable + task;
}

inline const uint32_t* SchGA::getBestGene() {
    return const_cast<const uint32_t*>(_bestGene);
}

inline uint32_t SchGA::getGeneLength() {
    return _geneLength;
}

#endif //GAHW_SCHGA_H
