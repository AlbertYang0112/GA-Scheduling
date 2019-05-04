#include "SchGA.h"
#include <values.h>


void SchGA::_generateInitGene() {
    uint32_t *pGene = _gene;

    // Generate the basic gene pattern
    _geneLength = _taskTable->totalNum + _initialFlightState->num;
    uint32_t baseGene[_geneLength];
    for(uint32_t i = 0; i < _taskTable->totalNum; i++) {
        baseGene[i] = i;
    }
    std::fill(baseGene + _taskTable->totalNum, baseGene + _geneLength, _taskTable->totalNum);

    // Shuffler
    std::random_device rd;
    std::mt19937 g(rd());

    for(uint32_t gene = 0; gene < _population; gene++) {
        std::shuffle(baseGene, baseGene + _geneLength, g);
        std::copy(baseGene, baseGene + _geneLength, pGene);
        pGene += _geneLength;
    }
}

void SchGA::_fitnessCal() {
    uint32_t *pGene = _gene;

    TASK_PARAMETER *pTask;

    FLIGHT_STATE flightState;
    flightState.num = _initialFlightState->num;
    flightState.flightState = new SINGLE_FLIGHT_STATE[flightState.num];

    DubinsPath path;
    double_t startPoint[3];
    double_t endPoint[3];

    double_t taskTime[_taskTable->totalNum];
    double_t totalTaskTime;
    double_t flyTime;

    double_t fitness[_population];
    double_t maxTime = 0;
    double_t minTime = MAXDOUBLE;


    for(uint32_t gene = 0; gene < _population; gene++) {

        // Initialize the flight state
        std::copy(_initialFlightState->flightState,
                _initialFlightState->flightState + _initialFlightState->num,
                flightState.flightState
                );

        std::fill(taskTime, taskTime + _taskTable->totalNum, 0);

        uint32_t *pTempGene = pGene;
        for(uint32_t flight = 0; flight < flightState.num; flight++) {
            totalTaskTime = 0;
            while(*pTempGene != _taskTable->totalNum) {
                pTask = _visitTask(*pTempGene);

                // Generate the path
                // Todo: Optimize the memory access.
                startPoint[0] = flightState.flightState->x;
                startPoint[1] = flightState.flightState->y;
                startPoint[2] = flightState.flightState->deg;
                endPoint[0] = pTask->x;
                endPoint[1] = pTask->y;
                endPoint[2] = pTask->deg;
                dubins_shortest_path(&path, startPoint, endPoint, _rho);

                // Todo: Collision detection

                flyTime = dubins_path_length(&path);
                totalTaskTime += flyTime;
                taskTime[*pTempGene] = totalTaskTime;
                pTempGene++;
            }
        }

        // Sequence Constrain
        double_t* pTaskTime = taskTime;
        TASK_QUEUE* pTaskQueue = _taskTable->taskQueue;
        bool isFeasibleSolution = true;
        for(uint32_t queue = 0; queue < _taskTable->queueNum & isFeasibleSolution; queue++) {
            totalTaskTime = 0;
            for(uint32_t task = 0; task < pTaskQueue->num; task++) {
                if(totalTaskTime > *pTaskTime) {
                    // Constrain Violation Occurs
                    isFeasibleSolution = false;
                    break;
                    // Todo: Add the soft-margin method
                }
                totalTaskTime = *pTaskTime++;
            }
            pTaskQueue++;
        }
        fitness[gene] = 0;
        if(isFeasibleSolution) {
            // Sum to get the total fitness
            pTaskTime = taskTime;
            for(uint32_t task = 0; task < _taskTable->totalNum; task++) {
                fitness[gene] += *pTaskTime++;
            }
            if(fitness[gene] > maxTime) {
                maxTime = fitness[gene];
            } else if(fitness[gene] < minTime) {
                minTime = fitness[gene];
                _bestGene = *pGene;
            }
        } else {
            // Set the fitness of the infeasible solution to -1.
            // We do not set the fitness directly to 0 in case of
            // there exist a solution whose fitness is exactly 0.
            fitness[gene] = -1;
        }
        pGene += _geneLength;
    }

    // Reverse the fitness, the higher the better.
    // Set the fitness of the infeasible solution to 0
    for(uint32_t gene = 0; gene < _population; gene++) {
        if(fitness[gene] >= 0) {
            fitness[gene] = maxTime - fitness[gene];
        } else {
            fitness[gene] = 0;
        }
    }
    delete [] flightState.flightState;
}


