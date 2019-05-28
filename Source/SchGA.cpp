#include "SchGA.h"
#include <cfloat>
#include <cassert>
#include "Descriptor.h"
#include "debug.h"


SchGA::SchGA(
        uint32_t population,
        FLIGHT_STATE *flights,
        TASK *taskTable,
        double_t rho,
        double_t crossRate,
        double_t mutationRate
        ):_rho(rho), _population(population) {
    assert(flights != nullptr);
    assert(taskTable != nullptr);

    // Copy the task table
    _taskTable = new TASK;
    _taskTable->queueNum = taskTable->queueNum;
    _taskTable->totalNum = taskTable->totalNum;
    _taskTable->taskQueue = new TASK_QUEUE[_taskTable->queueNum];
    _taskParamTable = new TASK_PARAMETER[_taskTable->totalNum];
    TASK_QUEUE* pSrcQueue = taskTable->taskQueue;
    TASK_QUEUE* pDstQueue = _taskTable->taskQueue;

    _taskTable->totalNum = 0;
    for(uint32_t queue = 0; queue < _taskTable->queueNum; queue++) {
        pDstQueue->num = pSrcQueue->num;

        pDstQueue->tasks = _taskParamTable + _taskTable->totalNum;
        _taskTable->totalNum += pDstQueue->num;
        std::copy(
                pSrcQueue->tasks,
                pSrcQueue->tasks + pDstQueue->num,
                pDstQueue->tasks);
        pSrcQueue++;
        pDstQueue++;
    }

    assert(_taskTable->totalNum == taskTable->totalNum);


    // Copy the initial flight state
    _initialFlightState = new FLIGHT_STATE;
    _initialFlightState->num = flights->num;
    _initialFlightState->flightState = new SINGLE_FLIGHT_STATE[_initialFlightState->num];
    std::copy(
            flights->flightState,
            flights->flightState + _initialFlightState->num,
            _initialFlightState->flightState);

    // Allocate the memory for the gene
    _gene = new uint32_t[
            population * (_taskTable->totalNum + _initialFlightState->num)
            ];
    _nextGene = new uint32_t[
            population * (_taskTable->totalNum + _initialFlightState->num)
            ];

    _bestGene = nullptr;
    
    // Calculate the length of the gene
    _geneLength = _taskTable->totalNum + _initialFlightState->num - 1;
    
    _fitness = new double_t[_population];

    // Allocate the space for the search engine
    _searchGene = new uint32_t[_geneLength];

    // Initialize the random generator
    std::random_device rd;
    _rng.seed(rd());

    // Scale the cross rate and the mutation rate
    assert(crossRate > 0 && crossRate < 1);
    assert(mutationRate > 0 && mutationRate < 1);
    _crossRate = crossRate * _rng.max();
    _mutationRate = mutationRate * _rng.max();
}

void SchGA::_generateInitGene() {
    uint32_t *pGene = _gene;

    // Generate the basic gene pattern
    uint32_t baseGene[_geneLength];
    for(uint32_t i = 0; i < _geneLength; i++) {
        baseGene[i] = i;
    }

    std::fill(_gene, _gene + _population * _geneLength, 1234);
    std::fill(_nextGene, _nextGene + _population * _geneLength, 4321);

    // Shuffle the basic gene to generate the population
    for(uint32_t gene = 0; gene < _population; gene++) {
        std::shuffle(baseGene, baseGene + _geneLength, _rng);
        std::copy(baseGene, baseGene + _geneLength, pGene);
        pGene += _geneLength;
    }
    uint32_t i = 0;
    pGene = _gene;
    for(uint32_t queueNo = 0; queueNo < _taskTable->queueNum; queueNo++) {
        for(uint32_t j = 0; j <_taskTable->taskQueue[queueNo].num; j++) {
            *pGene++ = j + i;
        }
        i += _taskTable->taskQueue[queueNo].num;
        if(queueNo < _initialFlightState->num - 1)
            *pGene++ = _taskTable->totalNum + queueNo;
    }
    DEBUG_BRIEF("%d\n", _geneLength);
    for(uint32_t i = 0; i < _geneLength; i++) {
        DEBUG_BRIEF("%d ", _gene[i]);
    }
    DEBUG_BRIEF("\n");
}

double_t SchGA::_singleFitnessCal(uint32_t* pGene) {
    TASK_PARAMETER *pTask;

    DubinsPath path;
    double startPoint[3];
    double endPoint[3];

    double_t taskTime[_taskTable->totalNum];
    double_t totalTaskTime;
    double_t flyTime;
    double_t tempMaxTime = 0;

    std::fill(taskTime, taskTime + _taskTable->totalNum, 0);

    uint32_t *pTempGene = pGene;
    for(uint32_t flight = 0; flight < _initialFlightState->num; flight++) {
        totalTaskTime = 0;
        startPoint[0] = _initialFlightState->flightState[flight].x;
        startPoint[1] = _initialFlightState->flightState[flight].y;
        startPoint[2] = _initialFlightState->flightState[flight].deg;
        while(*pTempGene < _taskTable->totalNum && pTempGene != pGene + _geneLength) {
            pTask = _visitTask(*pTempGene);

            // Generate the path
            endPoint[0] = pTask->x;
            endPoint[1] = pTask->y;
            endPoint[2] = pTask->deg;
            dubins_shortest_path(&path, startPoint, endPoint, _rho);

            startPoint[0] = endPoint[0];
            startPoint[1] = endPoint[1];
            startPoint[2] = endPoint[2];

            // Todo: Collision detection

            flyTime = dubins_path_length(&path);
            totalTaskTime += flyTime;
            taskTime[*pTempGene] = totalTaskTime;
            if(totalTaskTime > tempMaxTime) {
                tempMaxTime = totalTaskTime;
            }
            pTempGene++;
        }
        pTempGene++;
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
    if(isFeasibleSolution) {
        return tempMaxTime;
    } else {
        return DBL_MAX;
    }
}

void SchGA::_fitnessCal() {
    double_t maxTime = 0;
    double_t prevBestFitness = _bestFitness;
    auto minTime = DBL_MAX;

//#pragma omp parallel for num_threads(8)
    for(uint32_t gene = 0; gene < _population; gene++) {
        _fitness[gene] = _singleFitnessCal(_gene + gene * _geneLength);
    }

    _feasibleGeneCnt = 0;
    for(uint32_t gene = 0; gene < _population; gene++) {
        if(_fitness[gene] != DBL_MAX) {
            _feasibleGeneCnt++;
            if(_fitness[gene] > maxTime) {
                maxTime = _fitness[gene];
            }
            if(_fitness[gene] < minTime) {
                minTime = _fitness[gene];
                _bestGene = _gene + gene * _geneLength;
                _bestFitness = _fitness[gene];
            }
        }
    }
    _bestFitnessUpdated = _bestFitness != prevBestFitness;
}

void SchGA::_cross(uint32_t parentA, uint32_t parentB,
        uint32_t &childA, uint32_t &childB) {
    uint32_t *pGeneA = _gene + parentA * _geneLength;
    uint32_t *pGeneB = _gene + parentB * _geneLength;

    uint32_t *pNextGeneA = _nextGene + childA * _geneLength;
    uint32_t *pNextGeneB = _nextGene + childB * _geneLength;


    uint32_t sectionStart, sectionLen;
    sectionStart = static_cast<uint32_t>(_rng()) % _geneLength;
    // Todo: Define the length limit of the cross section as a hyper-parameter
    sectionLen = static_cast<uint32_t>(_rng()) % _geneLength + 1;

    uint32_t crossList[_geneLength];
    for(uint32_t i = 0; i < _geneLength; i++) {
        crossList[i] = i;
    }


    // Generate the swap table
    uint32_t pos = sectionStart;
    for(uint32_t i = 0; i < sectionLen; i++) {
        if( crossList[pGeneA[pos]] == pGeneA[pos] &&
            crossList[pGeneB[pos]] == pGeneB[pos] ) {
            // No swap collision
            crossList[pGeneA[pos]] = pGeneB[pos];
            crossList[pGeneB[pos]] = pGeneA[pos];
        }
        if(pos == _geneLength - 1) {
            pos = 0;
        } else {
            pos++;
        }
    }
    for(uint32_t i = 0; i < _geneLength; i++) {
        *pNextGeneA++ = crossList[*pGeneA++];
        *pNextGeneB++ = crossList[*pGeneB++];
    }

}
uint32_t SchGA::_mutation(uint32_t child) {
    uint32_t *pGene = _nextGene + _geneLength * child;
    uint32_t sectionStart, sectionLen;
    uint32_t swapSpace[_geneLength];
    sectionStart = static_cast<uint32_t>(_rng()) % _geneLength;
    // Todo: Define the length limit of the mutation section as a hyper-parameter
    sectionLen = static_cast<uint32_t>(_rng()) % _geneLength + 1;

    if(sectionStart + sectionLen <= _geneLength) {
        std::shuffle(pGene + sectionStart, pGene + sectionStart + sectionLen, _rng);
    } else {
        // Concatenate the segment into the swap space.
        std::copy(pGene + sectionStart, pGene + _geneLength, swapSpace);
        std::copy(pGene, pGene + sectionLen - (_geneLength - sectionStart),
                swapSpace + _geneLength - sectionStart);

        // Shuffle
        std::shuffle(swapSpace, swapSpace + sectionLen, _rng);

        // Write back
        std::copy(swapSpace, swapSpace + _geneLength - sectionStart, pGene + sectionStart);
        std::copy(swapSpace + _geneLength - sectionStart, swapSpace + sectionLen, pGene);
    }
    return 0;
}

void SchGA::_selectParents(uint32_t* parentsNo, uint32_t num) {
    uint32_t selA, selB;
    for(uint32_t i = 0; i < num; i++) {
        selA = _rng() % _population;
        selB = _rng() % _population;
        if(_rng() < _rng.max() * 0.9) {
            if (_fitness[selA] >= _fitness[selB]) {
                parentsNo[i] = selB;
            } else {
                parentsNo[i] = selA;
            }
        } else {
            if (_fitness[selA] >= _fitness[selB]) {
                parentsNo[i] = selA;
            } else {
                parentsNo[i] = selB;
            }
        }
    }
}

void SchGA::evaluate(
        uint32_t iterations,
        uint32_t &bestGene,
        double_t &bestFitness) {
    // Generate the initial gene
    _generateInitGene();

    // Initiate the search engine
    _search(true);

    bool searchUpdated;
    uint32_t parentsNo[2];
    for(uint32_t iter_cnt = 0; iter_cnt < iterations; iter_cnt++) {
        _fitnessCal();
        if(_feasibleGeneCnt != 0) {
            // If we found a feasible solution

            searchUpdated = _search(false);
            if(searchUpdated) {
                DEBUG_BRIEF("Found a better gene by search: %f -> %f\n", _bestFitness, _searchBestFitness);
            }
            if(iter_cnt % 10 == 0) {
                DEBUG_BRIEF("> GA %f Search %f\n", _bestFitness, _searchBestFitness);
            }
            if(_bestFitnessUpdated && searchUpdated) {
                // Put the gene with better fitness to the first place
                // Put the other gene to the secondary place
                if(_bestFitness < _searchBestFitness) {
                    std::copy(_bestGene, _bestGene + _geneLength, _nextGene);
                    std::copy(_searchGene, _searchGene + _geneLength, _nextGene + _geneLength);
                    _fitness[0] = _bestFitness;
                    _fitness[1] = _searchBestFitness;
                    // Update the search engine
                    _search(true);
                } else {
                    std::copy(_searchGene, _searchGene + _geneLength, _nextGene);
                    std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
                    _fitness[0] = _searchBestFitness;
                    _fitness[1] = _bestFitness;
                }
            } else if(searchUpdated) {
                std::copy(_searchGene, _searchGene + _geneLength, _nextGene);
                std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
                _fitness[0] = _searchBestFitness;
                _fitness[1] = _bestFitness;
            } else {
                // Simply Preserve the best gene
                std::copy(_bestGene, _bestGene + _geneLength, _nextGene);
                std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
                _fitness[0] = _bestFitness;
                _fitness[1] = _bestFitness;
            }
            if (_rng() < _mutationRate) {
                _mutation(1);
            }
        }

        for(uint32_t cross_cnt = (_feasibleGeneCnt == 0 ? 0 : 2); cross_cnt < _population; cross_cnt += 2) {
            _selectParents(parentsNo, 2);
            uint32_t cross_tmp = cross_cnt + 1;
            if(_fitness[parentsNo[0]] == DBL_MAX || _fitness[parentsNo[1]] == DBL_MAX) {
                std::copy(_gene + parentsNo[0] * _geneLength,
                          _gene + parentsNo[0] * _geneLength + _geneLength,
                          _nextGene + cross_cnt * _geneLength);
                std::copy(_gene + parentsNo[1] * _geneLength,
                          _gene + parentsNo[1] * _geneLength + _geneLength,
                          _nextGene + cross_cnt * _geneLength + _geneLength);
                if(_fitness[parentsNo[0]] == DBL_MAX || parentsNo[0] < 2) {
                    _mutation(cross_cnt);
                    _mutation(cross_cnt);
                    _mutation(cross_cnt);
                }
                if(_fitness[parentsNo[1]] == DBL_MAX || parentsNo[1] < 2) {
                    _mutation(cross_tmp);
                    _mutation(cross_tmp);
                    _mutation(cross_tmp);
                }
            } else {
                if (_rng() < _crossRate) {
                    _cross(parentsNo[0], parentsNo[1], cross_cnt, cross_tmp);
                } else {
                    std::copy(_gene + parentsNo[0] * _geneLength,
                              _gene + parentsNo[0] * _geneLength + _geneLength,
                              _nextGene + cross_cnt * _geneLength);
                    std::copy(_gene + parentsNo[1] * _geneLength,
                              _gene + parentsNo[1] * _geneLength + _geneLength,
                              _nextGene + cross_cnt * _geneLength + _geneLength);
                }
                if (_rng() < _mutationRate) {
                    _mutation(cross_cnt);
                }
                if (_rng() < _mutationRate) {
                    _mutation(cross_tmp);
                }
            }
        }

        if(iter_cnt % 10 == 0) {
            double_t avg = fitnessAverage();
            double_t var = fitnessVar();
            DEBUG_BRIEF("Iteration: %d\n", iter_cnt);
            DEBUG_BRIEF("Feasible Gene: %d\n", _feasibleGeneCnt);
            if(avg != DBL_MAX) {
                DEBUG_BRIEF("Fitness Average: %f\n", avg);
                DEBUG_BRIEF("Fitness Variance: %E\n", var);
            }
            if(_feasibleGeneCnt != 0) {
                DEBUG_BRIEF("Best Fitness %f\n", static_cast<double>(_bestFitness));
                DEBUG_BRIEF("Search Best Fitness %f\n", _searchBestFitness);
                //for(uint32_t i = 0; i < _geneLength; i++) {
                //    if(_bestGene[i] < _taskTable->totalNum) {
                //        DEBUG_BRIEF("%d ", _bestGene[i]);
                //    } else {
                //        DEBUG_BRIEF("S ");
                //    }
                //}
                //DEBUG_BRIEF("S\n");
            }
        }
        std::swap(_gene, _nextGene);
    }
    bestFitness = _bestFitness;
    bestGene = 0;
}

double_t SchGA::fitnessAverage() {
    double_t fitnessSum = 0;
    uint32_t cnt = 0;
    for(uint32_t i = 0; i < _population; i++) {
        if(_fitness[i] != DBL_MAX) {
            fitnessSum += _fitness[i];
            cnt++;
        }
    }
    if(cnt == 0) return DBL_MAX;
    return fitnessSum / cnt;
}

double_t SchGA::fitnessVar() {
    double_t fitnessAvg = fitnessAverage();

    if(fitnessAvg == DBL_MAX) return DBL_MAX;

    double_t fitnessVarSum = 0;
    double_t term;
    uint32_t cnt = 0;
    for(uint32_t i = 0; i < _population; i++) {
        if(_fitness[i] != DBL_MAX) {
            term = _fitness[i] - fitnessAvg;
            fitnessVarSum += term * term;
            cnt++;
        }
    }
    return sqrt(fitnessVarSum / cnt);
}

bool SchGA::_search(bool update) {
    const uint32_t SAMPLE_POINTS = 10;
    const uint32_t TABO_LIST_LEN = 5;
    static uint32_t taboListH[TABO_LIST_LEN];
    static uint32_t taboListL[TABO_LIST_LEN];
    static uint32_t taboListPos = 0;

    double_t fitness[SAMPLE_POINTS];
    double_t minFit = DBL_MAX;
    uint32_t minSwapPosA = 0;
    uint32_t minSwapPosB = 0;


    // Restart the search
    if(update) {
        taboListPos = 0;
        std::fill(taboListH, taboListH + TABO_LIST_LEN, UINT32_MAX);
        std::fill(taboListL, taboListL + TABO_LIST_LEN, UINT32_MAX);
        if(_bestGene != nullptr) {
            std::copy(_bestGene, _bestGene + _geneLength, _searchGene);
            _searchBestFitness = _bestFitness;
        } else {
            _searchBestFitness = DBL_MAX;
        }
        return false;
    }

    for(uint32_t sampIter = 0; sampIter < SAMPLE_POINTS; sampIter++) {
        uint32_t posA = _rng() % (_geneLength - 1);
        uint32_t posB = _rng() % (_geneLength - posA - 1) + posA + 1;
        std::swap(_searchGene[posA], _searchGene[posB]);
        fitness[sampIter] = _singleFitnessCal(_searchGene);
        std::swap(_searchGene[posA], _searchGene[posB]);
        uint32_t geneL = std::min(_searchGene[posA], _searchGene[posB]);
        uint32_t geneH = std::max(_searchGene[posA], _searchGene[posB]);
        uint32_t findPos = std::find(taboListL, taboListL + TABO_LIST_LEN, geneL) - taboListL;
        if(findPos != TABO_LIST_LEN && taboListH[findPos] == geneH) {
            if(fitness[sampIter] <= _bestFitness)
                fitness[sampIter] = DBL_MAX;
        } else {
            // Update the tabo list
            taboListH[taboListPos] = geneH;
            taboListL[taboListPos] = geneL;
            if(taboListPos < TABO_LIST_LEN - 1) {
                taboListPos++;
            } else {
                taboListPos = 0;
            }
        }

        if(fitness[sampIter] < minFit) {
            minSwapPosA = posA;
            minSwapPosB = posB;
            minFit = fitness[sampIter];
        }
    }

    std::swap(_searchGene[minSwapPosA], _searchGene[minSwapPosB]);

    if(minFit < _searchBestFitness) {
        _searchBestFitness = minFit;
    }

    return minFit < _bestFitness;
}

SchGA::~SchGA() {
    delete [] _fitness;
    delete [] _gene;
    delete [] _nextGene;
    delete [] _searchGene;
    delete [] _initialFlightState->flightState;
    delete _initialFlightState;
    delete [] _taskParamTable;
    delete [] _taskTable->taskQueue;
    delete _taskTable;
}