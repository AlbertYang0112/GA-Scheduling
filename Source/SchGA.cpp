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
        double_t mutationRate,
        const char* recorderName
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
    _bestFitness = DBL_MAX;
    
    // Calculate the length of the gene
    _geneLength = _taskTable->totalNum + _initialFlightState->num - 1;
    
    _fitness = new double_t[_population];

    // Initialize the random generator
    std::random_device rd;
    _rng.seed(rd());

#ifdef ADVANCED_FEATURE
    // Initialize the search engine array
    _searchEngines.insert(
        _searchEngines.begin(), 
        SEARCH_ENGINE_NUM, 
        SearchEngine(
            _geneLength, 
            static_cast<uint32_t>(0.5*_geneLength), 
            TABOO_LIST_LEN, 
            this)
        );
#endif

    // Scale the cross rate and the mutation rate
    assert(crossRate > 0 && crossRate < 1);
    assert(mutationRate > 0 && mutationRate < 1);
    _crossRate = crossRate * _rng.max();
    _mutationRate = mutationRate * _rng.max();

    // Open the recorder file
    _recorder.open(recorderName);
    _recorder << "Iter," << "FeasibleCnt," << "BestFitness," 
        << "FitnessAvg," << "FitnessVar," << "DistantAvg," 
       << "DistantVar," << "Saved" << std::endl;
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
    DEBUG_BRIEF("Gene Length: %d\n", _geneLength);
}

double_t SchGA::_singleFitnessCal(uint32_t* pGene) {
    TASK_PARAMETER *pTask;

    DubinsPath path;
    double startPoint[3];
    double endPoint[3];

    double_t taskTime[_taskTable->totalNum];
    uint32_t taskFlight[_taskTable->totalNum];
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
            taskFlight[*pTempGene] = flight;
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
#ifdef ADVANCED_FEATURE
                // Soft-margin
                tempMaxTime = _timeCompute(taskTime, taskFlight);
                isFeasibleSolution = tempMaxTime != DBL_MAX;
                if(isFeasibleSolution) _saveCnt++;
#else
                isFeasibleSolution = false;
#endif
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

    _saveCnt = 0;

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
#ifdef ADVANCED_FEATURE
        if(_rng() < _rng.max() * 0.85) {
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
#else
        if (_fitness[selA] >= _fitness[selB]) {
            parentsNo[i] = selB;
        } else {
            parentsNo[i] = selA;
        }
#endif
    }
}

uint32_t SchGA::_selectByDistant(uint32_t** geneList, uint32_t minDist, uint32_t num) {
    uint32_t selGene[num];
    std::fill(selGene, selGene + num, UINT32_MAX);

    double_t minFit;

    uint32_t* pGene;
    uint32_t genNum;
    for(genNum = 0; genNum < num; genNum++) {
        bool avail;
        minFit = DBL_MAX;
        pGene = _gene;
        for(uint32_t i = 0; i < _population; i++) {
            // Make sure the available gene is not visited.
            avail = std::find(selGene, selGene + num, i) == selGene + num;

            // Make sure the gene is outside of keepout region of previously selected gene.
            for(uint32_t j = 0; j < genNum && avail; j++) {
                avail = _distance(geneList[j], pGene) > minDist;
            }

            if(avail) {
                // Select the best available gene
                //DEBUG("GET ONE for %d result", genNum);
                if(minFit >= _fitness[i]) {
                    minFit = _fitness[i];
                    geneList[genNum] = pGene;
                    selGene[genNum] = i;
                }
            }
            pGene += _geneLength;
        }
        if(selGene[genNum] == DBL_MAX) {
            // No available gene in this iteration.
            break;
        }
    }
    return genNum;
}

void SchGA::evaluate(
        uint32_t iterations,
        uint32_t &bestGene,
        double_t &bestFitness) {
    // Generate the initial gene
    _generateInitGene();

    double_t searchBestFitness;
    bool searchUpdated;
    uint32_t parentsNo[2];
    for(uint32_t iter_cnt = 0; iter_cnt < iterations; iter_cnt++) {
        _fitnessCal();

#ifdef ADVANCED_FEATURE
        if(iter_cnt == 0) {
            _searchEngines[0].update(_bestGene, _bestFitness);
        }
#endif
        if(_feasibleGeneCnt != 0) {
            // If we found a feasible solution

#ifdef ADVANCED_FEATURE
            if(iter_cnt % 50 == 0) {
                DEBUG_BRIEF("Update the Boost Genes\n");
                // Update the preserved genes
                uint32_t* distGene[6];
                uint32_t distSelNum;
                distSelNum = _selectByDistant(distGene, 20, PRESERVED_SLOT - 2);
                uint32_t copyI;
                for(copyI = 0; copyI < distSelNum; copyI++) {
                    std::copy(distGene[copyI], distGene[copyI] + _geneLength, _gene + (2 + copyI) * _geneLength);
                }
                for(; copyI < PRESERVED_SLOT - 2; copyI++) {
                    std::copy(_bestGene, _bestGene + _geneLength, _gene + (2 + copyI) * _geneLength);
                }
                for(copyI = 1; copyI < SEARCH_ENGINE_NUM; copyI++) {
                    _searchEngines[copyI].update(_gene + (copyI + 2) * _geneLength);
                }
                DEBUG_BRIEF("Distant Sel Found %d Available Genes\n", distSelNum);
                for(uint32_t distA = 0; distA < distSelNum; distA++) {
                    DEBUG_BRIEF("From %d: ", distA);
                    for(uint32_t distB = 1; distB < distSelNum; distB++) {
                        DEBUG_BRIEF("%d ", _distance(distGene[distA], distGene[distB]));
                    }
                    DEBUG_BRIEF("\n");
                }
            }

            // Activate the Search Engine
            searchUpdated = false;
            uint32_t bestEngineNo = 0;
            searchBestFitness = DBL_MAX;

            for(uint32_t engineNo = 0; engineNo < SEARCH_ENGINE_NUM; engineNo++) {
                bool updated;
                uint32_t engineSearchIter = 0;
                do {
                    _searchEngines[engineNo].search(10);
                    engineSearchIter += 10;
                    updated = _searchEngines[engineNo].isBetter();
                    if(updated) searchUpdated = true;
                } while(updated && engineSearchIter < 100);
                if(_searchEngines[engineNo].bestFitness() < searchBestFitness) {
                    searchBestFitness = _searchEngines[engineNo].bestFitness();
                    bestEngineNo = engineNo;
                }
                DEBUG_BRIEF("%d ", engineSearchIter);
            }
            DEBUG_BRIEF("Iterations\n");
            for(uint32_t engineNo = 1; engineNo < SEARCH_ENGINE_NUM; engineNo++) {
                if(_searchEngines[engineNo].isBetter())
                    _searchEngines[engineNo].copyTheBestTo(_gene + (engineNo + 1) * _geneLength);
            }

            if(searchUpdated) {
                DEBUG_BRIEF("Engine %d Found a better gene: %f\n", bestEngineNo, (double)searchBestFitness);
            }

            if(_bestFitnessUpdated) {
                DEBUG_BRIEF("GA Found a better gene: %f\n", (double)_bestFitness);
            }


            // Gene preserved slot: | Best Gene | 2nd Best Gene | Distant Search Slot 1-6 |
            if(_bestFitnessUpdated && searchUpdated) {
                // Put the gene with better fitness to the first place
                // Put the other gene to the secondary place
                if(_bestFitness < searchBestFitness) {
                    std::copy(_bestGene, _bestGene + _geneLength, _nextGene);
                    _searchEngines[bestEngineNo].copyTheBestTo(_nextGene + _geneLength);
                    _fitness[0] = _bestFitness;
                    _fitness[1] = searchBestFitness;
                    // Update the search engine
                    _searchEngines[0].update(_bestGene, _bestFitness);
                } else {
                    std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
                    _searchEngines[bestEngineNo].copyTheBestTo(_nextGene);
                    _fitness[0] = searchBestFitness;
                    _fitness[1] = _bestFitness;
                    _bestFitness = searchBestFitness;
                    _bestGene = _nextGene;
                    if(bestEngineNo != 0) {
                        _searchEngines[0].update(_bestGene, _bestFitness);
                    }
                }
            } else if(searchUpdated) {
                std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
                _searchEngines[bestEngineNo].copyTheBestTo(_nextGene);
                _fitness[0] = searchBestFitness;
                _fitness[1] = _bestFitness;
                _bestFitness = searchBestFitness;
                _bestGene = _nextGene;
                if(bestEngineNo != 0) {
                    _searchEngines[0].update(_bestGene, _bestFitness);
                }
            } else {
                if(_bestFitnessUpdated) {
                    _searchEngines[0].update(_bestGene, _bestFitness);
                }
                // Simply Preserve the best gene
                std::copy(_bestGene, _bestGene + _geneLength, _nextGene);
                std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
                _fitness[0] = _bestFitness;
                _fitness[1] = _bestFitness;
            }
#else
            std::copy(_bestGene, _bestGene + _geneLength, _nextGene);
            std::copy(_bestGene, _bestGene + _geneLength, _nextGene + _geneLength);
            _fitness[0] = _bestFitness;
            _fitness[1] = _bestFitness;
            if (_rng() < _mutationRate) {
                _mutation(1);
            }
#endif
        }
        _bestFitnessUpdated = false;

        for(uint32_t cross_cnt = 2; cross_cnt < _population; cross_cnt += 2) {
            _selectParents(parentsNo, 2);
            uint32_t cross_tmp = cross_cnt + 1;
            if(_fitness[parentsNo[0]] == DBL_MAX || _fitness[parentsNo[1]] == DBL_MAX) {
                std::copy(_gene + parentsNo[0] * _geneLength,
                          _gene + parentsNo[0] * _geneLength + _geneLength,
                          _nextGene + cross_cnt * _geneLength);
                std::copy(_gene + parentsNo[1] * _geneLength,
                          _gene + parentsNo[1] * _geneLength + _geneLength,
                          _nextGene + cross_cnt * _geneLength + _geneLength);
#ifdef ADVANCED_FEATURE
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
#endif
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
            double_t distAvg = distantToBestAverage();
            double_t distVar = distantToBestVar();
            DEBUG_BRIEF("Iteration: %d\n", iter_cnt);
            DEBUG_BRIEF("Feasible Gene: %d\n", _feasibleGeneCnt);
            if(avg != DBL_MAX) {
                DEBUG_BRIEF("Fitness Average: %f\n", static_cast<double>(avg));
                DEBUG_BRIEF("Fitness Variance: %E\n", static_cast<double>(var));
            }
            DEBUG_BRIEF("Average Distance to the Best Gene: %f\n", static_cast<double>(distAvg)); 
            DEBUG_BRIEF("Distance Variance to the Best Gene: %f\n", static_cast<double>(distVar)); 
            DEBUG_BRIEF("Seq Saved: %d\n", _saveCnt);
            if(_feasibleGeneCnt != 0) {
                DEBUG_BRIEF("Best Fitness %f\n", static_cast<double>(_bestFitness));
            }
            _recorder << iter_cnt << ',' << _feasibleGeneCnt << ',' << _bestFitness
                << ',' << avg << ',' << var << ',' << distAvg << ',' << distVar << ',' << _saveCnt<< std::endl;
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

uint32_t SchGA::_distance(uint32_t* pGeneA, uint32_t* pGeneB) {
    uint32_t tempGene[_geneLength];
    std::copy(pGeneA, pGeneA + _geneLength, tempGene);

    uint32_t dist = 0;
    pGeneA = tempGene;
    uint32_t pairPos;
    for(uint32_t startPos = 0; startPos < _geneLength - 1; startPos++) {
        if(*pGeneA != *pGeneB) {
            pairPos = static_cast<uint32_t>(std::find(pGeneA, tempGene + _geneLength, *pGeneB) - tempGene);
            if(pairPos == _geneLength) {
                DEBUG("Distant calculation error\n");
                continue;
            }
            dist++;
            std::swap(*pGeneA, tempGene[pairPos]);
        }
        pGeneA++;
        pGeneB++;
    }
    return dist;
}

double_t SchGA::distantToBestAverage() {
    double_t distSum = 0;
    uint32_t cnt = 0;
    for(uint32_t i = 0; i < _population; i++) {
        if(_fitness[i] != DBL_MAX) {
            distSum += _distance(_bestGene, _gene + i * _geneLength);
            cnt++;
        }
    }

    distSum /= cnt;

    return distSum;
}

double_t SchGA::distantToBestVar() {
    uint32_t distAvg = distantToBestAverage();
    double_t var = 0;
    uint32_t cnt = 0;
    int32_t prod;
    for(uint32_t i = 0; i < _population; i++) {
        if(_fitness[i] != DBL_MAX) {
            prod = static_cast<int32_t>(
                _distance(_bestGene, _gene + i * _geneLength)
                ) - distAvg;
            var += prod * prod;
            cnt++;
        }
    }
    var = sqrt(var / cnt);
    return var;
}

bool check(uint32_t arg, double_t* stamp, uint32_t* division, uint32_t num){
    for(uint32_t i=0;i<num;i++){
        if(arg==division[i])
            return true;
    }
    bool t1;
    if(arg>0){
        t1=(stamp[arg-1]<=stamp[arg]);
    }
    else
        t1=true;
    return t1;
}

uint32_t* argsort(double_t* array, uint32_t size){
    uint32_t len=size;
    uint32_t* array_index=new uint32_t[len];
    for(uint32_t i=0;i<len;i++)
        array_index[i]=i;
    std::sort(array_index, array_index+len,
              [&array](uint32_t pos1, uint32_t pos2) {return (array[pos1] < array[pos2]);});

    return array_index;
}

double_t SchGA::_timeCompute(double_t* timestamp, uint32_t* flightNo, uint32_t max_step){
    uint32_t num=this->_taskTable->totalNum;
    uint32_t separator_num=this->_taskTable->queueNum;
    uint32_t* separator=new uint32_t[separator_num];
    uint32_t temp=0;
    bool* circle_detect=new bool[num];

    for(uint32_t i=0;i<num;i++)
        circle_detect[i]=false;
    for(uint32_t i=0;i<separator_num;i++){
        separator[i]=temp;
        temp+=this->_taskTable->taskQueue[i].num;
    }
    uint32_t* stamp_seq= argsort(timestamp, num);
    uint32_t step=0;
    for(uint32_t i=0;i<num;i++){
        if(step>max_step){
            //DEBUG_BRIEF("End by max_step\n");
            return DBL_MAX;
        }
        bool t=check(stamp_seq[i], timestamp, separator, separator_num);
        if(!t){
            if(circle_detect[stamp_seq[i]]==true){
                //DEBUG_BRIEF("End by circle\n");
                return DBL_MAX;
            }
            else
                circle_detect[stamp_seq[i]]=true;
            step++;
            double_t diff=timestamp[stamp_seq[i]-1]-timestamp[stamp_seq[i]];
            double_t current=timestamp[stamp_seq[i]];
            for(uint32_t j=i;j<num;j++){
                if((flightNo[stamp_seq[j]]==flightNo[stamp_seq[i]])&&
                   (timestamp[stamp_seq[j]]>=current))
                    timestamp[stamp_seq[j]]+=diff;
            }
            //print log
            //DEBUG_BRIEF("Change:\n");
            //for(uint32_t j=0;j<num;j++)
            //    DEBUG_BRIEF("%f ", double(timestamp[j]));
            //DEBUG_BRIEF("\n");
            delete[] stamp_seq;
            stamp_seq= argsort(timestamp, num);
            i-=1;
        }
    }
    double_t time=timestamp[stamp_seq[num-1]];
    delete[] stamp_seq;
    //DEBUG_BRIEF("End by default\n");
    return time;
}

SchGA::~SchGA() {
    _recorder.close();
    delete [] _fitness;
    delete [] _gene;
    delete [] _nextGene;
    delete [] _initialFlightState->flightState;
    delete _initialFlightState;
    delete [] _taskParamTable;
    delete [] _taskTable->taskQueue;
    delete _taskTable;
}