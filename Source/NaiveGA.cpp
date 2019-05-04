#include "../Include/NaiveGA.h"

void NaiveGA::_generateInitGene() {
    uint32_t* pGene = _gene;
    for(uint32_t i = 0; i < _num; i++) {
        *pGene++ = static_cast<uint32_t >(random()) % 100;
    }
}

void NaiveGA::_fitnessCal() {
    uint32_t* pGene = _gene;
    uint32_t* pFitness = _fitness;
    uint32_t* pCFitness = _cumulateFitness;
    for(uint32_t i = 0; i < _num; i++) {
        *pFitness = 1000 - (*pGene);
        if(*pFitness > 1000) {
            *pFitness = 0;
        }

        if(*pFitness > _maxFitness) {
            _maxFitness = *pFitness;
            _bestGene = *pGene;
        }

        if(i == 0) {
            *pCFitness = *pFitness;
        } else {
            *pCFitness = *(pCFitness - 1) + *pFitness;
        }
        pFitness++;
        pCFitness++;
        pGene++;
    }
}

void NaiveGA::_cross(uint32_t parentA, uint32_t parentB, uint32_t &childA, uint32_t &childB) {
    uint32_t crossPos = static_cast<uint32_t >(random()) % 31 + 1;  // cross in range [1, 31]
    uint32_t crossMask = (1U << crossPos) - 1U;
    childA = (parentA & crossMask) | (parentB & (~crossMask));
    childB = (parentB & crossMask) | (parentA & (~crossMask));
}

void NaiveGA::evaluate(uint32_t iterations, uint32_t &bestGene, double_t &bestFitness) {
    uint32_t children[_num];
    uint32_t childA, childB;
    for(uint32_t iter = 0; iter < iterations; iter++) {
        _fitnessCal();
        bestGene = _bestGene;
        bestFitness = _maxFitness;
#pragma omp parallel for num_threads(4)
        for (uint32_t child = 0; child < _num; child += 2) {
            uint32_t rand1 = random() % _cumulateFitness[_num - 1];
            uint32_t rand2 = random() % _cumulateFitness[_num - 1];
            uint32_t parent1 = _gene[
                    std::lower_bound(_cumulateFitness, _cumulateFitness + _num, rand1) - _cumulateFitness
            ];
            uint32_t parent2 = _gene[
                    std::lower_bound(_cumulateFitness, _cumulateFitness + _num, rand2) - _cumulateFitness
            ];
            if (random() < _crossRate) {
                _cross(parent1, parent2, childA, childB);
            } else {
                childA = parent1;
                childB = parent2;
            }
            if (random() < _mutationRate) {
                childA = _mutation(childA);
            }
            if (random() < _mutationRate) {
                childB = _mutation(childB);
            }
            children[child] = childA;
            children[child + 1] = childB;
        }
        std::copy(children, children + _num, _gene);
    }
}

uint32_t NaiveGA::_mutation(uint32_t child) {
    uint32_t bitSel = static_cast<uint32_t >(random()) % 32;
    return child ^ (1U << bitSel);
}
