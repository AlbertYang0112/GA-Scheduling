#include "SearchEngine.h"
#include "SchGA.h"
#include <random>
#include <cassert>
#include <algorithm>

SearchEngine::SearchEngine(
    uint32_t geneLen, 
    uint32_t samplePoints, 
    uint32_t tabooListLen,
    SchGA* baseGA):
    GENE_LEN(geneLen),
    SAMPLE_POINTS(samplePoints),
    TABOO_LIST_LEN(tabooListLen),
    _baseGA(baseGA){
    assert(baseGA != nullptr);
    _tabooListH = new uint32_t[TABOO_LIST_LEN];
    _tabooListL = new uint32_t[TABOO_LIST_LEN];
    _gene = new uint32_t[GENE_LEN];
    _bestGene = new uint32_t[GENE_LEN];
    _bestGeneUpdated = false;
    _bestFitness = DBL_MAX;
    _tabooListPos = 0;
    std::random_device rd;
    _rng.seed(rd());
}

SearchEngine::SearchEngine(const SearchEngine& se):
    GENE_LEN(se.GENE_LEN),
    SAMPLE_POINTS(se.SAMPLE_POINTS),
    TABOO_LIST_LEN(se.TABOO_LIST_LEN),
    _baseGA(se._baseGA) {
    _tabooListH = new uint32_t[TABOO_LIST_LEN];
    _tabooListL = new uint32_t[TABOO_LIST_LEN];
    _gene = new uint32_t[GENE_LEN];
    _bestGene = new uint32_t[GENE_LEN];
    _bestGeneUpdated = false;
    _bestFitness = DBL_MAX;
    _tabooListPos = 0;
    std::random_device rd;
    _rng.seed(rd());
}

void SearchEngine::update(uint32_t* pGene, double_t fitness) {
    std::fill(_tabooListL, _tabooListL + TABOO_LIST_LEN, UINT32_MAX);
    std::fill(_tabooListH, _tabooListH + TABOO_LIST_LEN, UINT32_MAX);
    _tabooListPos = 0;
    _bestFitness = DBL_MAX;
    std::copy(pGene, pGene + GENE_LEN, _gene);
    std::copy(_gene, _gene + GENE_LEN, _bestGene);
    if(fitness == DBL_MAX) {
        _bestFitness = _baseGA->_singleFitnessCal(_bestGene);
    }
    _bestGeneUpdated = false;
}

void SearchEngine::search(uint32_t iterations) {
    double_t fitness[SAMPLE_POINTS];
    double_t minFit;
    uint32_t minSwapPosA, minSwapPosB;
    for(uint32_t iter = 0; iter < iterations; iter++) {
        minSwapPosA = 0;
        minSwapPosB = 0;
        minFit = DBL_MAX;
        for(uint32_t sampIter = 0; sampIter < SAMPLE_POINTS; sampIter++) {
            uint32_t posA = _rng() % (GENE_LEN - 1);
            uint32_t posB = _rng() % (GENE_LEN - posA - 1) + posA + 1;

            std::swap(_gene[posA], _gene[posB]);
            fitness[sampIter] = _baseGA->_singleFitnessCal(_gene);
            std::swap(_gene[posA], _gene[posB]);

            if(_searchInTabooList(_gene[posA], _gene[posB]) && fitness[sampIter] >= _bestFitness) {
                fitness[sampIter] = DBL_MAX;
            } else {
                _appendTabooList(_gene[posA], _gene[posB]);
            }

            if(fitness[sampIter] < minFit) {
                minSwapPosA = posA;
                minSwapPosB = posB;
                minFit = fitness[sampIter];
            }
        }
        std::swap(_gene[minSwapPosA], _gene[minSwapPosB]);
        if(minFit < _bestFitness) {
            _bestFitness = minFit;
            _bestGeneUpdated = true;
            std::copy(_gene, _gene + GENE_LEN, _bestGene);
        }
    }
}

void SearchEngine::_appendTabooList(uint32_t elemA, uint32_t elemB) {
    if(elemA < elemB) {
        std::swap(elemA, elemB);
    }
    _tabooListH[_tabooListPos] = elemA;
    _tabooListL[_tabooListPos] = elemB;
    if(_tabooListPos < TABOO_LIST_LEN - 1) {
        _tabooListPos++;
    } else {
        _tabooListPos = 0;
    }
}

bool SearchEngine::_searchInTabooList(uint32_t elemA, uint32_t elemB) {
    if(elemA < elemB) {
        std::swap(elemA, elemB);
    }
    uint32_t findPos = static_cast<uint32_t>(std::find(
        _tabooListH, _tabooListH + TABOO_LIST_LEN, elemA) - _tabooListH);
    return findPos != TABOO_LIST_LEN && _tabooListL[findPos] == elemB;
}

SearchEngine::~SearchEngine() {
    delete [] _tabooListH;
    delete [] _tabooListL;
    delete [] _gene;
    delete [] _bestGene;
}
    