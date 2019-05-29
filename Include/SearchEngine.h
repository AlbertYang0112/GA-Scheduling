#ifndef SEARCH_ENGINE_H
#define SEARCH_ENGINE_H

#include <cstdint>
#include <cfloat>
#include <cmath>
#include <random>
#include <algorithm>

class SchGA;

typedef double_t (SchGA::*FIT_CAL)(uint32_t*);

class SearchEngine {
public:
    SearchEngine(
        uint32_t geneLen, 
        uint32_t samplePoints, 
        uint32_t tabooListLen,
        SchGA* baseGA
    );
    ~SearchEngine();
    void update(uint32_t* pGene, double_t fitness=DBL_MAX);
    void search(uint32_t iterations);
    bool isBetter();
    double_t bestFitness();
    SearchEngine(const SearchEngine& se);
    void copyTheBestTo(uint32_t* dst);
private:
    uint32_t SAMPLE_POINTS;
    uint32_t TABOO_LIST_LEN;
    uint32_t GENE_LEN;
    uint32_t* _tabooListH;
    uint32_t* _tabooListL;
    uint32_t _tabooListPos;
    uint32_t* _gene;
    uint32_t* _bestGene;
    double_t _bestFitness;
    std::mt19937 _rng;
    bool _bestGeneUpdated;
    SchGA* _baseGA;
    void _appendTabooList(uint32_t elemA, uint32_t elemB);
    bool _searchInTabooList(uint32_t elemA, uint32_t elemB);
};

inline bool SearchEngine::isBetter() {
    if(_bestGeneUpdated) {
        _bestGeneUpdated = false;
        return true;
    }
    return false;
}

inline double_t SearchEngine::bestFitness() {
    return _bestFitness;
}

inline void SearchEngine::copyTheBestTo(uint32_t* dst) {
    std::copy(_bestGene, _bestGene + GENE_LEN, dst);
}

#endif