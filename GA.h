#ifndef GAHW_GA_H
#define GAHW_GA_H

#include <cstdint>
#include <algorithm>

class GA {

private:
    uint32_t* _gene;
    uint32_t* _fitness;
    uint32_t _num;

    virtual void _generateInitGene() = 0;
    virtual void _fitnessCal() = 0;
    virtual void _cross() = 0;
    virtual void _mutation() = 0;

public:
    explicit GA(uint32_t num, uint32_t *initGene = nullptr);
    ~GA() {
        delete [] _gene;
        delete [] _fitness;
    }

    /*
     * Execute the GA algorithm.
     * Return the best gene and the best fitness through the
     *  reference.
     */
    virtual void evaluate(uint32_t iterations,
            uint32_t &bestGene, uint32_t &bestFitness) = 0;
};


#endif //GAHW_GA_H
