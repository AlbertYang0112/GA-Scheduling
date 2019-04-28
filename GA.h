#ifndef GAHW_GA_H
#define GAHW_GA_H

#include <cstdint>
#include <algorithm>

class GA {

private:
    virtual void _generateInitGene() = 0;
    virtual void _fitnessCal() = 0;
    virtual void _cross() = 0;
    virtual void _mutation() = 0;

public:
    /*
     * Execute the GA algorithm.
     * Return the best gene and the best fitness through the
     *  reference.
     */
    virtual void evaluate(uint32_t iterations,
            uint32_t &bestGene, uint32_t &bestFitness) = 0;
};


#endif //GAHW_GA_H
