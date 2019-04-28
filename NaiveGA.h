#ifndef GAHW_NAIVEGA_H
#define GAHW_NAIVEGA_H

#include "GA.h"

class NaiveGA: public GA {
    /*
     * A simple implementation of GA algorithm.
     * Finding the minimum value of the f(x) = x * x;
     */

private:
    uint32_t _num;      // The size of population
    uint32_t* _gene;    // Point to the memory storing all the gene
    double_t* _fitness;   // Point to the memory storing the corresponding fitness
    double_t* _cumulateFitness;   // Cumulation of the fitness
    double_t _maxFitness;
    uint32_t _bestGene;
    uint32_t _crossRate;        // Scaled cross rate in range [0, RAND_MAX]
    uint32_t _mutationRate;     // Scaled mutation rate in range [0, RAND_MAX]

    void _generateInitGene() override;
    void _fitnessCal() override;
    void _cross(uint32_t parentA, uint32_t parentB, uint32_t &childA, uint32_t &childB) override;
    uint32_t _mutation(uint32_t child) override;

public:
    void evaluate(uint32_t iterations, uint32_t &bestGene, double_t &bestFitness) override;
    NaiveGA(uint32_t num, double_t crossRate, double_t mutationRate) {
        _num = num;
        _gene = new uint32_t [num];
        _fitness = new double_t [num];
        _cumulateFitness = new double_t [num];
        _generateInitGene();
        _crossRate = static_cast<uint32_t>(RAND_MAX * crossRate);
        _mutationRate = static_cast<uint32_t >(RAND_MAX * mutationRate);
        _maxFitness = 0;
        _bestGene = 0;
    }
};


#endif //GAHW_NAIVEGA_H
