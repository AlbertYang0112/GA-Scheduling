#ifndef GAHW_GA_H
#define GAHW_GA_H

#include <cstdint>
#include <cmath>
#include <algorithm>

class GA {

private:
    /*
     * Generate the initial population.
     */
    virtual void _generateInitGene() = 0;
    /*
     * Update the fitness of the population.
     */
    virtual void _fitnessCal() = 0;
    /*
     * Implement the interbreed method here.
     * Note:        The gene transferred here must be interbred.
     * Parameter:   parentA, B <- gene to be interbred.
     * Return:      childA, childB <- offspring returned through reference.
     */
    virtual void _cross(uint32_t parentA, uint32_t parentB, uint32_t &childA, uint32_t &childB) = 0;
    /*
     * Implement the mutation method.
     * Note:        The gene transferred here must be mutated.
     * Parameter:   child <- gene to be mutated
     * Return:      Gene mutated.
     */
    virtual uint32_t _mutation(uint32_t child) = 0;

public:
    /*
     * Execute the GA algorithm.
     * Note:        This method should be the ONLY access to the GA algorithm.
     * Parameter:   iterations <- the generation to be evaluated;
     * Return:      bestGene <- the best gene in the last generation;
     *              bestFitness <- the fitness of the best gene
     */
    virtual void evaluate(uint32_t iterations,
            uint32_t &bestGene, double_t &bestFitness) = 0;
};


#endif //GAHW_GA_H
