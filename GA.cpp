#include "GA.h"

GA::GA(uint32_t num, uint32_t *initGene):_num(num) {
    _fitness = new uint32_t[num];
    if(initGene == nullptr) {
        _gene = new uint32_t[num];
        _generateInitGene();
    } else {
        _gene = initGene;
        std::copy(initGene, initGene + num, _gene);
    }
}

