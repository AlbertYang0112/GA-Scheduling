#include "SchGA.h"

extern "C" {
    #include "dubins.h"
}


void SchGA::_generateInitGene() {
    for(uint32_t i = 0; i < _taskTable->totalNum * _population; i++) {
        _gene[i] = static_cast<uint32_t>(random()) % _initialFlightState->num;
    }
}

void SchGA::_fitnessCal() {

}
