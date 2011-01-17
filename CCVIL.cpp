#include "CCVIL.h"

CCVIL::CCVIL(RunParameter* runParam){
    parents = new PopulationT<double>(runParam->NP,  ChromosomeT< double >(runParam->dimension));
	offsprings = new PopulationT<double>(runParam->NP,  ChromosomeT< double >(runParam->dimension));

}


void CCVIL::run(){
	parents->setMinimize();
	offsprings->setMinimize();
}

