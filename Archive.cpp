#include "Archive.h"

Archive::Archive(unsigned inNP, unsigned inD){
	MAX_NP = inNP; 
	dimension = inD;
	pop = new PopulationT<double>(0, ChromosomeT<double>(dimension));
}

Archive::~Archive(){
	delete pop;
}

unsigned Archive::getNP(){
	return (*pop).size();
}

PopulationT<double>* Archive::getPop(){
	return pop;
}
