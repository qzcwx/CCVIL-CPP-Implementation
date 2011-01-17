#ifndef _CCVIL_H
#define _CCVIL_H

#include "Benchmarks.h"
#include <EALib/PopulationT.h>

class CCVIL{
protect:
	Benchmarks* fp;
	vector<vector<unsigned>>  groupInfo;
	learningStage();
	optimizationStage();
	PopulationT<double>* parents, offsprings;

public:
	void setObjFunc(Benchmarks* fp){this.fp = fp};
	CCVIL(RunParameter* runParam);
	~CCVIL();
	run();
}
#endif
