#ifndef _CCVIL_H
#define _CCVIL_H

#include <EALib/PopulationT.h>
#include <vector>
#include "Benchmarks.h"

using namespace std;

class CCVIL{
protected:

	PopulationT<double> *parents,*offsprings;
	Benchmarks* fp;
	vector< vector<int> > groupInfo;

	void learningStage();
	void optimizationStage();

public:
	CCVIL(RunParameter* runParam);
	~CCVIL();
	void run();
	void setObjFunc(Benchmarks* infp){fp = infp;};
};
#endif
