#ifndef _CCVIL_H
#define _CCVIL_H

#include <EALib/PopulationT.h>
#include <EALib/IndividualT.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <climits>
#include <Rng/Normal.h>
#include <Rng/Uniform.h>
//#include <GlobalRng.h>

#include "Benchmarks.h"
#include "Archive.h"

using namespace std;

class CCVIL{
protected:
	vector< PopulationT<double> > pop;
	vector< vector<unsigned> > groupInfo;

	PopulationT<double>* bestCand;
	Benchmarks* fp;
	unsigned* p; // random permutation
	unsigned* lookUpGroup;
	unsigned MaxFitEval; // total fitness evaluation limitation for optimization
	unsigned curFitEval; // current fitness evaluation limitation for optimization
	unsigned* randPerm(unsigned N);
	unsigned lowerThreshold;
	unsigned upperThreshold;
	RunParameter* param;
	double* groupCR;
	double* groupF;
	long fes;

	void learningStage();
	void optimizationStage();
	bool sameGroup(unsigned v1, unsigned v2);
	void rJADECC(unsigned index, bool learnStageFlag);
	void captureInter(unsigned curDim, unsigned lastDim);
	void popInit(PopulationT<double> *inPop);
	void popInit();
	void popGenerate();
	void randFCR(unsigned NP, double CRm, double CRsigma, double Fm, double Fsigma, double *&F, double *&CR);
	double cauchyRnd(double mu, double delta);
	void printArray(double* a, unsigned D);
	void printArray(unsigned* a, unsigned D);
	void printPopulation(PopulationT<double> &printPop);
	void printWholePop();
	double sum(vector<double> vec);
	vector<double> dotMultiply(vector<double> v1, vector<double> v2);
	void gnR1R2(unsigned NP1, unsigned NP2, unsigned *r1, unsigned *r2);
	PopulationT<double> combinePopulation(PopulationT<double> p1, PopulationT<double> p2);
	void findPbestIndex(PopulationT<double> inPop, unsigned pNP, unsigned* indBest);

public:
	CCVIL(RunParameter* runParam);
	~CCVIL();
	void run();
	void setObjFunc(Benchmarks* infp){fp = infp;};
};
#endif
