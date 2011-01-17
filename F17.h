#ifndef _F17_H
#define _F17_H

#include "Benchmarks.h"

class F17:public Benchmarks{
protected:
public:
	F17(RunParameter runParam);
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F17();
};

#endif

