#ifndef _F14_H
#define _F14_H

#include "Benchmarks.h"

class F14:public Benchmarks{
protected:
public:
	F14(RunParameter runParam);
	double compute(double* x) ;
	double compute(vector<double> x) ;
	~F14();
};

#endif

