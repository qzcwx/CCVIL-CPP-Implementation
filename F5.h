#ifndef _F5_H
#define _F5_H

#include "Benchmarks.h"

class F5:public Benchmarks{
protected:
	static const int minX = -5;
	static const int maxX = 5;
public:
	F5(RunParameter runParam);
	double compute(double* x) ;
	~F5();
};

#endif
