#ifndef _F20_H
#define _F20_H

#include "Benchmarks.h"

class F20:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F20(RunParameter runParam);
	double compute(double* x) ;
	~F20();
};

#endif

