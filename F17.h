#ifndef _F17_H
#define _F17_H

#include "Benchmarks.h"

class F17:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F17(RunParameter runParam);
	double compute(double* x) ;
	~F17();
};

#endif

