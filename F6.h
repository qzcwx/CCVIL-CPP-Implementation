#ifndef _F6_H
#define _F6_H

#include "Benchmarks.h"

class F6:public Benchmarks{
protected:
	static const int minX = -32;
	static const int maxX = 32;
public:
	F6(RunParameter runParam);
	double compute(double* x) ;
	~F6();
};

#endif
