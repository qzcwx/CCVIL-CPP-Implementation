#ifndef _F14_H
#define _F14_H

#include "Benchmarks.h"

class F14:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F14(RunParameter runParam);
	double compute(double* x) ;
	~F14();
};

#endif

