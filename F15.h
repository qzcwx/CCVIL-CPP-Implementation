#ifndef _F15_H
#define _F15_H

#include "Benchmarks.h"

class F15:public Benchmarks{
protected:
	static const int minX = -5;
	static const int maxX = 5;
public:
	F15(RunParameter runParam);
	double compute(double* x) ;
	~F15();
};

#endif
