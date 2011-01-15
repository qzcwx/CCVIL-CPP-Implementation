#ifndef _F10_H
#define _F10_H

#include "Benchmarks.h"

class F10:public Benchmarks{
protected:
	static const int minX = -5;
	static const int maxX = 5;
public:
	F10(RunParameter runParam);
	double compute(double* x) ;
	~F10();
};

#endif
