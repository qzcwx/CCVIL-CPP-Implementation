#ifndef _F16_H
#define _F16_H

#include "Benchmarks.h"

class F16:public Benchmarks{
protected:
	static const int minX = -32;
	static const int maxX = 32;
public:
	F16(RunParameter runParam);
	double compute(double* x) ;
	~F16();
};

#endif

