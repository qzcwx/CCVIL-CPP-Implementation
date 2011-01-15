#ifndef _F4_H
#define _F4_H

#include "Benchmarks.h"

class F4:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F4(RunParameter runParam);
	double compute(double* x) ;
	~F4();
};

#endif
