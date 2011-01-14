
#ifndef _F7_H
#define _F7_H

#include "Benchmarks.h"

class F7:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F7(RunParameter runParam);
	double compute(double* x) ;
	~F7();
};

#endif
