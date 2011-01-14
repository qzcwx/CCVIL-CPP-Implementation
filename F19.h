#ifndef _F19_H
#define _F19_H

#include "Benchmarks.h"

class F19:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F19(RunParameter runParam);
	double compute(double* x) ;
	~F19();
};

#endif


