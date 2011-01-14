#ifndef _F13_H
#define _F13_H

#include "Benchmarks.h"

class F13:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F13(RunParameter runParam);
	double compute(double* x) ;
	~F13();
};

#endif

