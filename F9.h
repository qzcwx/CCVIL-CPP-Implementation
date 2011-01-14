#ifndef _F9_H
#define _F9_H

#include "Benchmarks.h"

class F9:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F9(RunParameter runParam);
	double compute(double* x) ;
	~F9();
};

#endif
