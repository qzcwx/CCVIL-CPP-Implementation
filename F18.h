#ifndef _F18_H
#define _F18_H

#include "Benchmarks.h"

class F18:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F18(RunParameter runParam);
	double compute(double* x) ;
	~F18();
};

#endif

