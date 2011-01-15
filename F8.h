#ifndef _F8_H
#define _F8_H

#include "Benchmarks.h"

class F8:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F8(RunParameter runParam);
	double compute(double* x) ;
	~F8();
};

#endif
