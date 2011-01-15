#ifndef _F11_H
#define _F11_H

#include "Benchmarks.h"

class F11:public Benchmarks{
protected:
	static const int minX = -32;
	static const int maxX = 32;
public:
	F11(RunParameter runParam);
	double compute(double* x) ;
	~F11();
};

#endif
