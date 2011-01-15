#ifndef _F3_H
#define _F3_H

#include "Benchmarks.h"


class F3:public Benchmarks{
protected:
	static const int minX = -32;
	static const int maxX = 32;

public:
	F3(RunParameter runParam);
	double compute(double* x) ;
	~F3();
};
#endif
