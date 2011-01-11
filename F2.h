#ifndef _F2_H
#define _F2_H

#include "Benchmarks.h"


class F2:public Benchmarks{
protected:
	static const int minX = -5;
	static const int maxX = 5;
	double *Ovector;

public:
	F2(RunParameter runParam);
	double compute(double* x) ;
	~F2();
};
#endif
