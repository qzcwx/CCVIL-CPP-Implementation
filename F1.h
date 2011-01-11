#ifndef _F1_H
#define _F1_H

#include "Benchmarks.h"


class F1:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
	double *Ovector;

public:
	F1(RunParameter runParam);
	double compute(double* x) ;
	~F1();
};
#endif
