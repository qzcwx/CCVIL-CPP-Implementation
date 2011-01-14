#ifndef _F12_H
#define _F12_H

#include "Benchmarks.h"

class F12:public Benchmarks{
protected:
	static const int minX = -100;
	static const int maxX = 100;
public:
	F12(RunParameter runParam);
	double compute(double* x) ;
	~F12();
};

#endif
