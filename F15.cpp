#include "F15.h"
#include <stdio.h>

/**
 * Single-group Shifted and m-rotated Elliptic Function
 *
 * as defined in "Benchmark Functions for the CEC'2010 Special Session
 * and Competition on Large-Scale Global Optimization" by Ke Tang,
 * Xiaodong Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise
 * published as technical report on January 8, 2010 at Nature Inspired
 * Computation and Applications Laboratory (NICAL), School of Computer
 * Science and Technology, University of Science and Technology of China,
 * Hefei, Anhui, China.
 */

F15::F15(RunParameter runParam):Benchmarks(runParam){
	cout<<"F15 Class initialization"<<endl;
	dimension = runParam.dimension;
	m_havenextGaussian=0;
	Ovector = NULL;
}

F15::~F15(){
	delete[] Ovector;
	delete[] Pvector;
	delete[] RotMatrix;
	cout<<"F15 Class destroyed"<<endl;
}

double F15::compute(double*x){
	int i,k;
	double result=0.0;

	if(Ovector==NULL)
	{
		Ovector=createShiftVector(dimension,minX,maxX);
		Pvector=createPermVector(dimension);
		RotMatrix=createRotMatrix1D(nonSeparableGroupSize);
	}
	for(i=0;i<dimension;i++)
	{
		anotherz[i]=x[i]-Ovector[i];
	}
	for(k=1;k<=dimension/(nonSeparableGroupSize);k++)
	{
		result+=rot_rastrigin(anotherz,nonSeparableGroupSize,k);

	}
	return(result);
}

