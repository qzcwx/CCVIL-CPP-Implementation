#include "F16.h"
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

F16::F16(RunParameter runParam):Benchmarks(runParam){
	cout<<"F16 Class initialization"<<endl;
	dimension = runParam.dimension;
	m_havenextGaussian=0;
	Ovector = NULL;
}

F16::~F16(){
	delete[] Ovector;
	delete[] Pvector;
	// delete 2D array
	int i;
	for(i=0;i<dimension/(nonSeparableGroupSize);i++){
		delete[] MultiRotMatrix1D[i];
	}
	delete[] MultiRotMatrix1D;

	cout<<"F16 Class destroyed"<<endl;
}

double F16::compute(double*x){
	int i,k;
	double result=0.0;

	if(Ovector==NULL)
	{
		Ovector=createShiftVector(dimension,minX,maxX);
		Pvector=createPermVector(dimension);
		MultiRotMatrix1D=createMultiRotateMatrix1D(nonSeparableGroupSize,dimension/(nonSeparableGroupSize));
	}
	for(i=0;i<dimension;i++)
	{
		anotherz[i]=x[i]-Ovector[i];
	}
	for(k=1;k<=dimension/(nonSeparableGroupSize);k++)
	{
		result+=rot_ackley(anotherz,nonSeparableGroupSize,k);
	}
	return(result);
}

