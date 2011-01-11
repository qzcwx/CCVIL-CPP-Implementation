#include "F4.h"

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

F4::F4(RunParameter runParam):Benchmarks(runParam){
	cout<<"F4 Class initialization"<<endl;
	dimension = runParam.dimension;
	m_havenextGaussian=0;
	lookup  = NULL;
	Ovector = NULL;
}

F4::~F4(){
 	delete[] Ovector;
	cout<<"F4 Class destroyed"<<endl;
}

double F4::compute(double*x){
  int    i;
  double result = 0.0;

  if(Ovector == NULL) {
    Ovector   = createShiftVector(dimension,minX,maxX);
    Pvector   = createPermVector(dimension);
    RotMatrix = createRotMatrix1D(nonSeparableGroupSize);
    lookupprepare();
  }

  for(i = 0; i < dimension; i++) {
    anotherz[i] = x[i] - Ovector[i];
  }

  for(i = 0; i < nonSeparableGroupSize; i++) {
    anotherz1[i] = anotherz[Pvector[i]];
  }

  for(i = nonSeparableGroupSize; i < dimension; i++) {
    anotherz2[i - nonSeparableGroupSize] = anotherz[Pvector[i]];
  }

  result = rot_elliptic(anotherz1,nonSeparableGroupSize) * 1e6 + elliptic(
    anotherz2,dimension - nonSeparableGroupSize);
  return(result);
}
