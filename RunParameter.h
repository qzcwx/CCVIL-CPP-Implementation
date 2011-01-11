#ifndef _RUNPARAMETER_H
#define _RUNPARAMETER_H

#include "Header.h"
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>

class RunParameter{
public:
	// Dimension of problems
	int dimension;

	// the IDes of benchmark functions to be tested in the experiment
	vector<int> functionToRun;

	// the amount of independent run
	int numOfRun;

	// initial population size
	int NP;

	// initial Group Size
	int initialGroupSize;

	// Fitness check point
	vector<int> fitnessCheckPoint;

	// Sampling interval for plotting the convergence curve
	int samplingInterval;

	// initialized random seed
	int initRandomSeed;

	// group size for non-separable part of function
	int nonSeparableGroupSize;

	// default constructor
	RunParameter();

	// default destructor
	~RunParameter();
};
#endif
