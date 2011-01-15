#ifndef _RUNPARAMETER_H
#define _RUNPARAMETER_H

#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;

class RunParameter{
public:
	// Dimension of problems
	int dimension;

	// the amount of independent run
	int numOfRun;

	// initial population size
	int NP;

	// initial Group Size
	int initialGroupSize;


	// Sampling interval for plotting the convergence curve
	int samplingInterval;

	// initialized random seed
	int initRandomSeed;

	// group size for non-separable part of function
	int nonSeparableGroupSize;
	// Fitness check point
	vector<int> fitnessCheckPoint;

	// the IDes of benchmark functions to be tested in the experiment
	vector<int> functionToRun;

	// default constructor
	RunParameter();

	// default destructor
	~RunParameter();

};
#endif
