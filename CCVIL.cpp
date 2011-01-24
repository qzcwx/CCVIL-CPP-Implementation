#include "CCVIL.h"

CCVIL::CCVIL(RunParameter* runParam){
	param = runParam;
	p = randPerm(param->dimension);

	Rng::seed(param->initRandomSeed);

	lookUpGroup = new unsigned[param->dimension];

	MaxFitEval = runParam->fitnessCheckPoint[(*runParam).fitnessCheckPoint.size()-1];
	cout<<"Max Fitness Evaluation = "<<MaxFitEval<<endl;

	lowerThreshold = runParam->lowerThreshold;
	upperThreshold = min(round(MaxFitEval*0.6/(runParam->dimension*((1+1)*(3)+1))), 800.0);
	cout<<"Lower threshold = "<<lowerThreshold<<", Upper threshold = "<<upperThreshold<<endl;

	//	cout<<"Init Grouping"<<endl;
	// initialize the groupInfo
	for (unsigned i = 0; i<runParam->dimension; i++){
		vector<unsigned> tempVec;
		tempVec.push_back(i);
		lookUpGroup[i] = i;
		//		printf("%d\n",p[i]);
		groupInfo.push_back(tempVec);
	}
	printf("Finish Initialization on Group");
}

CCVIL::~CCVIL(){
	delete bestCand;
	delete[] p;
	delete[] lookUpGroup;
}

void CCVIL::run(){
	learningStage();
	// optimizationStage();
}

/* 
 * procedure of learning stage, update the "groupInfo"
 */
void CCVIL::learningStage(){
	cout<<"Learning Stage ... "<<endl;
	unsigned cycle = 1, lastCycleIndex = 0;
	bool learnStageFlag, needCapture, isSameGroup, separableFunc = true; // assume every benchmark function is separable at the first beginning

	bestCand = new PopulationT<double>(1, ChromosomeT<double>(param->dimension));
	bestCand->setMinimize();
	(*bestCand)[0][0].initialize(fp->getMinX(), fp->getMaxX());

	popGenerate();

	while ( (learnStageFlag = !(cycle > upperThreshold 
					|| (cycle > lowerThreshold && groupInfo.size() == param->dimension) 
					|| groupInfo.size() == 1))==true){
		// start a new cycle
		printf("Cycle = %d\n", cycle);
		for (unsigned i=0; i<param->dimension; i++) {
			if (i == 0){
				// start a new phase for each cycle, each phase checking one dimension, i.e., p[i]
				printf("Population Re-initialization & Issue Random Permutation\n");

				popInit();
				p = randPerm(param->dimension);
				printf("Random Permutation p:\n");
				printArray(p, param->dimension);
				lastCycleIndex = 0;
			}

			needCapture = groupInfo.size()!=1 && ((cycle<= lowerThreshold) ||(separableFunc == false && cycle <= upperThreshold)) && lastCycleIndex!=0;
			//	printf("Need Capture ? %d\n", needCapture);
			
			if (lastCycleIndex!=0){
				//	to decide whether current dimesion are in the same group with last dimension
				isSameGroup = sameGroup(p[i], lastCycleIndex);
			}

			if (lastCycleIndex == 0 || isSameGroup == false){
				printf("i = %d, p[i] = %d\n", i, p[i]);
				// if current dimension and last optimized index are in the different group, then optimize on current dimension
				rJADECC(p[i], true /*learnStage = true*/);
			} else{
				needCapture = false;
			}

			// begin interaction capture process
			// TODO: Need to deal with possible inaccurate fitness value influencing the capturing result
			if (needCapture == true){
				captureInter(p[i],lastCycleIndex);
				separableFunc = false;
			}

			// record last optimized dimension
			if ( lastCycleIndex==0 || isSameGroup == false ){
				lastCycleIndex = p[i];
			}
		}
		cycle++;
	}

	/*
	// initialize population
	for (unsigned i = 0; i < parents->size(); ++i) {
		(*parents)[ i ][ 0 ].initialize(fp->getMinX(), fp->getMaxX());
	}

	for (unsigned i = 0; i < parents->size(); ++i){
		printf("pop %d = %1.20E\n", i, fp->compute((*parents)[ i ][ 0 ]));
		(*parents)[ i ].setFitness(fp->compute((*parents)[ i ][ 0 ]));
	}
	*/
}

/* 
 * rJADECC: the internal optimizer of CCVIL
 * Optimize on one specific dimension for each run
 *
 * index: the ID of group
 */
void CCVIL::rJADECC(unsigned index, bool learnStageFlag){
	/*
	   printf("Optimizing on Dimension ");
	   for (unsigned i=0; i<index.size(); i++){
	   if (i!=index.size()-1)
	   printf("%d,\t",index[i]);
	   else
	   printf("%d\n",index[i]);
	   }
	   */

	/*
	 * LB: 	Lower Bound
	 * UB: 	Upper Bound
	 * D: 	Dimension of current group
	 * NP:	Number of Population
	 * G:	Generation Limit for optimizing current group
	 */
	unsigned LB = fp->getMinX(), UB = fp->getMaxX(), D = pop[index][0].size(), NP= pop[index].size(), G, g;
	double Fm, CRm, c = param->c, p = param->p;
	double *F, *CR;

	if (learnStageFlag == true){
		G = 1;
	} else if (pop.size() == 1){
		G = INT_MAX;
	} else{
		G = min(D+5, (unsigned)500);
	}

	// printf("LB = %d, UB = %d, D = %d, NP = %d, G = %d\n", LB, UB, D, NP, G);

	if ( learnStageFlag == true ){
		Fm = 0.5;
		CRm = 0.95;
	} else {
		CRm = groupCR[index];
		Fm = groupF[index];
	}

	F = new double[NP];
	CR = new double[NP];

	randFCR(NP, CRm, 0.1, Fm, 0.1, F, CR);

	if (param->Afactor > 0){
		//TODO: Archive initialization
	}

	/***************************** Population Initialization and evaluation **************************/
	g = 1;
	PopulationT<double> parents(NP, ChromosomeT<double>(param->dimension));
	parents.setMinimize();

	printf("The whole population\n");
	printWholePop();

	printf("Best Candidate at the beginning of each phase\n");
	printPopulation(*bestCand);

	for(unsigned i=0; i < parents.size(); i++){
		// run on individuals
		if (learnStageFlag == true)  {
			// For learning stage, groupInfo may merge internally, while pop always maintain one-dimensional group
			// jth dimension's group belongs to the same group as current optimized group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// run on dimensions
				printf("i = %d, j = %d, index = %d\n", i, j, index);
				if (j == index)  {
					parents[i][0][j] = (pop[index])[i][0][0];
				}else if (learnStageFlag == false && lookUpGroup[j] == index){
					// For optimization stage,

				} else {
					parents[i][0][j] = (*bestCand)[0][0][j];
				}
			}
		} else{
			// For optimization stage, use goupInfo to enumerate all elements in the same group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// not belongs to the current optimized group
				if (lookUpGroup[j] != index){
					parents[i][0][j] = (*bestCand)[0][0][j];
				}

				// deal with the current optimized dimensions
				for (unsigned j=0; j< groupInfo[index].size(); j++){
					unsigned I = groupInfo[index][j];
					parents[i][0][I] = (pop[index])[i][0][j];
				}
			}
		}
	}

	printf("Parents Population\n");
	printPopulation(parents);
	
	for (unsigned i=0; i<parents.size(); i++){
		parents[i].setFitness( fp->compute(parents[i][0]) );
		printf("indiv %d, fit = %f\n", i, parents[i].fitnessValue()) ;
	}

	unsigned bestIndex = parents.bestIndex();
	printf(" Best Index = %d, Value =  %f\n",parents.bestIndex(), parents.best().fitnessValue());

	// after iterations, update the optimized optimized dimensions to bestCand and pop simultaneously
	
	if (learnStageFlag == true) {
		for (unsigned j=0; j<parents[0][0].size(); j++){
			// update bestCand
			if (j == index) {
				(*bestCand)[0][0][j] = parents[bestIndex][0][j];
			}
			// update pop, the entire external population
			for(unsigned i=0; i<NP; i++){
				(pop[index])[i][0][0]=parents[i][0][index];
			}
		}
	} else {
		for (unsigned j=0; j< groupInfo[index].size(); j++){
			unsigned I = groupInfo[index][j];
			(*bestCand)[0][0][I] = parents[bestIndex][0][I];
			for (unsigned i=0; i<NP; i++){
				(pop[index])[i][0][j] = parents[i][0][I];
			}
		}
	}

	printf("Best Candidate after update\n");
	printPopulation(*bestCand);

	delete[] F;
	delete[] CR;
}

void CCVIL::printArray(double* a, unsigned D){
	for(unsigned i = 0; i<D; i++){
		printf("%f\n",a[i]);
	}
}

void CCVIL::printArray(unsigned* a, unsigned D){
	for(unsigned i = 0; i<D; i++){
		printf("%d\n",a[i]);
	}
}

void CCVIL::printPopulation(PopulationT<double> &printPop){
	printf("Dimension of printed Population = %d\n", (printPop)[0][0].size());
	for (unsigned j=0; j<printPop.size(); j++){
		for (unsigned k=0; k<printPop[j][0].size(); k++){
			printf("%f\t",(printPop)[j][0][k]);
		}
		printf("\n");
	}
	cout<<endl;
}

void CCVIL::printWholePop(){
	printf("Whole Population\n");
	for (unsigned i=0; i<pop.size(); i++){
		printf("********** Subpopulation %d **********\n", i);
		for (unsigned j=0; j<(pop[i]).size(); j++){
			for (unsigned k=0; k<(pop[i])[j][0].size(); k++){
				printf("%f\t",(pop[i])[j][0][k]);
			}
			printf("\n");
		}
	}
}

void CCVIL::randFCR(unsigned NP, double CRm, double CRsigma, double Fm, double Fsigma, double* &F, double* &CR){
	//	norRnd->mean(NP);
	//	norRnd->variance(1);
	Normal norRnd;
	for (unsigned i=0; i<NP; i++){
		CR[i] = CRm + CRsigma * (norRnd)(); 
		CR[i] = min( 1.0 , max( 0.0, CR[i] ) ); // truncated to [0 1]
	}

	for (unsigned i=0; i<NP; i++){
		F[i] = cauchyRnd(Fm, Fsigma);
		F[i] = min(1.0, F[i]); // truncated to [-INF, 1]
	}

	// instead of truncation for dealing with the lower bound -1, we regenerate invalid random numbers
	vector<unsigned> pos1, pos2;
	pos1.clear();
	pos2.clear();
	for (unsigned i=0; i<NP; i++){
		if (F[i]<=0){
			pos1.push_back(i);
		}
		else if (F[i]>=1){
			pos2.push_back(i);
		}
	}

	while ( !pos1.empty() || !pos2.empty() ) {
		for (unsigned i=0; i<pos1.size(); i++){
			F[pos1[i]] = cauchyRnd(Fm, Fsigma);
			F[pos1[i]] = min(1.0, F[pos1[i]]); // truncated to [-INF, 1]
		}

		for (unsigned i=0; i<pos2.size(); i++){
			F[pos2[i]] = cauchyRnd(Fm, Fsigma);
			F[pos2[i]] = min(1.0, F[pos2[i]]); // truncated to [-INF, 1]
		}

		pos1.clear();
		pos2.clear();
		for (unsigned i=0; i<NP; i++){
			if (F[i]<=0){
				pos1.push_back(i);
			}
			else if (F[i]>=1){
				pos2.push_back(i);
			}
		}
	}

	for (unsigned i=0; i<NP; i++){
		if (F[i]>=1){
			cerr<<"ERROR: F = 1"<<endl;
			exit(-1);
		}
	}
}

double CCVIL::cauchyRnd(double mu, double delta){
	Uniform uniRnd;
	return ( mu + delta * tan(PI * uniRnd()) );
}

void CCVIL::popGenerate(){
	for (unsigned i=0; i<param->dimension; i++){
		PopulationT<double> tempPop(param->NP, ChromosomeT<double>(param->initialGroupSize));
		tempPop.setMinimize();
		pop.push_back(tempPop);
	}
}

void CCVIL::popInit(){
	for (unsigned i=0; i<pop.size(); i++){
		for(unsigned j=0; j<pop[i].size(); j++){
			(pop[i])[j][0].initialize(fp->getMinX(), fp->getMaxX());
		}
	}
}

/*
void CCVIL::popInit(PopulationT<double> *inPop){
	// initialize population
	for (unsigned i = 0; i < inPop->size(); ++i) {
		(*inPop)[ i ][ 0 ].initialize(fp->getMinX(), fp->getMaxX());
	}
}
*/



void CCVIL::captureInter(unsigned curDim, unsigned lastDim){
	printf("Capture Interaction Between %d & %d\n", curDim, lastDim);
}

/*
 * check whether the two variable belong to the same group or not
 */
bool CCVIL::sameGroup(unsigned v1, unsigned v2){
	if (lookUpGroup[v1] == lookUpGroup[v2]){
		return true;
	}else{
		return false;
	}
}

/*
 * generate an random permutation with the length N
 */
unsigned* CCVIL::randPerm(unsigned N)
{
	unsigned* p = new unsigned[N];
	for (unsigned i = 0; i < N; ++i) {
		int j = rand() % (i + 1);
		p[i] = p[j];
		p[j] = i;
	}
	return p;
}
