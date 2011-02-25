/*
 * =====================================================================================
 *
 *       Filename:  CCVIL.cpp
 *
 *    Description:  The main source file of CCVIL algorithm's implementation, it includes
 *    				two stages:
 *    					1) Learning Stage
 *    					2) Optimization Stage
 *
 *        Version:  1.0
 *        Created:  02/24/2011 07:56:20 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wenxiang Chen (http://cs-chen.net), chenwx.ustc@gmail.com
 *        Company:  Nature Inspired Computation and Application Laboratory (NICAL), USTC
 *
 * =====================================================================================
 */
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

	// initialize the groupInfo
	for (unsigned i = 0; i<runParam->dimension; i++){
		vector<unsigned> tempVec;
		tempVec.push_back(i);
		lookUpGroup[i] = i;
		groupInfo.push_back(tempVec);
	}
	printf("Finish Initialization on Group");

	fes = 0;
}

CCVIL::~CCVIL(){
	delete bestCand;
	delete[] p;
	delete[] lookUpGroup;
}

void CCVIL::run(){
	learningStage();
	optimizationStage();
}


/* 
 * procedure of learning stage, update the "groupInfo"
 */
void CCVIL::learningStage(){
	cout<<"Learning Stage ... "<<endl;
	unsigned cycle = 1; 
 	int lastCycleIndex = -1;
	bool learnStageFlag, needCapture, isSameGroup, separableFunc = true; // assume every benchmark function is separable at the first beginning

	bestCand = new IndividualT<double>(ChromosomeT<double>(param->dimension));
	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());

	popGenerate(true);

	while ( (learnStageFlag = !(cycle > upperThreshold 
					|| (cycle > lowerThreshold && groupInfo.size() == param->dimension) 
					|| groupInfo.size() == 1))==true ){
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
				lastCycleIndex = -1;
			}

			needCapture = groupInfo.size()!=1 && ((cycle<= lowerThreshold) ||(separableFunc == false && cycle <= upperThreshold)) && lastCycleIndex!=-1;
			//	printf("Need Capture ? %d\n", needCapture);
			
			if (lastCycleIndex!=-1){
				//	to decide whether current dimesion are in the same group with last dimension
				isSameGroup = sameGroup(p[i], lastCycleIndex);
			}

			if (lastCycleIndex == -1 || isSameGroup == false){
				printf("i = %d, p[i] = %d\n", i, p[i]);
				// if current dimension and last optimized index are in the different group, then optimize on current dimension
				JADECC(p[i], true /*learnStage = true*/);
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
			if ( lastCycleIndex==-1 || isSameGroup == false ){
				lastCycleIndex = p[i];
			}
		}
		cycle++;
		printf("F %d, Learning Cycle =%d, GroupAmount = %d, BestVal = %.6e \n",
			fp->getID(), 	cycle, 			groupInfo.size(), bestCand->fitnessValue());
	}
}

/*
 * procedure of optimization stage, based on the groupInfo to group the entire population
 */
void CCVIL::optimizationStage(){
	unsigned cycle = 0, groupAmount = groupInfo.size();
	bool learnStageFlag = false;
	printf("Optimization Stage...\n");

	//	learnStageFlag = false, generate new population with regard to the groupInfo
	popGenerate(learnStageFlag);

	while (fes<MaxFitEval){
		printf("Cycle = %d\n", ++cycle);

	
		for (unsigned i=0; i<groupAmount; i++ ) {

			if ( cycle==0 ) {
				//	initialization for the first cycle of optimization stage		
				groupCR = new double[groupAmount];
				groupF = new double[groupAmount];
			}
			//left for restart strategy in optimization stage
			else {
				// <-ELSE PART->
			}

			JADECC(i,learnStageFlag);

		}

		printf("F %d, Optimization Cycle =%d, GroupAmount = %d, BestVal = %.6e \n",
			fp->getID(), 	cycle, 			groupInfo.size(), bestCand->fitnessValue());

	}

	delete[] groupCR;
	delete[] groupF;
}

/* 
 * JADECC: the internal optimizer of CCVIL
 * Optimize on one specific dimension for each run
 *
 * index: the ID of group
 */
void CCVIL::JADECC(unsigned index, bool learnStageFlag){
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
	/***************************** Parameters Setting **************************/
	unsigned D = pop[index][0].size(), NP= pop[index].size(), G, g, *r1, *r2;
	int LB = fp->getMinX(), UB = fp->getMaxX();
	double Fm, CRm, c = param->c, p = param->p;
	double *F, *CR;
	vector<double> goodCR, goodF;
	vector<unsigned> vecIndex;
	Archive* archive = NULL;

	if (learnStageFlag == true){
		G = 1;
	} else if (pop.size() == 1){
		G = INT_MAX;
	} else{
		G = min(D+5, (unsigned)500);
	}

	//  build up the vector of index, depending on whether it is learning stage or not
	if (learnStageFlag == true){
		vecIndex.push_back(index);
	}else{
		vecIndex = groupInfo[index];
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

	if (param->Afactor > 0) {
		// define and initialize the archive
		archive = new Archive(NP * param->Afactor, param->dimension);
		cout<<"create archive"<<endl;
		r1 = new unsigned[NP];
		r2 = new unsigned[NP+archive->getNP()];
	} else {
		r1 = new unsigned[NP];
		r2 = new unsigned[NP];
	}

	/***************************** Population Initialization and evaluation **************************/
	PopulationT<double> parents(NP, ChromosomeT<double>(param->dimension));
	PopulationT<double> offsprings(NP, ChromosomeT<double>(param->dimension));
	PopulationT<double> popAll(0, ChromosomeT<double>(param->dimension));
	parents.setMinimize();
	offsprings.setMinimize();

	printf("The whole population\n");
	print2Dvector(pop);

	printf("Best Candidate at the beginning of each phase\n");
	printPopulation((*bestCand));

	for(unsigned i=0; i < parents.size(); i++){
		// run on individuals
		if (learnStageFlag == true)  {
			// For learning stage, groupInfo may merge internally, while pop always maintain one-dimensional group
			// jth dimension's group belongs to the same group as current optimized group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// run on dimensions
			//	printf("i = %d, j = %d, index = %d\n", i, j, index);
				if (j == index)  {
					parents[i][0][j] = (pop[index])[i][0][0];
				}else if (learnStageFlag == false && lookUpGroup[j] == index){
					//TODO: For optimization stage,
				} else {
					parents[i][0][j] = (*bestCand)[0][j];
				}
			}
		} else{
			// For optimization stage, use goupInfo to enumerate all elements in the same group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// not belongs to the current optimized group
				if (lookUpGroup[j] != index){
					parents[i][0][j] = (*bestCand)[0][j];
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
	}
	printFitness(parents);

	unsigned bestIndex = parents.bestIndex();
	printf("Best Index = %d, Value =  %f\n", parents.bestIndex(), parents.best().fitnessValue());
	
	g = 1;
	fes = fes + NP;

	/***************************** Iterations **************************/
	while ( g<=G && fes < MaxFitEval){
		
		printf("Generation %d\n", g);

		offsprings = parents;

		/*
		   CRm = (1-c)*CRm+c*sum(goodCR)/goodCR.size();
		   Fm = (1-c)*Fm + c*sum(goodF)/goodF.size();     // Lehmer mean
		   */

		// Generate CR according to a normal distribution with mean CRm, and std 0.1
		// Generate F according to a cauchy distribution with location parameter Fm & scale parameter 0.1
		randFCR(NP, CRm, 0.1, Fm, 0.1, F, CR);

		popAll.insert(popAll.size(),offsprings);
		if (param->Afactor == 0){
			// without archive (it is actually a special case of JADE with archive)
			gnR1R2(NP, NP, r1, r2);
		} else {
			// with archive
			// popAll = combinePopulation(offsprings, *(archive->getPop()));
			//popAll.append(*(archive->getPop()));
			popAll.insert(popAll.size(),*(archive->getPop()));
			gnR1R2(NP, NP + archive->getNP(), r1, r2);
		}

		printf("offsprings\n");
		printPopulation(offsprings);

		printf("popAll size = %d\n", popAll.size());
		printf("popAll size = %d, dimension = %d\n", popAll.size(), popAll[0][0].size());
		printPopulation(popAll);

		// Find the p-best solutions
		cout<<"Find the p-best solutions\n"<<endl;
		unsigned pNP = max(round(p*NP), 2.0);
		unsigned* randIndex = new unsigned[NP];
		unsigned* indBest = new unsigned[pNP];
		PopulationT<double> pBestIndiv(NP, ChromosomeT<double>(param->dimension));
		findPbestIndex(offsprings, pNP, indBest);
		for (unsigned i=0; i<NP; i++) {
			randIndex[i] = floor(Rng::uni()*pNP);
			pBestIndiv[i] = offsprings[indBest[randIndex[i]]];
		}

		/*
		printf("randIndex \n");
		printArray(randIndex, NP);

		printf("indBest \n");
		printArray(indBest, pNP);

		printf("pBestIndiv \n");
		printPopulation(pBestIndiv);
		*/

		printf("Begin Mutation\n");
		//************************************ Mutation ************************************//
		PopulationT<double> vi(NP, ChromosomeT<double>(param->dimension));
		vi = offsprings;
		printf("length of vi = %d, offsprings = %d, pBsetIndiv = %d, popAll = %d\n", vi.size(), offsprings.size(), pBestIndiv.size(), popAll.size());
		for (unsigned i=0; i<NP; i++){
			// for each individual
			printf("r1[i] = %d, r2[i] = %d\n", r1[i], r2[i]);
			for (unsigned j=0; j<vecIndex.size(); j++){
				printf("vector index = %d\n", vecIndex[j]);
				// for each dimension to be optimized
				vi[i][0][vecIndex[j]] = offsprings[i][0][vecIndex[j]] + F[i] * ( pBestIndiv[i][0][vecIndex[j]] - offsprings[i][0][vecIndex[j]] + offsprings[r1[i]][0][vecIndex[j]] - popAll[r2[i]][0][vecIndex[j]] );
			}
		}

		printf("Offsprings population\n");
		printPopulation(offsprings);

		printf("vi before bound constrain\n");
		printPopulation(vi);
		boundConstrain(vi, offsprings, LB, UB, vecIndex);

		printf("vi after bound constrain\n");
		printPopulation(vi);

		//************************************ Crossover ************************************//
		printf("Crossover\n");
		PopulationT<double> ui(NP, ChromosomeT<double>(param->dimension));
		for (unsigned i=0; i<NP; i++){
			for (unsigned j=0; j<param->dimension; j++){
				if (j == floor(Rng::uni()*param->dimension) || Rng::uni() < CR[j]){
					ui[i][0][j] = vi[i][0][j];
				}else{
					ui[i][0][j] = offsprings[i][0][j];
				}
			}
		}

		/*
		For j = 1 to D
			If j = jrand or rand(0, 1)< CRi
				u j,i,g = v j,i,g
			Else
				u j,i,g = x j,i,g
			End If
		End For
		*/

		//************************************ Selection ************************************//
		printf("Selection\n");
		for (unsigned i=0; i<NP; i++){
			ui[i].setFitness(fp->compute(ui[i][0]));
			fes += NP;
			if (offsprings[i].fitnessValue() > ui[i].fitnessValue()){
				// improved mutation is saved to offspring, archive the failed solution
				archive->addToArchive(offsprings[i]);

				offsprings[i] = ui[i];

				// save CR & F
				goodCR.push_back(CR[i]);
				goodF.push_back(F[i]);
			}
		}

		printf("Update Archive\n");
		archive->removeDuplicateElem();
		printf("Update Archive\n");
		archive->truncateArchive();
		printf("Update Archive\n");

		// update CRm and Fm
		printf("update CRm and Fm, goodCR size = %d, goodF size = %d\n", (int)goodCR.size(), (int)goodF.size());
		CRm = (1-c)*CRm + c*sum(goodCR)/goodCR.size();
		Fm = (1-c)*Fm + c*sum(dotMultiply(goodF, goodF))/sum(goodF);  // Lehmer Mean 

		printf("Population after Mutation\n");
		printFitness(ui);

		/*
		   If f (xi,g ) ≤ f (ui,g )
		   xi,g+1 = xi,g
		   Else
		   xi,g+1 = ui,g ; xi,g →A; CRi → SCR, Fi → SF
		   End If
		   */

		delete[] randIndex;
		delete[] indBest;
		g++;
	}// iterations

	/***************************** After Iterations Process **************************/
	// update the optimized optimized dimensions to bestCand and pop simultaneously
	if (learnStageFlag == true) {
		for (unsigned j=0; j<offsprings[0][0].size(); j++){
			// update bestCand
			if (j == index) {
				(*bestCand)[0][j] = offsprings[bestIndex][0][j];
			}
			// update pop, the entire external population
			for(unsigned i=0; i<NP; i++){
				(pop[index])[i][0][0] = offsprings[i][0][index];
			}
		}
	} else {
		for (unsigned j=0; j< groupInfo[index].size(); j++){
			unsigned I = groupInfo[index][j];
			(*bestCand)[0][I] = offsprings[bestIndex][0][I];
			for (unsigned i=0; i<NP; i++){
				(pop[index])[i][0][j] = offsprings[i][0][I];
			}
		}
	}

	printf("Best Candidate after update\n");
	printPopulation(*bestCand);
	bestCand->setFitness(offsprings.best().getFitness());

	delete[] F;
	delete[] CR;
	delete[] r1;
	delete[] r2;

	if (param->Afactor > 0){
		delete archive;
	}
}

void CCVIL::printFitness(PopulationT<double> printPop){
	for (unsigned i=0; i < printPop.size(); i++){
		printf("Individual %d, Fitness = %f\n", i, printPop[i].fitnessValue());
	}
}

// if the boundary constraint is violated, set the value to be the middle of
// the previous value and the bound
void CCVIL::boundConstrain(PopulationT<double> &vi, PopulationT<double> offsprings, int LB, int UB, vector<unsigned> vecIndex){
	for (unsigned i=0; i<vi.size(); i++){
		for (unsigned j=0; j<vecIndex.size(); j++){
			printf("vi = %f, LB = %f, UB = %f\n", vi[i][0][vecIndex[j]], (double)LB, (double)UB);
			if ( vi[i][0][vecIndex[j]] < LB){
				printf("Smaller than lower bound\n");
				vi[i][0][vecIndex[j]] = (offsprings[i][0][vecIndex[j]] + (double)LB)/2;
			} else if (vi[i][0][vecIndex[j]]> UB){
				printf("Larger than upper bound\n");
				vi[i][0][vecIndex[j]] = (offsprings[i][0][vecIndex[j]] + (double)UB)/2;
			}
		}
	}
}

// Find and return the indexes of the P% best individuals
void CCVIL::findPbestIndex(PopulationT<double> inPop, unsigned pNP, unsigned* indBest){
	IndividualT<double> temp;
	unsigned iMin, i, tempI, *index = new unsigned[inPop.size()];

	for (i=0; i<inPop.size(); i++){
		index[i] = i;
	}

	for (unsigned iPos=0; iPos<pNP; iPos++){
		iMin = iPos;
		for (i=iPos+1; i<inPop.size(); i++){
			if (inPop[i].getFitness()<inPop[iMin].getFitness()){
				iMin = i;
			}
		}

		if (iMin != iPos){
			// swap individuals in population
			temp = inPop[iPos];
			inPop[iPos] = inPop[iMin];
			inPop[iMin] = temp;

			// swap index
			tempI = index[iPos];
			index[iPos] = index[iMin];
			index[iMin] = tempI;
		}
	}

	for (i=0; i<pNP; i++){
		indBest[i] = index[i];
	}

	if (indBest[0]==indBest[1]){
		cerr<<"ERROR: p best index should be distinct"<<endl;
		exit(-1);
	}

	delete[] index;
}

// combine two populations for mutation process
PopulationT<double> CCVIL::combinePopulation(PopulationT<double> p1, PopulationT<double> p2){
	unsigned D1 = p1[0][0].size(), NP1 = p1.size();
	unsigned NP2 = p2.size();

	// one-on-one copying method
	PopulationT<double> popAll(NP1+NP2, ChromosomeT<double>(D1));
	for (unsigned i=0; i<popAll.size(); i++) {
		for (unsigned j=0; j<D1; j++){
			if (i < NP1){
				popAll[i][0][j] = p1[i][0][j];
			}
			else{
				popAll[i][0][j] = p2[i-NP1][0][j];
			}
		}
	}

	/*
	// try using internal interface "append"
	PopulationT<double> popAll(0, ChromosomeT<double>(D1));
	for (unsigned i=0; i<p1.size() + p2.size(); i++){
		cout<<"size of popAll = "<<popAll.size()<<endl;
		if (i < NP1){
			popAll.append(p1[i]);
		} else {
			popAll.append(p2[i-NP1]);
		}
	}
	*/
	return popAll;
}

void CCVIL::gnR1R2(unsigned NP1, unsigned NP2, unsigned *r1, unsigned *r2){
	unsigned* r0 = new unsigned[NP1];
	bool r2ReGenerateFlag;

	printf("NP1 = %d, NP2 = %d\n", NP1, NP2);
	
	// initialize r0 array
	for (unsigned i=0; i<NP1; i++){
		r0[i] = i;
	}

	for (unsigned i=0; i<NP1; i++){
		r1[i] = floor(Rng::uni()*NP1);
		for (unsigned j=0; j<INT_MAX; j++){
			if (j > 1000){
				cerr<<"Can not genrate r2 in 1000 iterations"<<endl;
				exit(-1);
			}

			if (r1[i]!=r0[i]){
				// distinct
				break;
			}else{
			//	printf("Re-generate NP1\n");
				r1[i] = floor(Rng::uni()*NP1);
			}
		}
		printf("r1[%d] = %d\t", i, r1[i]);
	}
	// eliminate the duplication
	cout<<endl;

	for (unsigned i=0; i<NP2; i++){
		r2[i] = floor(Rng::uni()*NP2);
		for (unsigned j=0; j<INT_MAX ; j++){
			r2ReGenerateFlag = (r2[i]==r0[i] || r2[i]==r1[i]);
			if (j > 1000){
				cerr<<"Can not genrate r1 in 1000 iterations"<<endl;
				exit(-1);
			} 

		//	printf("Begin checking ending criterion, r2ReGenerateFlag = %d\n", r2ReGenerateFlag);
			if (r2ReGenerateFlag == false){
				// distinct
				break;
			}

			if (r2ReGenerateFlag == true){
				r2[i] = floor(Rng::uni()*NP2);
			}
		}
		printf("r2[%d] = %d\t", i, r2[i]);
	}
	cout<<endl;
	delete[] r0;
}

double CCVIL::sum(vector<double> vec){
	double s = 0;
	for (unsigned i=0; i<vec.size(); i++){
		s += vec[i];
	}
	return s;
}

vector<double> CCVIL::dotMultiply(vector<double> v1, vector<double> v2){
	vector<double> vec;
	if (v1.size()!=v2.size()){
		printf("ERROR: Dimensions of vector are not match!");
		exit(-1);
	}else{
		for(unsigned i=0; i < v1.size(); i++){
			vec.push_back(v1[i] * v2[i]);
		}
	}
	return vec;
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

void CCVIL::printPopulation(PopulationT<double> printPop){
	if (printPop.size()>0){
		printf("Dimension of printed Population = %d\n", (printPop)[0][0].size());
		for (unsigned j=0; j<printPop.size(); j++){
			for (unsigned k=0; k<printPop[j][0].size(); k++){
				printf("%f\t",(printPop)[j][0][k]);
			}
			printf("\n");
		}
		cout<<endl;
	}
	else
		printf("size of population is zero!");
}

void CCVIL::printPopulation(IndividualT<double> printIndiv){
		printf("Dimension of printed individual = %d\n", (printIndiv)[0].size());
		for (unsigned j=0; j<(printIndiv)[0].size(); j++){
				printf("%f\t",(printIndiv)[0][j]);
		}
		cout<<endl;
}

void CCVIL::print2Dvector(vector< PopulationT<double> > vector2D){
	for (unsigned i=0; i<vector2D.size(); i++){
		printf("********************\n");
		for (unsigned j=0; j<(vector2D[i]).size(); j++){
			for (unsigned k=0; k<(vector2D[i])[j][0].size(); k++){
				printf("%f\t",(vector2D[i])[j][0][k]);
			}
			printf("\n");
		}
	}
}

void CCVIL::print2Dvector(vector< vector<unsigned> > vector2D){
	for (unsigned i=0; i<vector2D.size(); i++){
		printf("********************\n");
		for (unsigned k=0; k<(vector2D[i]).size(); k++){
			printf("%d\t",(vector2D[i])[k]);
		}
		printf("\n");
	}
}

void CCVIL::printVector(vector<unsigned> v){
	for (unsigned i=0; i<v.size(); i++){
		printf("%d\t", v[i]);
	}
	cout<<endl;
}

void CCVIL::printVector(vector<double> v){
	for (unsigned i=0; i<v.size(); i++){
		printf("%f\t", v[i]);
	}
	cout<<endl;
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

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  popGenerate
 *  Description:  generate the basic structure of the entire population, the pre-process
 *  for population initialization
 * =====================================================================================
 */
void CCVIL::popGenerate(bool learnStageFlag){
	if (learnStageFlag == true){
		for (unsigned i=0; i<param->dimension; i++){
			PopulationT<double> tempPop(param->NP, ChromosomeT<double>(param->initialGroupSize));
			tempPop.setMinimize();
			pop.push_back(tempPop);
		}
	}else{
		//		population generation in optimization stage
		for (unsigned i=0; i<groupInfo.size(); i++){
			unsigned NP = groupInfo[i].size()+10;
			PopulationT<double> tempPop(NP, ChromosomeT<double>(groupInfo[i].size()));
			tempPop.setMinimize();
			pop.push_back(tempPop);
		}
	}
}

/*
 * Initialize the population
 */
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
	double bestCandValue = bestCand->getFitness();
	unsigned NP = pop[lastDim].size(), randi, counter=0;
	IndividualT<double> randIndiv = (*bestCand);
	printf("Capture Interaction Between %d & %d\n best Candidate Fitness Value = %f\n", curDim, lastDim, bestCandValue);
	printf("The population size of last optimized dimension = %d\n", NP);
	
	while (true) {
		randi = Rng::uni()*NP;
		if ( pop[lastDim][randi][0][0] != (*bestCand)[0][lastDim] ){
			break;
		}else if(counter > 1000){
			printf("ERROR: Fail to generate random index within 1000 trying\n");	
		}else{
			printf("Re-generate random index\n");
		}
	}

	randIndiv[0][lastDim] = pop[lastDim][randi][0][0];
	randIndiv.setFitness(fp->compute(randIndiv[0]));
	fes++;

	// if there is any interaction detected, combine the groupInfo
	if (randIndiv.getFitness()<bestCand->getFitness()){
		printf("Interaction Detected\n");
		unsigned group1, group2;
		vector<unsigned> rmVec;


		printf("============= Before Merge =============\n");
		printf("Look Up Group Table:\n");
		printArray(lookUpGroup, param->dimension);

		printf("groupInfo:\n");
		print2Dvector(groupInfo);

		group1 = lookUpGroup[curDim];
		group2 = lookUpGroup[lastDim];
		printf("group1 = %d, group2 =%d\n", group1, group2);

		// through comparison, always join the latter one into the previous
		if(group1 < group2){
			// join group2 into group1
			rmVec = groupInfo[group2];
			groupInfo.erase(groupInfo.begin() + group2);
			for (unsigned i=0; i<rmVec.size(); i++){
				groupInfo[group1].push_back(rmVec[i]);
			}
			lookUpGroup[lastDim] = group1;

			for (unsigned i=0; i<param->dimension; i++){
				if (lookUpGroup[i]>group2){
					lookUpGroup[i]--;
				}
			}


		}else{// group2 < group1
			// join group1 into group2
			rmVec = groupInfo[group1];
			groupInfo.erase(groupInfo.begin() + group1);
			for (unsigned i=0; i<rmVec.size(); i++){
				groupInfo[group2].push_back(rmVec[i]);
			}
			lookUpGroup[curDim] = group2;

			for (unsigned i=0; i<param->dimension; i++){
				if (lookUpGroup[i]>group1){
					lookUpGroup[i]--;
				}
			}
		}

		printf("============= After Merge =============\n");
		printf("Look Up Group Table:\n");
		printArray(lookUpGroup, param->dimension);
		printf("groupInfo:\n");
		print2Dvector(groupInfo);
	}
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
