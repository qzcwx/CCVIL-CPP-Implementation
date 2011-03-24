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
}

CCVIL::~CCVIL(){
	delete bestCand;
	delete[] p;
	delete[] lookUpGroup;
}

void CCVIL::run(){
	double startT, stopT, runTime;


	//
	//	printf ( "lookUpGroup\n" );
	//	printArray(lookUpGroup, runParam->dimension);

	// check if folders "result" and "trace" exist or not
	mkdir ("result", O_CREAT|S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
	mkdir ("trace", O_CREAT|S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	// in result folder
	string resultStr("result/resF");
	resultStr += itos(fp->getID());
	resultStr += ".txt";
	printf("resultStr = %s", resultStr.c_str());
	printf("\n");
	resultFP = fopen(resultStr.c_str(), "w");

	string timeStr("result/timeF");
	timeStr += itos(fp->getID());
	timeStr += ".txt";
	printf("timeStr = %s", timeStr.c_str());
	printf("\n");
	timeFP = fopen(timeStr.c_str(), "w");


	for (unsigned i=0; i < param->numOfRun; i++){
		printf ( "\n\n\n========================== F %d, Run %d ========================\n\n\n", fp->getID(), i+1 );

		/************************* re-initialize the sampling points *************************/
		groupInfo.clear();
		// initialize the groupInfo
		for (unsigned j = 0; j<param->dimension; j++){
			vector<unsigned> tempVec;
			tempVec.push_back(j);
			lookUpGroup[j] = j;
			groupInfo.push_back(tempVec);
		}
//		printf ( "groupInfo\n" );
//		print2Dvector(groupInfo);

		samplingPoints.clear();
		//	printf ( "Sampling Points\n" );
		for (unsigned j=0 ; j<= param->samplingPoint; j++){
			samplingPoints.push_back(j*param->samplingInterval);
			//		printf("%d\n",samplingPoints.back());
		}

		/************************* in trace folder *************************/
		// store grouping information
		string groupStr("trace/groupF");
		groupStr += itos(fp->getID());
		groupStr += "-";
		groupStr += itos(i+1);
		groupStr += ".txt";
		printf("groupStr = %s", groupStr.c_str());
		printf("\n");
		groupFP = fopen(groupStr.c_str(), "w");


		string fesStr("trace/fesF");
		fesStr += itos(fp->getID());
		fesStr += "-";
		fesStr += itos(i+1);
		fesStr += ".txt";
		printf("fesStr = %s", fesStr.c_str());
		printf("\n");
		fesFP = fopen(fesStr.c_str(), "w");

		string valStr("trace/valF");
		valStr += itos(fp->getID());
		valStr += "-";
		valStr += itos(i+1);
		valStr += ".txt";
		printf("valStr = %s", valStr.c_str());
		printf("\n");
		valFP = fopen(valStr.c_str(), "w");



		/* algorithm runing part: start */
		startT = clock();
		fes = 0;
		bestFit = DBL_MAX;
		learningStage();
		optimizationStage();
		stopT = clock();
		/* algorithm runing part: end */


		runTime = (stopT - startT)/CLOCKS_PER_SEC;
		printf ( "Result = %.8e, Running Time = %fs\n", bestFit, runTime );

		resultRec.push_back(bestFit);
		timeRec.push_back(runTime);

		printf ( "\n\n\n========================================================\n\n\n" );

		for (unsigned i=0; i<groupRec.size(); i++){
			fprintf(groupFP, "%d\n", groupRec[i]);
		}
		fclose(groupFP);
		groupRec.clear();

		for (unsigned i=0; i<fesRec.size(); i++){
			fprintf(fesFP, "%d\n", fesRec[i]);
		}
		fclose(fesFP);
		fesRec.clear();

		for (unsigned i=0; i<valRec.size(); i++){
			fprintf(valFP, "%.8e\n", valRec[i]);
		}
		fclose(valFP);
		valRec.clear();
	}

	// delete all file pointers

	// result
	for (unsigned i=0; i<resultRec.size(); i++){
		fprintf(resultFP, "%.8e\n", resultRec[i]);
	}
	fclose(resultFP);
	resultRec.clear();

	// time
	for (unsigned i=0; i<timeRec.size(); i++){
		fprintf(timeFP, "%.8e\n", timeRec[i]);
	}
	fclose(timeFP);
	timeRec.clear();
}
/* 
 * procedure of learning stage, update the "groupInfo"
 */
void CCVIL::learningStage(){
	cout<<"Learning Stage ... "<<endl;
	cycle = 1; 
	int lastCycleIndex = -1;
	bool learnStageFlag, needCapture, isSameGroup = false, separableFunc = true; // assume every benchmark function is separable at the first beginning

	bestCand = new IndividualT<double>(ChromosomeT<double>(param->dimension));
//	(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());

	//	printf ( "Best Cand\n" );
	//	printPopulation((*bestCand));

	popGenerate(true);

	//	printf ( "Compute Learn Stage Flag = %d\n", !(cycle > upperThreshold || (cycle > lowerThreshold && groupInfo.size() == param->dimension) || groupInfo.size() == 1));

	while ( (learnStageFlag = !(cycle > upperThreshold 
					|| (cycle > lowerThreshold && groupInfo.size() == param->dimension) 
					|| groupInfo.size() == 1))==true ){
		// start a new cycle
		//		printf("===================================================\n===================================================\n");
		printf("Cycle = %d\n", cycle);

		for (unsigned i=0; i<param->dimension; i++) {
			//			printf("=========================================================================\nPhase = %d, Cycle = %d\n",i, cycle);
			if (i == 0){
				// start a new phase for each cycle, each phase checking one dimension, i.e., p[i]
				//	printf("Population Re-initialization & Issue Random Permutation\n");
				popInit();
//				printf ( "beforeInitBestCand\n" );
//				print2Dvector(pop);
//				printPopulation((*bestCand));
				initBestCand();
//				printf ( "afterInitBestCand\n" );
//				print2Dvector(pop);
//				printPopulation((*bestCand));
				p = randPerm(param->dimension);
				//				printf("Random Permutation p:\n");
				//				printArray(p, param->dimension);
				lastCycleIndex = -1;
			}

			needCapture = groupInfo.size()!=1 && ((cycle<= lowerThreshold) ||(separableFunc == false && cycle <= upperThreshold)) && lastCycleIndex!=-1;
			//			printf("Need Capture ? %d: between %d and %d\n", needCapture, p[i], lastCycleIndex);

			if (lastCycleIndex!=-1){
				//	to decide whether current dimesion are in the same group with last dimension
				isSameGroup = sameGroup(p[i], lastCycleIndex);
				//				printf("In the Same group? %d\n", isSameGroup);
			}

			if (lastCycleIndex == -1 || isSameGroup == false){
				//				printf("i = %d, p[i] = %d\n", i, p[i]);
				// if current dimension and last optimized index are in the different group, then optimize on current dimension
				JADECC(p[i], true); 	/*learnStage = true*/
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
		printf("F %d, Learning Cycle =%d, GroupAmount = %d, BestVal = %.8e \n",
				fp->getID(), 	cycle, 			(int)groupInfo.size(), bestCand->fitnessValue());
	}
}

/*
 * procedure of optimization stage, based on the groupInfo to group the entire population
 */
void CCVIL::optimizationStage(){
	unsigned	groupAmount = groupInfo.size();
	bool learnStageFlag = false;
	printf("\nOptimization Stage, groupAmount = %d\n", groupAmount);

	cycle = 0;
	//	learnStageFlag = false, generate new population with regard to the groupInfo
	popGenerate(learnStageFlag);

//	printf ( "Group info\n" );
//	print2Dvector(groupInfo);

	popInit();
	//	popInitZeros();

	groupCR = new double[groupAmount];
	groupF = new double[groupAmount];
	failCounter = new unsigned[groupAmount];

	// initialize the control parameters for each group
	for (unsigned i = 0; i< groupAmount; i++){
		groupF[i] = 0.5;
		groupCR[i] = 0.9;
		failCounter[i] = 0;
	}


	double lastCycleBestVal = 0, improveRate=0;
	unsigned innerImprove;

	while (fes<MaxFitEval){
		//		printf("===================================================\n\n\n\n===================================================\n");
		++cycle;
		printf ( "Cycle %d\n", cycle );

//		for (unsigned i=0; i<pop.size(); i++){
//			printf("%d:\tGroupSize = %d,\tPopSize = %d\n", i, groupInfo[i].size(), pop[i].size());
//		}

//		printf("F %d, Optimization Cycle =%d, GroupAmount = %d, fes = %ld, Improved FES = %f, BestVal = %.8e\n",
//				fp->getID(), 	cycle, 			(int)groupInfo.size(), fes, improveRate, bestCand->fitnessValue());

		for (unsigned i=0; i<groupAmount; i++ ) {
			//			printf ( "Phase = %d\n" , i);
			if (failCounter[i] <= param->failThreshold){
				//only optimize on current group, if no a single improvement in the past "failThreshold" successive cycle
				innerImprove = JADECC(i,learnStageFlag);
				if (innerImprove==0){
					failCounter[i] = 0;
				}else{
					failCounter[i] += innerImprove;
				}
			}
			if (sum(failCounter,groupAmount)>=((param->failThreshold+1)*groupAmount)){
				printf ( "*** Restart as no group can be optimized ***\n" );
				popSizeVary(3.0);
				popInit();
//				(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());
//				printf ( "beforeInitBestCand\n" );
//				print2Dvector(pop);
//				printPopulation((*bestCand));
				initBestCand();
//				printf ( "afterInitBestCand\n" );
//				print2Dvector(pop);
//				printPopulation((*bestCand));

				// initialize the control parameters for each group
				for (unsigned i = 0; i< groupAmount; i++){
					groupF[i] = 0.5;
					groupCR[i] = 0.9;
					failCounter[i] = 0;
				}
			}
		}

		improveRate = abs((lastCycleBestVal-bestCand->fitnessValue())/bestCand->fitnessValue());

		if (bestFit==0){
			printf ( "Reach Optmal\n" );
		}

		printf("F %d, Optimization Cycle =%d, GroupAmount = %d, fes = %ld, Improved FES = %f, BestVal = %.8e\n",
				fp->getID(), 	cycle, 			(int)groupInfo.size(), fes, improveRate, bestCand->fitnessValue());

		if (cycle>1 && improveRate<0.01){
			printf ( "*** Restart as non-improvement ***\n" );
			popSizeVary(3.0);
			popInit();
//			(*bestCand)[0].initialize(fp->getMinX(), fp->getMaxX());	
//			printf ( "beforeInitBestCand\n" );
//				print2Dvector(pop);
//			printPopulation((*bestCand));
			initBestCand();
//			printf ( "afterInitBestCand\n" );
//				print2Dvector(pop);
//			printPopulation((*bestCand));

			// initialize the control parameters for each group
			for (unsigned i = 0; i< groupAmount; i++){
				groupF[i] = 0.5;
				groupCR[i] = 0.9;
				failCounter[i] = 0;
			}
		}

		lastCycleBestVal = bestCand->fitnessValue();
	}

	delete[] groupCR;
	delete[] groupF;
	delete[] failCounter;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  popSizeVary
 *  Description:  Change the structure of global population 'pop' according to the factor
 * =====================================================================================
 */
	void
CCVIL::popSizeVary ( double factor )
{
	vector<unsigned> newPopSize;
	
	//	compute the new size of subpopualtion size and store them in vector
	for (unsigned i=0; i<pop.size(); i++){
		newPopSize.push_back(round(factor*pop[i].size()));
	}

	pop.clear();

	//	population generation in optimization stage
	for (unsigned i=0; i<groupInfo.size(); i++){
		PopulationT<double> tempPop(newPopSize[i], ChromosomeT<double>(groupInfo[i].size()));
		tempPop.setMinimize();
		pop.push_back(tempPop);
	}
}		/* -----  end of function popSizeVary  ----- */

/* 
 * JADECC: the internal optimizer of CCVIL
 * Optimize on one specific dimension for each run
 *
 * index: the ID of group
 */
unsigned CCVIL::JADECC(unsigned index, bool learnStageFlag){
	/*********************************************************** 
	 * LB: 	Lower Bound
	 * UB: 	Upper Bound
	 * D: 	Dimension of current group
	 * NP:	Number of Population
	 * G:	Generation Limit for optimizing current group
	 ************************************************************/

	/***************************** Parameters Setting **************************/
	unsigned D = pop[index][0].size(), NP = pop[index].size(), G, g, *r1, *r2;
	int LB = fp->getMinX(), UB = fp->getMaxX();
	double Fm, CRm, c = param->c, p = param->p, preBestVal;
	double *F, *CR;
	// f_rec stands for the fitness value improvement, which is necessary for the adaptation strategy in rJADE
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

//	printf ( "Indices for optimization in this phase\n" );
//	printVector(vecIndex);

//	printf("LB = %d, UB = %d, D = %d, NP = %d, G = %d, index = %d\n", LB, UB, D, NP, G, index);
	if ( learnStageFlag == true ){
		Fm = 0.5;
		CRm = 0.95;
	} else {
		CRm = groupCR[index];
		Fm = groupF[index];
	}

	F = new double[NP];
	CR = new double[NP];

	if (param->Afactor > 0) {
		// define and initialize the archive
		archive = new Archive( (NP * param->Afactor), param->dimension);
		r1 = new unsigned[NP];
		r2 = new unsigned[NP+ archive->getCapacity()];
	} else {
		r1 = new unsigned[NP];
		r2 = new unsigned[NP];
	}

	/***************************** Population Initialization and evaluation **************************/
	PopulationT<double> parents(NP, ChromosomeT<double>(param->dimension));
	PopulationT<double> offsprings(NP, ChromosomeT<double>(param->dimension));

	parents.setMinimize();
	offsprings.setMinimize();

//	printf("The whole population\n");
//	/* Print the struture of the entire population */
//	printf ( "The amount of subpopulation = %d\n", (int)pop.size() );
	
//	printf("The whole population\n");
//	print2Dvector(pop);
//
//	printf("Best Candidate at the beginning of each phase\n");
//	printPopulation((*bestCand));

	for (unsigned i=0; i < parents.size(); i++){
		// run on individuals
		if (learnStageFlag == true)  {
			// For learning stage, groupInfo may merge internally, while pop always maintain one-dimensional group
			// jth dimension's group belongs to the same group as current optimized group
			for (unsigned j=0; j<parents[i][0].size(); j++){
				// run on dimensions
				//	printf("i = %d, j = %d, index = %d\n", i, j, index);
				if (j == index)  {
					parents[i][0][j] = (pop[index])[i][0][0];
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

	
	for (unsigned i=0; i<parents.size(); i++){
		parents[i].setFitness( fp->compute(parents[i][0]) );
	}
//	printFitness(parents);

	unsigned bestIndex = parents.bestIndex();
//	printf("Best Index = %d, Value =  %f\n", parents.bestIndex(), parents.best().fitnessValue());
	
	preBestVal = parents.best().fitnessValue();
	if (preBestVal < bestFit){
		// improve on bestFit
		bestFit = preBestVal;
	}
	
	g = 1;
	fes = fes + NP;
	sampleInfo(preBestVal);
	

	/***************************** Iterations **************************/
	while ( g<=G && fes < MaxFitEval){
//		printf("***************************** Iterations **************************\n");
//		printf("Generation %d, fes %ld\n", g, fes);
	// vector<double> goodCR, goodF; 
		vector<double> goodCR, goodF, f_rec; 

		offsprings = parents;
		PopulationT<double> popAll(parents);

//		printf("Parents Population\n");
//		printPopulation(parents);
		
		// Generate F according to a cauchy distribution with location parameter Fm & scale parameter 0.1
//		printf ( "\nrandFCR, CRm = %f, Fm = %f\n", CRm, Fm );

		randFCR(NP, CRm, 0.1, Fm, 0.1, F, CR);

//		printf ( "F array:\n" );
//		printArray(F, NP);
//		printf ( "CR array:\n" );
//		printArray(CR, NP);

//		printf ( "\ngnR1R2\n" );
		if (param->Afactor == 0){
			// without archive (it is actually a special case of JADE with archive)
			gnR1R2(NP, NP, r1, r2);
		} else {
			// with archive
			//	for (unsigned i=0; i<archive->getNP(); i++){
			//		popAll.append((*archive->getPop())[i]);
			//	}

			if (archive->getNP()>0){
				popAll.insert(popAll.size(), (*archive->getPop()));
			}
			
			gnR1R2(NP, (int)(NP + archive->getNP()), r1, r2);
		}
//		printf ( "r1\n" );
//		printArray(r1, NP);
//		printf ( "r2\n" );
//		printArray(r2, (int)(NP + archive->getNP()));


//		printf("offsprings\n");
//		printPopulation(offsprings);

//		printf("popAll size = %d, dimension = %d\n", popAll.size(), popAll[0][0].size());
//		printPopulation(popAll);

		// Find the p-best solutions
//		printf("Find the p-best solutions\n");
		unsigned pNP = max(round(p*NP), 2.0);
		unsigned* randIndex = new unsigned[NP];
		unsigned* indBest = new unsigned[pNP];
		PopulationT<double> pBestIndiv(NP, ChromosomeT<double>(param->dimension));
		findPbestIndex(offsprings, pNP, indBest);
		for (unsigned i=0; i<NP; i++) {
			randIndex[i] = floor(Rng::uni()*pNP);
			pBestIndiv[i] = offsprings[indBest[randIndex[i]]];
		}

//		printf ( "Fitness of offsprings\n" );
//		printFitness(offsprings);
//
//		printf("randIndex \n");
//		printArray(randIndex, NP);
//
//		printf("indBest \n");
//		printArray(indBest, pNP);
//
//		printf("pBestIndiv \n");
//		printPopulation(pBestIndiv);


//		printf("\nBegin Mutation\n");
		//************************************ Mutation ************************************//
		PopulationT<double> vi(NP, ChromosomeT<double>(param->dimension));
		vi = offsprings;
		for (unsigned i=0; i<NP; i++){
			// for each individual
			//printf("r1[i] = %d, r2[i] = %d\n", r1[i], r2[i]);
			for (unsigned j=0; j<vecIndex.size(); j++){
				//				printf("vector index = %d\n", vecIndex[j]);
				// for each dimension to be optimized
				vi[i][0][vecIndex[j]] = offsprings[i][0][vecIndex[j]] + F[i] * ( pBestIndiv[i][0][vecIndex[j]] - offsprings[i][0][vecIndex[j]] + offsprings[r1[i]][0][vecIndex[j]] - popAll[r2[i]][0][vecIndex[j]] );
			}
		}

//		printf("vi \n");
//		printPopulation(vi);

		boundConstrain(vi, offsprings, LB, UB, vecIndex);

//		printf("vi, after bound constrain \n");
//		printPopulation(vi);

//		printf("\nCrossover\n");
		//************************************ Crossover ************************************//
		PopulationT<double> ui(offsprings);
		ui.setMinimize();
		for (unsigned i=0; i<NP; i++){
			unsigned randIndexCrossover = floor(Rng::uni()*param->dimension);
			for (unsigned j=0; j<vecIndex.size(); j++){
//				printf ( "NP %d, D %d, Inherit mutation: ", i, j );
				if (vecIndex[j] == randIndexCrossover || Rng::uni() < CR[i]){
//					printf ( "Y\n" );
					ui[i][0][vecIndex[j]] = vi[i][0][vecIndex[j]];
				}else{
//					printf ( "N\n" );
				}
			}
		}

//		printf("\nSelection\n");
		//************************************ Selection ************************************//
		for (unsigned i=0; i<NP; i++){
			ui[i].setFitness(fp->compute(ui[i][0]));
			double fitImprv = parents[i].fitnessValue() - ui[i].fitnessValue();
			if ( fitImprv > 0){
				// save CR & F & f_rec
				goodCR.push_back(CR[i]);
				goodF.push_back(F[i]);
				f_rec.push_back(fitImprv);

				// improved mutation is saved to offspring, archive the failed solution
//				printf ( "NP %d, Improved and Updated\n", i );
				archive->addToArchive(parents[i]);
				parents.replace(i, ui[i]);
			}
		}

		if (parents.best().getFitness() < bestFit){
			// improve on bestFit
			bestFit = parents.best().getFitness();
		}
		fes += NP;
		sampleInfo(parents.best().getFitness());

//		printf("Update Archive...\n");
//		printf("Remove Duplicate Elememt\n");
		archive->removeDuplicateElem();

//		printf("Truncate Archive\n");
		archive->truncateArchive();

		// update CRm and Fm
		//printf("update CRm and Fm, goodCR size = %d, goodF size = %d\n", (int)goodCR.size(), (int)goodF.size());

//		/* Adaptation in original JADE */
////		printf("CRm = %f, Fm = %f\n", CRm, Fm);
//		if (goodCR.size()>0 && sum(goodF)>0){
//			CRm = (1-c)*CRm + c*sum(goodCR)/goodCR.size();
//			Fm = (1-c)*Fm + c*sum(dotMultiply(goodF, goodF))/sum(goodF);
////			printf ( "after adaptation, CRm = %f, Fm = %f\n", CRm, Fm );
//		}

		/* Adaptation in rJADE by [Fei Peng al. et. CEC 2009] */
//		printf("rJADE adaptation: CRm = %f, Fm = %f\n", CRm, Fm);
		if (goodCR.size()>0 && sum(goodF)>0){
			CRm = (1-c)*CRm + c*sum(dotMultiply(f_rec,goodCR))/sum(f_rec);
			Fm = (1-c)*Fm + c*sum(dotMultiply(goodF, goodF))/sum(goodF);
//			printf ( "after adaptation, CRm = %f, Fm = %f\n", CRm, Fm );
		}

		//		printf("Population after Mutation, population size of it = %d \n", ui.size());
		//		printPopulation(ui);
		//
		//		printf("Fitness of ui\n");
		//		printFitness(ui);


		if (g % (20000/vecIndex.size()) == 0){
			if ( learnStageFlag == true ) {
				printf("LearningStage, GroupNum = %d, D = %d, groupIndex = %d, Cycle = %d, G = %d, fes = %ld, Best Fitness = %.8e\n", (int)groupInfo.size(), (int)vecIndex.size(), index, cycle, g, fes, parents.best().fitnessValue() );
			}else{
				printf("OptimizationStage, GroupNum = %d, D = %d, groupIndex = %d, Cycle = %d, G = %d, fes = %ld, Best Fitness = %.8e\n", (int)groupInfo.size(), (int)vecIndex.size(), index, cycle, g, fes, parents.best().fitnessValue());
			}
		}
		delete[] randIndex;
		delete[] indBest;

		bestIndex = parents.bestIndex();

		g++;
	}// iterations

	/***************************** After Iterations Process **************************/


	//	printf("goodCR size = %d, goodF size = %d\n", (int)goodCR.size(), (int)goodF.size());

//	// update the optimized optimized dimensions to bestCand and pop simultaneously
//	if (learnStageFlag == true) {
//		for (unsigned j=0; j<offsprings[0][0].size(); j++){
//			// update bestCand
//			if (j == index) {
//				(*bestCand)[0][j] = offsprings[bestIndex][0][j];
//			}
//			// update pop, the entire external population
//			for(unsigned i=0; i<NP; i++){
//				(pop[index])[i][0][0] = offsprings[i][0][index];
//			}
//		}
//	} else {
//		for (unsigned j=0; j< groupInfo[index].size(); j++){
//			unsigned I = groupInfo[index][j];
//			(*bestCand)[0][I] = offsprings[bestIndex][0][I];
//			for (unsigned i=0; i<NP; i++){
//				(pop[index])[i][0][j] = offsprings[i][0][I];
//			}
//		}
//		groupCR[index] = CRm ;
//		groupF[index] = Fm;
//	}
//
//	//	printf("Best Candidate after update\n");
//	//	printPopulation(*bestCand);
//	bestCand->setFitness(offsprings.best().getFitness());
	

	// update the optimized optimized dimensions to bestCand and pop simultaneously
	if (learnStageFlag == true) {
		for (unsigned j=0; j<parents[0][0].size(); j++){
			// update bestCand
			if (j == index) {
				(*bestCand)[0][j] = parents[bestIndex][0][j];
			}
			// update pop, the entire external population
			for(unsigned i=0; i<NP; i++){
				(pop[index])[i][0][0] = parents[i][0][index];
			}
		}
	} else {
		for (unsigned j=0; j< groupInfo[index].size(); j++){
			unsigned I = groupInfo[index][j];
			(*bestCand)[0][I] = parents[bestIndex][0][I];
			for (unsigned i=0; i<NP; i++){
				(pop[index])[i][0][j] = parents[i][0][I];
			}
		}
		groupCR[index] = CRm ;
		groupF[index] = Fm;
	}

	//	printf("Best Candidate after update\n");
	//	printPopulation(*bestCand);
	bestCand->setFitness(parents.best().getFitness());

	delete[] F;
	delete[] CR;
	delete[] r1;
	delete[] r2;

	if (param->Afactor > 0){
		delete archive;
	}

	if (preBestVal - parents.best().getFitness()>0){
		return 0;
	}else{
		return 1;
	}
}

void CCVIL::printFitness(PopulationT<double> printPop){
	for (unsigned i=0; i < printPop.size(); i++){
		printf("Individual %d, Fitness = %0.8e\n", i, printPop[i].fitnessValue());
	}
}

// if the boundary constraint is violated, set the value to be the middle of
// the previous value and the bound
void CCVIL::boundConstrain(PopulationT<double> &vi, PopulationT<double> offsprings, int LB, int UB, vector<unsigned> vecIndex){
	for (unsigned i=0; i<vi.size(); i++){
		for (unsigned j=0; j<vecIndex.size(); j++){
			//			printf("vi = %f, LB = %f, UB = %f\n", vi[i][0][vecIndex[j]], (double)LB, (double)UB);
			if ( vi[i][0][vecIndex[j]] < LB){
				//				printf("Smaller than lower bound\n");
				vi[i][0][vecIndex[j]] = (offsprings[i][0][vecIndex[j]] + (double)LB)/2;
			} else if (vi[i][0][vecIndex[j]]> UB){
				//				printf("Larger than upper bound\n");
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

//	printf("NP1 = %d, NP2 = %d\n", NP1, NP2);
	
	// initialize r0 array
	for (unsigned i=0; i<NP1; i++){
		r0[i] = i;
	}

	for (unsigned i=0; i<NP1; i++){
		r1[i] = floor(Rng::uni()*NP1);
		for (unsigned j=0; j<INT_MAX; j++){
			if (j > 100){
				cerr<<"Can not genrate r1 in 100 iterations"<<endl;
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
//		printf("r1[%d] = %d\t", i, r1[i]);
	}
	// eliminate the duplication

	for (unsigned i=0; i<NP2; i++){
		r2[i] = floor(Rng::uni()*NP2);
		for (unsigned j=0; j<INT_MAX ; j++){
			r2ReGenerateFlag = (r2[i]==r0[i] || r2[i]==r1[i]);
			if (j > 1000){
				cerr<<"Can not genrate r2 in 1000 iterations"<<endl;
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
//		printf("r2[%d] = %d\t", i, r2[i]);
	}
	delete[] r0;
}

double CCVIL::sum(vector<double> vec){
	double s = 0;
	for (unsigned i=0; i<vec.size(); i++){
		s += vec[i];
	}
	return s;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sum
 *  Description:  
 * =====================================================================================
 */
	unsigned
CCVIL::sum ( unsigned* arr, unsigned N )
{
	unsigned totalSum = 0;
	for (unsigned i=0; i<N; i++){
		totalSum += arr[i];
	}
	return totalSum;
}		/* -----  end of function sum  ----- */

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
//	printf ( "Inside the print POP\n" );
	if (printPop.size()>0){
		printf("Dimension of printed Population = %d\n", (printPop)[0][0].size());
		for (unsigned j=0; j<printPop.size(); j++){
			for (unsigned k=0; k<printPop[j][0].size(); k++){
				printf("%f\t",(printPop)[j][0][k]);
			}
			printf("\n=============================================\n");
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
//			F[pos1[i]] = cauchyRnd(Fm, Fsigma);
//			F[pos1[i]] = min(1.0, F[pos1[i]]); // truncated to [-INF, 1]
			F[pos1[i]] = min(1.0, cauchyRnd(Fm, Fsigma)); // truncated to [-INF, 1]
		}

		for (unsigned i=0; i<pos2.size(); i++){
//			F[pos2[i]] = cauchyRnd(Fm, Fsigma);
//			F[pos2[i]] = min(1.0, F[pos2[i]]); // truncated to [-INF, 1]
			double R = cauchyRnd(Fm, Fsigma);
			F[pos2[i]] = min(1.0, R); // truncated to [-INF, 1]
		}

//		printf ( "size of pos1 = %d, pos2 = %d, before clear\n", size1, size2);
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
//	return (mu + delta * tan(PI * (uniRnd()-0.5)));
	return (mu + delta * tan(PI * (uniRnd())));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  popGenerate
 *  Description:  generate the basic structure of the entire population, the pre-process
 *  for population initialization
 * =====================================================================================
 */
void CCVIL::popGenerate(bool learnStageFlag){

	pop.clear();

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
 * ===  FUNCTION  ======================================================================
 *         Name:  popInitZeros
 *  Description:  Initial all genes with 0
 *   Dependency:  popGenerate
 * =====================================================================================
 */
	void
CCVIL::popInitZeros ()
{
	for (unsigned i=0; i<pop.size(); i++){
		for(unsigned j=0; j<pop[i].size(); j++){
			(pop[i])[j][0].initialize(fp->getMinX(), fp->getMaxX());

			//	TODO: Remove Begin
			for (unsigned k=0; k<pop[i][0].size(); k++){
				(pop[i])[j][0][k] = 0;
			}
			//	TODO: Remove End
		}
	}

//	//	TODO: Remove Begin
//	printf ( "Print the entire population\n" );
//	print2Dvector(pop);

	for (unsigned i = 0; i < (*bestCand)[0].size(); i++){
		(*bestCand)[0][i] = 0;
	}

	bestCand->setFitness(fp->compute((*bestCand)[0]));

	printf ( "Fitness of bestCand = %f\n", bestCand->getFitness());

//	printf ( "Print the best Candidate\n" );
//	printPopulation(*bestCand);
//	//	TODO: Remove End
		
}		/* -----  end of function popInitZeros  ----- */

void CCVIL::captureInter(unsigned curDim, unsigned lastDim){
	unsigned NP = pop[lastDim].size(), randi, counter=0;
	IndividualT<double> randIndiv = (*bestCand);
//	printf("Capture Interaction Between %d & %d\n", curDim, lastDim);
//	printf("The population  size of last optimized dimension = %d\n", NP);

	while (true) {
		randi = Rng::uni()*NP;
		if ( pop[lastDim][randi][0][0] != (*bestCand)[0][lastDim] ){
			break;
		}else if(counter > 1000){
			printf("ERROR: Fail to generate random index within 1000 trying\n");	
		}else{
//			printf("Re-generate random index\n");
		}
	}

	randIndiv[0][lastDim] = pop[lastDim][randi][0][0];
	randIndiv.setFitness(fp->compute(randIndiv[0]));
	fes++;

	// if there is any interaction detected, combine the groupInfo
	if (randIndiv.getFitness()<bestCand->getFitness()){
//		printf("Interaction Detected between %d & %d\n", curDim, lastDim);
		unsigned group1, group2;
//		vector<unsigned> rmVec;

//		printf("============= Before Merge =============\n");
		//		printf("Look Up Group Table:\n");
		//		printArray(lookUpGroup, param->dimension);
		
		group1 = lookUpGroup[curDim];
		group2 = lookUpGroup[lastDim];
//		print2Dvector(groupInfo);
//		printf("groupInfo:\n");
//		printf("group1 = %d, group2 =%d\ncurrent D = %d, last D = %d\n", group1, group2, curDim, lastDim);

		// through comparison, always join the latter one into the previous
		if(group1 < group2){
			// join group2 into group1
//			rmVec = groupInfo[group2];
			for (unsigned i=0; i<groupInfo[group2].size(); i++){
//				groupInfo[group1].push_back(rmVec[i]);
				groupInfo[group1].push_back(groupInfo[group2][i]);
			}
//			groupInfo[group2].clear();

			for (unsigned i = 0; i<groupInfo[group2].size(); i++){
				lookUpGroup[groupInfo[group2][i]] = group1;
			}

			groupInfo.erase(groupInfo.begin() + group2);
			for (unsigned i=0; i<param->dimension; i++){
				if (lookUpGroup[i]>group2){
					lookUpGroup[i]--;
				}
			}
		}else{// group2 < group1
			// join group1 into group2
			// rmVec = groupInfo[group1];
			for (unsigned i=0; i<groupInfo[group1].size(); i++){
				groupInfo[group2].push_back(groupInfo[group1][i]);
				//		groupInfo[group2].push_back(rmVec[i]);
			}
			//	groupInfo[group1].clear();

			for (unsigned i = 0; i<groupInfo[group1].size(); i++){
				lookUpGroup[groupInfo[group1][i]] = group2;
			}

			groupInfo.erase(groupInfo.begin() + group1);

			for (unsigned i=0; i<param->dimension; i++){
				if (lookUpGroup[i]>group1){
					lookUpGroup[i]--;
				}
			}
		}

//		printf("============= After Merge =============\n");
		//		printf("Look Up Group Table:\n");
		//		printArray(lookUpGroup, param->dimension);
//		print2Dvector(groupInfo); 
//		printf("groupInfo:\n");
//		printf("group1 = %d, group2 =%d\ncurrent D = %d, last D = %d\n", lookUpGroup[curDim], lookUpGroup[lastDim], curDim, lastDim);
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


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  itos
 *  Description:  
 * =====================================================================================
 */
	string
CCVIL::itos ( int i )
{
	stringstream s;
	s << i;
	return s.str();
}		/* -----  end of function itos  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sampleInfo
 *  Description:  
 * =====================================================================================
 */
	void
CCVIL::sampleInfo ( double curFit )
{
	if (!samplingPoints.empty()&& fes>=samplingPoints.front()){
		samplingPoints.erase(samplingPoints.begin());
		groupRec.push_back(groupInfo.size());
		fesRec.push_back(fes);
		valRec.push_back(curFit);
	}
}
/* -----  end of function sampleInfo  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initBestCand
 *  Description:  initialize the best candidate according to global population
 * =====================================================================================
 */
	void
CCVIL::initBestCand (  )
{
	unsigned curDim=0, randi = 0;
	// assign bestCand from pop in a group by group fashion
	for (unsigned i=0; i<groupInfo.size(); i++){
		for (unsigned j=0; j<groupInfo[i].size(); j++){
			curDim = groupInfo[i][j];
			randi = Rng::uni()*pop[i].size();
			(*bestCand)[0][curDim] = pop[i][randi][0][j] ;
		}
	}
}		/* -----  end of function initBestCand  ----- */
